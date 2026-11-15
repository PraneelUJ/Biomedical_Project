%% Independent denoising experiment: Oslo L7 CIRS point targets
% Different acquisition/probe from oslo_adaptive_das.m. Noise is synthetic,
% white and channel-independent; it is not a model of actual scanner noise.
clear; clc; close all;
root=fileparts(mfilename('fullpath'));
dataDir=fullfile(root,'oslo','C_verasonics_phantom');
rawFile=fullfile(dataDir,'L7_FI_Verasonics_CIRS_points_frame1.h5');
metaFile=fullfile(dataDir,'L7_FI_Verasonics_CIRS_points_metadata.h5');
if ~isfile(rawFile) || ~isfile(metaFile)
    error('Run prepare_oslo_l7_points.py first to export the Blosc source for MATLAB.');
end
R=single(h5read(rawFile,'/raw')); % tx x time samples x receive element
[Ntx,Ns,Nrx]=size(R);
fs=double(h5read(metaFile,'/sampling_frequency'));
fc=double(h5read(metaFile,'/center_frequency'));
c=double(h5read(metaFile,'/sound_speed'));
t0=double(h5read(metaFile,'/initial_times')); t0=t0(:);
td=double(h5read(metaFile,'/t0_delays'));
ang=double(h5read(metaFile,'/polar_angles')); ang=ang(:);
orig=double(h5read(metaFile,'/transmit_origins'));
if size(orig,1)==3, orig=orig.'; end
pos=double(h5read(metaFile,'/probe_geometry'));
if size(pos,1)==3, pos=pos.'; end
assert(isequal(size(td),[Ntx Nrx]) && size(pos,1)==Nrx && numel(t0)==Ntx, ...
    'RF and exported scan metadata dimensions do not match.');
fprintf('L7 CIRS points: %d tx x %d samples x %d rx | fs %.4f MHz, fc %.4f MHz\n', ...
    Ntx,Ns,Nrx,fs/1e6,fc/1e6);

%% Cartesian image grid and DAS delay geometry
Nz=121; Nx=131; z=linspace(5e-3,100e-3,Nz);
xLimit=max(abs(pos(:,1)))+25e-3; x=linspace(-xLimit,xLimit,Nx);
origin=[median(orig(:,1)) median(orig(:,3))];
[X,Z]=meshgrid(x+origin(1),z+origin(2));
theta=atan2(X-origin(1),Z-origin(2));
valid=Z>origin(2) & theta>=min(ang) & theta<=max(ang);
xp=X(:); zp=Z(:); thetaP=theta(:); validP=valid(:);
elx=pos(:,1).'; elz=pos(:,3).'; rxTime=hypot(xp-elx,zp-elz)/c;
lambda=c/fc; aperture=max(elx)-min(elx); txHalf=max(2.5*lambda/aperture,0.5*pi/180);
G=struct('fs',fs,'c',c,'t0',t0,'td',td,'angle',ang,'elx',elx,'elz',elz, ...
    'rxTime',rxTime,'x',xp,'z',zp,'theta',thetaP,'valid',validP, ...
    'txHalf',txHalf,'imageSize',size(X),'UP',3);

fprintf('DAS clean reference...\n'); imgClean=dasL7(R,G);
%% Repeatable 5 dB image-domain additive-noise experiment
targetImageSNR=5; rng(42,'twister'); noiseUnit=randn(size(R),'single');
imgUnit=dasL7(noiseUnit,G);
signalLevel=mean(abs(imgClean(validP))); noiseLevel=mean(abs(imgUnit(validP)));
sigma=signalLevel/(10^(targetImageSNR/20)*max(noiseLevel,eps));
Rnoisy=R+sigma*noiseUnit; imgNoisy=dasL7(Rnoisy,G);
fprintf('Noise target %.1f dB; measured %.2f dB\n',targetImageSNR, ...
    20*log10(signalLevel/(sigma*noiseLevel+eps)));

%% Three-band passband followed by noise-calibrated post-DAS Wiener
bandEdges=fc*[0.42 0.74; 0.70 1.08; 1.02 1.62];
RcleanBand=bandpassL7(R,fs,fc); RnoisyBand=bandpassL7(Rnoisy,fs,fc);
imgCleanBand=dasL7(RcleanBand,G); imgNoisyBand=dasL7(RnoisyBand,G);
noiseMapDraws=4; noisePower=zeros(size(imgCleanBand));
fprintf('Calibrating band-limited noise map (%d independent draws)...\n',noiseMapDraws);
for k=1:noiseMapDraws
    rng(2000+k,'twister'); n=randn(size(R),'single');
    noiseImage=dasL7(bandpassL7(n,fs,fc),G);
    noisePower=noisePower+abs(noiseImage).^2;
end
noisePower=sigma^2*noisePower/noiseMapDraws;
localWin=[5 5]; kernel=ones(localWin)/prod(localWin);
obsPower=conv2(abs(imgNoisyBand).^2,kernel,'same');
localNoise=conv2(noisePower,kernel,'same');
gain=max(0,1-localNoise./max(obsPower,eps)); gain(~valid)=0;
imgWiener=imgNoisyBand.*gain;
baseDb=errorSnr(imgCleanBand,imgNoisyBand,validP);
outDb=errorSnr(imgCleanBand,imgWiener,validP);
fprintf('Matched-band error SNR: noisy %.2f dB | post-DAS Wiener %.2f dB (gain %+0.2f dB)\n', ...
    baseDb,outDb,outDb-baseDb);

% Point response widths are reported in pixels so sub-grid widths are clear.
env=abs(imgCleanBand); env(~valid)=0;
[~,pk]=max(env(:)); [pr,pc]=ind2sub(size(env),pk);
[cleanAxPix,cleanLatPix]=widthPixels(env,pr,pc);
[denAxPix,denLatPix]=widthPixels(abs(imgWiener),pr,pc);
fprintf('Strongest target: peak ratio %.3f | axial width %d/%d px | lateral width %d/%d px (clean/filtered)\n', ...
    abs(imgWiener(pr,pc))/(env(pr,pc)+eps),cleanAxPix,denAxPix,cleanLatPix,denLatPix);

%% Independent noise-seed check using the same calibrated noise map
evalSeeds=43:45; gainsDb=zeros(size(evalSeeds));
for k=1:numel(evalSeeds)
    rng(evalSeeds(k),'twister'); n=randn(size(R),'single');
    noisyImage=dasL7(bandpassL7(R+sigma*n,fs,fc),G);
    p=conv2(abs(noisyImage).^2,kernel,'same');
    g=max(0,1-localNoise./max(p,eps)); g(~valid)=0;
    filteredImage=noisyImage.*g;
    gainsDb(k)=errorSnr(imgCleanBand,filteredImage,validP)- ...
        errorSnr(imgCleanBand,noisyImage,validP);
    fprintf('Noise seed %d: post-DAS Wiener gain %+0.2f dB\n',evalSeeds(k),gainsDb(k));
end
fprintf('Independent seed gain: %+0.2f +/- %.2f dB (%d seeds)\n', ...
    mean(gainsDb),std(gainsDb),numel(gainsDb));

%% Display with identical dB scale
ref=max(abs(imgClean(validP))); imgs={imgClean,imgNoisy,imgCleanBand,imgWiener};
titles={'Original DAS','Noisy DAS','Band-limited clean reference','Post-DAS Wiener'};
figure('Color','w','Name','L7 CIRS points: denoising experiment','Position',[80 120 1450 520]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
for k=1:4
    q=20*log10(max(abs(imgs{k})/ref,1e-6)); q(~valid)=NaN;
    ax=nexttile; imagesc(ax,x*1e3,z*1e3,q); axis(ax,'image');
    set(ax,'YDir','reverse'); clim(ax,[-55 0]); colormap(ax,gray(256));
    title(ax,titles{k}); xlabel(ax,'Lateral (mm)'); ylabel(ax,'Depth (mm)');
end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='dB relative to clean peak';
save(fullfile(root,'oslo_l7_points_experiment_results.mat'), ...
    'imgClean','imgNoisy','imgCleanBand','imgNoisyBand','imgWiener','noisePower', ...
    'gain','x','z','valid','fs','fc','c','targetImageSNR','sigma', ...
    'baseDb','outDb','gainsDb','evalSeeds','cleanAxPix','denAxPix', ...
    'cleanLatPix','denLatPix','bandEdges','-v7.3');

function rfOut=bandpassL7(rf,fs,fc)
[T,Ns,Nrx]=size(rf); rfOut=zeros(size(rf),'like',rf);
pad=128; nfft=2^nextpow2(Ns+2*pad); f=(0:nfft-1)'*fs/nfft; f=min(f,fs-f);
edge=fc*[0.42 0.54 1.48 1.62]; pass=zeros(nfft,1);
q=f>=edge(1)&f<edge(2); pass(q)=0.5-0.5*cos(pi*(f(q)-edge(1))/(edge(2)-edge(1)));
pass(f>=edge(2)&f<=edge(3))=1;
q=f>edge(3)&f<=edge(4); pass(q)=0.5+0.5*cos(pi*(f(q)-edge(3))/(edge(4)-edge(3)));
centers=fc*[0.62 0.96 1.32]; widths=fc*[0.32 0.34 0.36]; W=zeros(nfft,3);
for b=1:3, W(:,b)=exp(-0.5*((f-centers(b))/widths(b)).^2); end
W=W./max(sum(W,2),eps); W=W.*pass; crop=pad+(1:Ns);
for tx=1:T
    y=double(squeeze(rf(tx,:,:)));
    y=[flipud(y(1:pad,:)); y; flipud(y(end-pad+1:end,:))];
    Y=fft(y,nfft,1); row=zeros(Ns,Nrx);
    for b=1:3
        band=real(ifft(Y.*W(:,b),[],1)); row=row+band(crop,:);
    end
    rfOut(tx,:,:)=single(row);
end
end

function img=dasL7(rf,G)
[T,~,Nrx]=size(rf); P=numel(G.x); accum=complex(zeros(P,1)); wsum=zeros(P,1);
for tx=1:T
    analytic=hilbert(resample(double(squeeze(rf(tx,:,:))),G.UP,1));
    txTime=min(G.td(tx,:)+G.rxTime,[],2);
    beamW=max(0,1-abs(G.theta-G.angle(tx))/G.txHalf);
    pix=find(G.valid & beamW>0); if isempty(pix), continue; end
    local=complex(zeros(numel(pix),1)); cnt=zeros(numel(pix),1);
    for rx=1:Nrx
        idx=(txTime(pix)+G.rxTime(pix,rx)-G.t0(tx))*G.fs*G.UP+1;
        ok=idx>=1 & idx<size(analytic,1); if ~any(ok), continue; end
        ids=find(ok); i=floor(idx(ok)); frac=idx(ok)-i; trace=analytic(:,rx);
        local(ids)=local(ids)+(1-frac).*trace(i)+frac.*trace(i+1); cnt(ids)=cnt(ids)+1;
    end
    q=cnt>0; local(q)=local(q)./cnt(q);
    accum(pix)=accum(pix)+beamW(pix).*local; wsum(pix)=wsum(pix)+beamW(pix).*q;
    if mod(tx,32)==0 || tx==T, fprintf('  DAS transmit %d/%d\n',tx,T); end
end
img=zeros(P,1); q=wsum>0; img(q)=accum(q)./wsum(q); img=reshape(img,G.imageSize);
end

function value=errorSnr(reference,estimate,mask)
value=10*log10(sum(abs(reference(mask)).^2)/(sum(abs(estimate(mask)-reference(mask)).^2)+eps));
end

function [axialPx,lateralPx]=widthPixels(envelope,row,col)
level=envelope(row,col)*10^(-6/20);
a=abs(envelope(:,col))>=level; l=abs(envelope(row,:))>=level;
lo=row; while lo>1 && a(lo-1), lo=lo-1; end
hi=row; while hi<numel(a) && a(hi+1), hi=hi+1; end
left=col; while left>1 && l(left-1), left=left-1; end
right=col; while right<numel(l) && l(right+1), right=right+1; end
axialPx=hi-lo+1; lateralPx=right-left+1;
end
