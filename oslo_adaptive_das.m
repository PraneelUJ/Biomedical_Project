%% Oslo P4 cyst phantom: DAS before/after adaptive RF denoising
% RF comes from MAT; timing and element geometry come from the matching HDF5.
clear; clc; close all;
root=fileparts(mfilename('fullpath'));
datasetName='points'; % 'points' gives bright targets; use 'cyst' for the cyst phantom
if strcmpi(datasetName,'points')
    dataFile=fullfile(root,'oslo','C_verasonics_phantom','FI_P4_point_scatterers_frame1.h5');
    metaFile=fullfile(root,'oslo','C_verasonics_phantom','FI_P4_point_scatterers_metadata.h5');
    R=h5read(dataFile,'/raw');
    if size(R,1)==64 && size(R,3)==128, R=permute(R,[3 2 1]); end
else
    matFile=fullfile(root,'oslo','oslo_rf_frame.mat');
    metaFile=fullfile(root,'oslo','C_verasonics_phantom','FI_P4_cysts_center_metadata.h5');
    S=load(matFile,'raw'); R=squeeze(S.raw);
end
R=single(R); [Ntx,Ns,Nrx]=size(R);

% Acquisition metadata read from the actual HDF5 file.
fs=double(h5read(metaFile,'/sampling_frequency'));
fc=double(h5read(metaFile,'/center_frequency'));
c=double(h5read(metaFile,'/sound_speed'));
t0=double(h5read(metaFile,'/initial_times')); t0=t0(:);
td=double(h5read(metaFile,'/t0_delays'));
if size(td,1)==64, td=td.'; end
ang=double(h5read(metaFile,'/polar_angles')); ang=ang(:);
orig=double(h5read(metaFile,'/transmit_origins'));
if size(orig,1)==3, orig=orig.'; end
pos=double(h5read(metaFile,'/probe_geometry'));
if size(pos,1)==3, pos=pos.'; end
if size(td,1)==Nrx, td=td.'; end
assert(isequal(size(td),[Ntx Nrx]) && isequal(size(pos,1),Nrx) && numel(t0)==Ntx);
fprintf('RF %d tx x %d samples x %d rx | fs %.5f MHz, fc %.5f MHz\n', ...
    Ntx,Ns,Nrx,fs/1e6,fc/1e6);

%% Sector image grid and DAS geometry
Nr=150; Nx=151; r=linspace(5e-3,115e-3,Nr);
xmax=max(r)*sin(max(abs(ang))); x=linspace(-xmax,xmax,Nx); z=linspace(0,max(r),Nr);
origin=median(orig(:,[1 3]),1);
[X0,Z0]=meshgrid(x,z); X=X0+origin(1); Z=Z0+origin(2);
theta=atan2(X-origin(1),Z-origin(2)); radius=hypot(X-origin(1),Z-origin(2));
valid=radius>=min(r) & radius<=max(r) & theta>=min(ang) & theta<=max(ang);
xp=X(:); zp=Z(:); thetaP=theta(:); validP=valid(:);
elx=pos(:,1).'; elz=pos(:,3).'; rxTime=hypot(xp-elx,zp-elz)/c;
lambda=c/fc; aperture=max(elx)-min(elx); txHalf=max(2.5*lambda/aperture,0.5*pi/180);
G=struct('fs',fs,'c',c,'t0',t0,'td',td,'angle',ang,'origin',orig, ...
    'elx',elx,'elz',elz,'rxTime',rxTime,'x',xp,'z',zp,'theta',thetaP, ...
    'valid',validP,'txHalf',txHalf,'imageSize',size(X),'UP',4);

fprintf('Beamforming original-data reference...\n'); imgClean=sectorDAS(R,G);
%% Controlled, repeatable noise level set in the reconstructed image
rng(42,'twister'); noiseUnit=randn(size(R),'single');
imgUnit=sectorDAS(noiseUnit,G); targetImageSNR=5;
signalLevel=mean(abs(imgClean(validP))); unitNoiseLevel=mean(abs(imgUnit(validP)));
sigma=signalLevel/(10^(targetImageSNR/20)*max(unitNoiseLevel,eps));
Rnoisy=R+sigma*noiseUnit; imgNoisy=sectorDAS(Rnoisy,G);
fprintf('Added-noise target SNR %.1f dB; calibrated %.2f dB\n',targetImageSNR, ...
    20*log10(signalLevel/(sigma*unitNoiseLevel+eps)));
noiseMapDraws=8; % independent draws for a stable spatial noise-power estimate

%% Earlier three-band RF filter-bank configuration
% Restore the setup that previously gave a small standalone improvement.
% The intervals describe the overlapping bands; the Gaussian masks and
% passband taper are defined in adaptiveBank below.
edges=fc*[0.42 0.74; 0.70 1.08; 1.02 1.62];
win=max(9,2*floor(1.5*fs/fc)+1);
WienerStrength=0.25; % earlier, gentler RF-domain shrinkage
fprintf('Adaptive filter bank (%d overlapping bands)...\n',size(edges,1));
Rpass=adaptiveBank(R,noiseUnit,fs,fc,sigma^2,win,false,WienerStrength,edges);
RnoisyPass=adaptiveBank(Rnoisy,noiseUnit,fs,fc,0,win,false,WienerStrength,edges);
Rden=adaptiveBank(Rnoisy,noiseUnit,fs,fc,sigma^2,win,true,WienerStrength,edges);
fprintf('Beamforming adaptive RF output...\n'); imgDenoised=sectorDAS(Rden,G);
imgPass=sectorDAS(Rpass,G); imgNoisyPass=sectorDAS(RnoisyPass,G);
snrIn=10*log10(sum(abs(imgPass(validP)).^2)/(sum(abs(imgNoisyPass(validP)-imgPass(validP)).^2)+eps));
snrOut=10*log10(sum(abs(imgPass(validP)).^2)/(sum(abs(imgDenoised(validP)-imgPass(validP)).^2)+eps));
fprintf('Matched-band error SNR: noisy %.2f dB | adaptive %.2f dB\n',snrIn,snrOut);
if snrOut<=snrIn
    warning('Adaptive RF filter did not improve matched-band error SNR; use the band-limited DAS baseline.');
end

%% Post-DAS noise-calibrated complex Wiener shrinkage candidate
% Noise power is measured through the identical RF bandpass and DAS path,
% using independent noise draws so the test noise realization is not reused.
fprintf('Estimating beamformed noise-power map (%d draws)...\n',noiseMapDraws);
noisePowerMap=zeros(size(imgPass));
for k=1:noiseMapDraws
    rng(1000+k,'twister'); nk=randn(size(R),'single');
    nk=adaptiveBank(nk,noiseUnit,fs,fc,0,win,false,WienerStrength,edges);
    ink=sectorDAS(nk,G);
    noisePowerMap=noisePowerMap+abs(ink).^2;
end
noisePowerMap=noisePowerMap/noiseMapDraws*sigma^2;
% Apply exactly the same linear bandpass to clean and noisy RF references.
imgNoisyMatched=sectorDAS(RnoisyPass,G);
imgCleanMatched=imgPass;
localWin=[5 5];
localObs=conv2(abs(imgNoisyMatched).^2,ones(localWin)/prod(localWin),'same');
localNoise=conv2(noisePowerMap,ones(localWin)/prod(localWin),'same');
postGain=max(0,1-localNoise./max(localObs,eps));
postGain(~valid)=0;
imgPostDenoised=imgNoisyMatched.*postGain;
postSNRin=10*log10(sum(abs(imgCleanMatched(validP)).^2)/ ...
    (sum(abs(imgNoisyMatched(validP)-imgCleanMatched(validP)).^2)+eps));
postSNRout=10*log10(sum(abs(imgCleanMatched(validP)).^2)/ ...
    (sum(abs(imgPostDenoised(validP)-imgCleanMatched(validP)).^2)+eps));
fprintf('Post-DAS matched-band error SNR: noisy %.2f dB | Wiener %.2f dB (gain %.2f dB)\n', ...
    postSNRin,postSNRout,postSNRout-postSNRin);

% Identity check: the same gain on the noise-free matched-band image should
% not substantially alter strong point echoes. Report peak amplitude ratio.
imgPostClean=imgCleanMatched.*postGain;
peakMask=valid & abs(imgCleanMatched)>=0.5*max(abs(imgCleanMatched(validP)));
identityPeakRatio=max(abs(imgPostClean(peakMask)))/max(abs(imgCleanMatched(peakMask))+eps);
fprintf('Clean identity check: peak amplitude ratio %.3f (1.0 means unchanged)\n',identityPeakRatio);

% Point-target preservation at the strongest clean target location.
cleanEnvelope=abs(imgCleanMatched); cleanEnvelope(~valid)=0;
[~,targetLin]=max(cleanEnvelope(:)); [targetRow,targetCol]=ind2sub(size(cleanEnvelope),targetLin);
cleanTargetPeak=cleanEnvelope(targetRow,targetCol);
denTargetPeak=abs(imgPostDenoised(targetRow,targetCol));
[cleanAxialFWHM,cleanLateralFWHM]=targetFwhm(cleanEnvelope,targetRow,targetCol,z,x);
[denAxialFWHM,denLateralFWHM]=targetFwhm(abs(imgPostDenoised),targetRow,targetCol,z,x);
fprintf('Strongest target: peak ratio %.3f | axial FWHM %.2f/%.2f mm | lateral FWHM %.2f/%.2f mm (clean/filtered)\n', ...
    denTargetPeak/(cleanTargetPeak+eps),cleanAxialFWHM,denAxialFWHM, ...
    cleanLateralFWHM,denLateralFWHM);

% Independent test noise draws. The Wiener noise map above is not trained on
% these seeds; summarize mean and spread instead of relying on seed 42.
evalSeeds=43:47; seedGainDb=zeros(size(evalSeeds)); seedInputDb=zeros(size(evalSeeds));
for q=1:numel(evalSeeds)
    rng(evalSeeds(q),'twister'); nEval=randn(size(R),'single');
    noisyEvalRF=R+sigma*nEval;
    noisyEvalRF=adaptiveBank(noisyEvalRF,noiseUnit,fs,fc,0,win,false,WienerStrength,edges);
    noisyEval=sectorDAS(noisyEvalRF,G);
    obsEval=conv2(abs(noisyEval).^2,ones(localWin)/prod(localWin),'same');
    gainEval=max(0,1-localNoise./max(obsEval,eps)); gainEval(~valid)=0;
    denEval=noisyEval.*gainEval;
    seedInputDb(q)=10*log10(sum(abs(imgCleanMatched(validP)).^2)/ ...
        (sum(abs(noisyEval(validP)-imgCleanMatched(validP)).^2)+eps));
    seedOutputDb=10*log10(sum(abs(imgCleanMatched(validP)).^2)/ ...
        (sum(abs(denEval(validP)-imgCleanMatched(validP)).^2)+eps));
    seedGainDb(q)=seedOutputDb-seedInputDb(q);
    fprintf('Validation noise seed %d: input %.2f dB, filter gain %+0.2f dB\n', ...
        evalSeeds(q),seedInputDb(q),seedGainDb(q));
end
fprintf('Independent-seed filter gain: mean %+0.2f +/- %.2f dB (%d seeds)\n', ...
    mean(seedGainDb),std(seedGainDb),numel(evalSeeds));

%% Clear three-panel display, same reference and sector coordinates
ref=max(abs(imgClean(validP)));
bClean=20*log10(max(abs(imgClean)/ref,1e-6));
bNoisy=20*log10(max(abs(imgNoisy)/ref,1e-6));
bDenoised=20*log10(max(abs(imgPostDenoised)/ref,1e-6));
bPass=20*log10(max(abs(imgPass)/ref,1e-6));
bClean(~valid)=NaN; bNoisy(~valid)=NaN; bDenoised(~valid)=NaN; bPass(~valid)=NaN;
figure('Color','w','Name','Oslo cyst phantom: DAS denoising','Position',[100 100 1300 520]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
images={bClean,bNoisy,bPass,bDenoised};
labels={'Original-data DAS','DAS with added noise','Band-limited reference','Post-DAS Wiener + DAS'};
for k=1:4
    ax=nexttile; imagesc(ax,x*1e3,z*1e3,images{k}); axis(ax,'image');
    set(ax,'YDir','reverse'); clim(ax,[-55 0]); colormap(ax,gray(256));
    title(ax,labels{k}); xlabel(ax,'Lateral (mm)'); ylabel(ax,'Depth (mm)');
end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='dB relative to clean peak';
save(fullfile(root,'oslo_adaptive_das_results.mat'),'imgClean','imgNoisy', ...
    'imgPass','imgNoisyPass','imgDenoised','imgPostDenoised','noisePowerMap','postGain', ...
    'x','z','valid','fs','fc','c','targetImageSNR','sigma','snrIn','snrOut', ...
    'postSNRin','postSNRout','identityPeakRatio','noiseMapDraws','seedGainDb', ...
    'seedInputDb','evalSeeds','cleanAxialFWHM','denAxialFWHM', ...
    'cleanLateralFWHM','denLateralFWHM','edges','-v7.3');

function out=adaptiveBank(noisy,noiseUnit,fs,fc,noiseVar,win,denoise,strength,edges)
[T,Ns,Nrx]=size(noisy); out=zeros(size(noisy),'like',noisy);
pad=128; nfft=2^nextpow2(Ns+2*pad);
f=(0:nfft-1)'*fs/nfft; f=min(f,fs-f);
nband=3;
edge=fc*[0.42 0.54 1.48 1.62]; pass=zeros(nfft,1);
q=f>=edge(1)&f<edge(2);
pass(q)=0.5-0.5*cos(pi*(f(q)-edge(1))/(edge(2)-edge(1)));
pass(f>=edge(2)&f<=edge(3))=1;
q=f>edge(3)&f<=edge(4);
pass(q)=0.5+0.5*cos(pi*(f(q)-edge(3))/(edge(4)-edge(3)));
centers=fc*[0.62 0.96 1.32]; widths=fc*[0.32 0.34 0.36]; W=zeros(nfft,nband);
for b=1:nband, W(:,b)=exp(-0.5*((f-centers(b))/widths(b)).^2); end
W=W./max(sum(W,2),eps); W=W.*pass;
for tx=1:T
    y=double(squeeze(noisy(tx,:,:))); n=double(squeeze(noiseUnit(tx,:,:)));
    yp=[flipud(y(1:pad,:)); y; flipud(y(end-pad+1:end,:))];
    np=[flipud(n(1:pad,:)); n; flipud(n(end-pad+1:end,:))];
    Y=fft(yp,nfft,1); N=fft(np,nfft,1); crop=pad+(1:Ns);
    row=zeros(Ns,Nrx);
    for b=1:nband
        yb=real(ifft(Y.*W(:,b),[],1)); yb=yb(crop,:);
        if denoise
            nb=real(ifft(N.*W(:,b),[],1)); nb=nb(crop,:);
            % Estimate noise power from the known injected-noise realization.
            pn=noiseVar*mean(nb(:).^2);
            p=movmean(yb.^2,win,1);
            % Restore the earlier partial-Wiener gain used for the ~13 dB run.
            wiener=min(max(1-pn./(p+eps),0),1);
            g=1-strength*(1-wiener); g=movmean(g,9,1);
            yb=yb.*g;
        end
        row=row+yb;
    end
    out(tx,:,:)=single(row);
    if mod(tx,32)==0 || tx==T, fprintf('  filtered transmit %d/%d\n',tx,T); end
end
end

function img=sectorDAS(rf,G)
[T,~,Nrx]=size(rf); P=numel(G.x); accum=complex(zeros(P,1)); wsum=zeros(P,1);
for tx=1:T
    analytic=hilbert(resample(double(squeeze(rf(tx,:,:))),G.UP,1));
    % Minimum element arrival uses the recorded per-element transmit delays.
    txTime=min(G.td(tx,:)+G.rxTime,[],2);
    beamW=max(0,1-abs(G.theta-G.angle(tx))/G.txHalf);
    pix=find(G.valid & beamW>0); if isempty(pix), continue; end
    local=complex(zeros(numel(pix),1)); cnt=zeros(numel(pix),1);
    for rx=1:Nrx
        idx=(txTime(pix)+G.rxTime(pix,rx)-G.t0(tx))*G.fs*G.UP+1;
        ok=idx>=1 & idx<size(analytic,1);
        if ~any(ok), continue; end
        ids=find(ok); i=floor(idx(ok)); f=idx(ok)-i;
        trace=analytic(:,rx);
        local(ids)=local(ids)+(1-f).*trace(i)+f.*trace(i+1); cnt(ids)=cnt(ids)+1;
    end
    q=cnt>0; local(q)=local(q)./cnt(q);
    accum(pix)=accum(pix)+beamW(pix).*local;
    wsum(pix)=wsum(pix)+beamW(pix).*q;
    if mod(tx,32)==0 || tx==T, fprintf('  beamformed transmit %d/%d\n',tx,T); end
end
img=zeros(P,1); q=wsum>0; img(q)=accum(q)./wsum(q);
img=reshape(img,G.imageSize);
end

function [axialMm,lateralMm]=targetFwhm(envelope,row,col,z,x)
% Contiguous -6 dB width through the specified target pixel.
level=envelope(row,col)*10^(-6/20);
axial=abs(envelope(:,col))>=level; lateral=abs(envelope(row,:))>=level;
lo=row; while lo>1 && axial(lo-1), lo=lo-1; end
hi=row; while hi<numel(z) && axial(hi+1), hi=hi+1; end
left=col; while left>1 && lateral(left-1), left=left-1; end
right=col; while right<numel(x) && lateral(right+1), right=right+1; end
axialMm=(z(hi)-z(lo))*1e3;
lateralMm=(x(right)-x(left))*1e3;
end
