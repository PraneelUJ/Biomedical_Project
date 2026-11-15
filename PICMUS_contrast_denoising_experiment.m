%% PICMUS simulated cyst/speckle channel-data denoising experiment
% Separate from all Oslo scripts. Adds synthetic white RF-channel noise,
% beamforms the same way for all cases, then applies the original
% noise-calibrated local complex Wiener gain after plane-wave compounding.
% Do not apply transmit-coherence weighting here: it reduced matched error SNR.
clear; clc; close all;
root=fileparts(mfilename('fullpath'));
file=fullfile(root,'PICMUS_simulation_contrast_speckle_matlab.h5');
if ~isfile(file)
    error('Run prepare_picmus_contrast.py first.');
end
rf=single(h5read(file,'/rf')); % MATLAB layout: transmit x time x receive
[Ntx,Ns,Nrx]=size(rf);
fs=double(h5read(file,'/sampling_frequency'));
c=double(h5read(file,'/sound_speed'));
t0=double(h5read(file,'/initial_time'));
elx=double(h5read(file,'/element_x')); elz=double(h5read(file,'/element_z'));
angles=double(h5read(file,'/angles')); angles=angles(:);
xAxis=double(h5read(file,'/x_axis')); xAxis=xAxis(:);
zAxis=double(h5read(file,'/z_axis')); zAxis=zAxis(:);
refReal=double(h5read(file,'/reference_real'));
refImag=double(h5read(file,'/reference_imag'));
officialRef=complex(refReal,refImag);
assert(numel(angles)==Ntx && numel(elx)==Nrx && numel(elz)==Nrx, ...
    'Exported channel data and metadata dimensions do not match.');
fprintf('PICMUS contrast/speckle: %d plane waves x %d samples x %d rx | fs %.4f MHz\n', ...
    Ntx,Ns,Nrx,fs/1e6);

%% Downsample the official Cartesian scan grid for a tractable MATLAB run
Nz=171; Nx=151;
iz=unique(round(linspace(1,numel(zAxis),Nz)));
ix=unique(round(linspace(1,numel(xAxis),Nx)));
z=zAxis(iz).'; x=xAxis(ix).';
[X,Z]=meshgrid(x,z); valid=Z>0;
xp=X(:); zp=Z(:); validP=valid(:);
rxTime=hypot(xp-elx(:).',zp-elz(:).')/c;
G=struct('fs',fs,'c',c,'t0',t0,'rxTime',rxTime,'x',xp,'z',zp, ...
    'angles',angles,'valid',validP,'imageSize',size(X),'UP',2);

fprintf('Beamforming clean reference...\n'); imgClean=picmusDAS(rf,G);
officialOnGrid=interp2(xAxis,zAxis,officialRef,X,Z,'linear',0);
beamCorr=abs(sum(conj(officialOnGrid(valid)).*imgClean(valid)))/ ...
    sqrt(sum(abs(officialOnGrid(valid)).^2)*sum(abs(imgClean(valid)).^2)+eps);
fprintf('Custom DAS vs bundled USTB reference normalized correlation: %.4f\n',beamCorr);
refEnv=abs(officialOnGrid(valid)); cleanEnv=abs(imgClean(valid));
refEnv=refEnv-mean(refEnv); cleanEnv=cleanEnv-mean(cleanEnv);
envelopeCorr=sum(refEnv.*cleanEnv)/sqrt(sum(refEnv.^2)*sum(cleanEnv.^2)+eps);
fprintf('Envelope correlation to bundled reference: %.4f\n',envelopeCorr);
%% Fixed 5 dB image-domain additive white noise
targetImageSNR=5; rng(42,'twister'); noiseUnit=randn(size(rf),'single');
imgNoiseUnit=picmusDAS(noiseUnit,G);
sigma=mean(abs(imgClean(validP)))/(10^(targetImageSNR/20)* ...
    max(mean(abs(imgNoiseUnit(validP))),eps));
rfNoisy=rf+sigma*noiseUnit;
imgNoisy=picmusDAS(rfNoisy,G);
fprintf('Added-noise target SNR %.1f dB; calibrated %.2f dB\n',targetImageSNR, ...
    20*log10(mean(abs(imgClean(validP)))/(sigma*mean(abs(imgNoiseUnit(validP)))+eps)));

%% Beamformed noise map and local complex Wiener shrinkage
noiseMapDraws=4; noisePower=zeros(size(imgClean));
fprintf('Estimating DAS noise power (%d independent draws)...\n',noiseMapDraws);
for k=1:noiseMapDraws
    rng(3000+k,'twister'); n=randn(size(rf),'single'); ni=picmusDAS(n,G);
    noisePower=noisePower+abs(ni).^2;
end
noisePower=sigma^2*noisePower/noiseMapDraws;
localWin=[3 3]; kernel=ones(localWin)/prod(localWin);
observedPower=conv2(abs(imgNoisy).^2,kernel,'same');
localNoisePower=conv2(noisePower,kernel,'same');
gain=max(0,1-localNoisePower./max(observedPower,eps)); gain(~valid)=0;
imgDenoised=imgNoisy.*gain;
snrNoisy=errorSNR(imgClean,imgNoisy,validP);
snrDenoised=errorSNR(imgClean,imgDenoised,validP);
fprintf('Matched DAS error SNR: noisy %.2f dB | Wiener %.2f dB (gain %+0.2f dB)\n', ...
    snrNoisy,snrDenoised,snrDenoised-snrNoisy);

%% Repeat on independent noise draws
evalSeeds=43:45; seedGainDb=zeros(size(evalSeeds));
for k=1:numel(evalSeeds)
    rng(evalSeeds(k),'twister'); n=randn(size(rf),'single');
    noisy=picmusDAS(rf+sigma*n,G);
    p=conv2(abs(noisy).^2,kernel,'same'); g=max(0,1-localNoisePower./max(p,eps));
    g(~valid)=0; den=noisy.*g;
    seedGainDb(k)=errorSNR(imgClean,den,validP)-errorSNR(imgClean,noisy,validP);
    fprintf('Noise seed %d: Wiener %+0.2f dB\n',evalSeeds(k),seedGainDb(k));
end
fprintf('Independent-seed gain: %+0.2f +/- %.2f dB (%d seeds)\n', ...
    mean(seedGainDb),std(seedGainDb),numel(seedGainDb));

%% Display and save a dataset-specific result file
ref=max(abs(imgClean(validP)));
figure('Color','w','Name','PICMUS Wiener-only denoising','Position',[80 120 1500 560]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
images={officialRef,imgClean,imgNoisy,imgDenoised};
imageX={xAxis,x,x,x}; imageZ={zAxis,z,z,z};
labels={'Bundled USTB ref','Custom clean DAS','Noisy DAS','Post-DAS Wiener'};
for k=1:4
    localRef=ref; if k==1, localRef=max(abs(officialRef(:))); end
    dB=20*log10(max(abs(images{k})/localRef,1e-6));
    if k>1, dB(~valid)=NaN; end
    ax=nexttile; imagesc(ax,imageX{k}*1e3,imageZ{k}*1e3,dB); axis(ax,'image');
    set(ax,'YDir','reverse'); clim(ax,[-55 0]); colormap(ax,gray(256));
    title(ax,labels{k}); xlabel(ax,'Lateral (mm)'); ylabel(ax,'Depth (mm)');
end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='dB relative to clean peak';
save(fullfile(root,'PICMUS_contrast_denoising_results.mat'), ...
    'officialRef','beamCorr','envelopeCorr','imgClean','imgNoisy','imgDenoised', ...
    'noisePower','gain','x','z','valid','fs','c','targetImageSNR','sigma', ...
    'snrNoisy','snrDenoised','seedGainDb','evalSeeds','-v7.3');

function img=picmusDAS(rf,G)
[Ntx,~,Nrx]=size(rf); nPix=numel(G.x); accum=complex(zeros(nPix,1));
validCount=zeros(nPix,1); up=G.UP;
for tx=1:Ntx
    traces=hilbert(resample(double(squeeze(rf(tx,:,:))),up,1));
    txTime=(G.x*sin(G.angles(tx))+G.z*cos(G.angles(tx)))/G.c;
    idx=(txTime+G.rxTime-G.t0)*G.fs*up+1;
    ok=idx>=1 & idx<size(traces,1);
    i=max(1,min(floor(idx),size(traces,1)-1)); frac=idx-i;
    offsets=(0:Nrx-1)*size(traces,1); lin=i+offsets;
    vals=(1-frac).*traces(lin)+frac.*traces(lin+1); vals(~ok)=0;
    count=sum(ok,2); local=sum(vals,2)./max(count,1);
    use=G.valid & count>0;
    accum(use)=accum(use)+local(use); validCount(use)=validCount(use)+1;
    if mod(tx,15)==0 || tx==Ntx, fprintf('  plane wave %d/%d\n',tx,Ntx); end
end
img=zeros(nPix,1); q=validCount>0; img(q)=accum(q)./validCount(q);
img=reshape(img,G.imageSize);
end

function s=errorSNR(reference,estimate,mask)
s=10*log10(sum(abs(reference(mask)).^2)/(sum(abs(estimate(mask)-reference(mask)).^2)+eps));
end
