%% Reproducible PICMUS denoising method sweep
% Compares post-DAS Wiener settings with RF frequency-bank shrinkage settings.
% This is a separate experiment; it does not modify the original PICMUS script.
% Scores are matched to each method's own clean reference and evaluated over
% independent channel-noise seeds. The sweep is finite and interpretable.
clear; clc; close all;
root=fileparts(mfilename('fullpath'));
file=fullfile(root,'PICMUS_simulation_contrast_speckle_matlab.h5');
if ~isfile(file), error('Missing PICMUS MATLAB HDF5 export: %s',file); end

%% Load channel RF and acquisition geometry
rf=single(h5read(file,'/rf')); % transmit x time x receive
[Ntx,Ns,Nrx]=size(rf);
fs=double(h5read(file,'/sampling_frequency'));
c=double(h5read(file,'/sound_speed'));
t0=double(h5read(file,'/initial_time'));
elx=double(h5read(file,'/element_x')); elz=double(h5read(file,'/element_z'));
angles=double(h5read(file,'/angles')); angles=angles(:);
xAxis=double(h5read(file,'/x_axis')); xAxis=xAxis(:);
zAxis=double(h5read(file,'/z_axis')); zAxis=zAxis(:);
try
    fc=double(h5read(file,'/center_frequency'));
catch
    % The exporter does not currently store fc; use the verified PICMUS value.
    fc=2.71739e6;
end
if isempty(fc) || ~isscalar(fc) || ~isfinite(fc), fc=2.71739e6; end
assert(numel(angles)==Ntx && numel(elx)==Nrx && numel(elz)==Nrx, ...
    'Channel data and acquisition metadata dimensions do not agree.');
fprintf('PICMUS sweep: %d tx x %d samples x %d rx | fs %.4f MHz\n', ...
    Ntx,Ns,Nrx,fs/1e6);

%% Common downsampled image grid and DAS geometry
Nz=171; Nx=151;
iz=unique(round(linspace(1,numel(zAxis),Nz)));
ix=unique(round(linspace(1,numel(xAxis),Nx)));
z=zAxis(iz).'; x=xAxis(ix).';
[X,Z]=meshgrid(x,z); valid=Z>0; validP=valid(:);
xp=X(:); zp=Z(:);
rxTime=hypot(xp-elx(:).',zp-elz(:).')/c;
G=struct('fs',fs,'c',c,'t0',t0,'rxTime',rxTime,'x',xp,'z',zp, ...
    'angles',angles,'valid',validP,'imageSize',size(X),'UP',2);

%% Clean image, noise calibration, and common evaluation seeds
fprintf('Beamforming clean DAS reference...\n');
imgClean=picmusDAS(rf,G);
targetImageSNR=5; rng(42,'twister'); nCal=randn(size(rf),'single');
imgNCal=picmusDAS(nCal,G);
sigma=mean(abs(imgClean(validP)))/(10^(targetImageSNR/20)* ...
    max(mean(abs(imgNCal(validP))),eps));
evalSeeds=42:45; nSeeds=numel(evalSeeds);
noisyImages=cell(nSeeds,1); noisySNR=zeros(nSeeds,1);
noiseRealizations=cell(nSeeds,1);
for s=1:nSeeds
    rng(evalSeeds(s),'twister'); noiseRealizations{s}=randn(size(rf),'single');
    noisyImages{s}=picmusDAS(rf+sigma*noiseRealizations{s},G);
    noisySNR(s)=errorSNR(imgClean,noisyImages{s},validP);
end
fprintf('Image-SNR target %.1f dB; noisy DAS error SNR %.2f +/- %.2f dB (%d seeds)\n', ...
    targetImageSNR,mean(noisySNR),std(noisySNR),nSeeds);

%% Monte-Carlo beamformed noise-power map for post-DAS candidates
noiseMapDraws=4; noisePower=zeros(size(imgClean));
fprintf('Estimating post-DAS noise map (%d independent draws)...\n',noiseMapDraws);
for k=1:noiseMapDraws
    rng(3000+k,'twister'); n=randn(size(rf),'single'); ni=picmusDAS(n,G);
    noisePower=noisePower+abs(ni).^2;
end
noisePower=sigma^2*noisePower/noiseMapDraws;

%% Candidate grid: local complex-image Wiener shrinkage
% Gain = max(0,1-beta*Pn/P)^gamma. Vary neighborhood and shrinkage shape.
windows=[3 3;5 5;7 7]; betas=[0.50 0.75 1.00 1.25]; gammas=[0.50 1.00];
results=struct('name',{},'method',{},'settings',{},'meanSNR',{},'stdSNR',{}, ...
    'meanGain',{},'meanEnvCorr',{},'meanAmpRatio',{},'scores',{}, ...
    'envCorr',{},'ampRatio',{});
caseNo=0; bestWienerDb=-Inf; bestWiener=[]; bestWienerName='';
bestOverallGain=-Inf; bestOverallImage=[]; bestOverallName='';
fprintf('Testing %d post-DAS Wiener settings...\n', ...
    size(windows,1)*numel(betas)*numel(gammas));
for iw=1:size(windows,1)
    win=windows(iw,:); kernel=ones(win)/prod(win);
    pn=conv2(noisePower,kernel,'same');
    for ib=1:numel(betas)
        beta=betas(ib);
        for ig=1:numel(gammas)
            gamma=gammas(ig); scores=zeros(nSeeds,1); ec=scores; ar=scores;
            estimates=cell(nSeeds,1);
            for s=1:nSeeds
                y=noisyImages{s}; p=conv2(abs(y).^2,kernel,'same');
                g=max(0,1-beta*pn./max(p,eps)).^gamma; g(~valid)=0;
                d=y.*g; estimates{s}=d; scores(s)=errorSNR(imgClean,d,validP);
                ec(s)=envCorr(imgClean,d,validP);
                ar(s)=mean(abs(d(validP)))/(mean(abs(imgClean(validP)))+eps);
            end
            caseNo=caseNo+1;
            name=sprintf('POST_WIENER_%dx%d_b%.2f_g%.2f',win(1),win(2),beta,gamma);
            results(caseNo)=packResult(name,'post-DAS Wiener',scores,noisySNR,ec,ar); %#ok<SAGROW>
            if results(caseNo).meanGain>bestOverallGain
                bestOverallGain=results(caseNo).meanGain;
                bestOverallImage=estimates{1}; bestOverallName=name;
            end
            if mean(scores)>bestWienerDb
                bestWienerDb=mean(scores); bestWienerName=name;
                y=noisyImages{1}; p=conv2(abs(y).^2,kernel,'same');
                g=max(0,1-beta*pn./max(p,eps)).^gamma; g(~valid)=0;
                bestWiener=y.*g;
            end
        end
    end
end

%% Candidate grid: overlapping RF filter banks followed by DAS
% A raised-cosine passband limits RF to 0.40-1.60 fc. Gaussian overlapping
% subbands partition this passband. Per-band Wiener gains are blended toward
% identity by shrinkStrength. Each bank is scored against its matched clean
% band-limited DAS image; the no-shrink version is its noisy baseline.
bankCounts=[2 3 4 6 8]; strengths=[0.20 0.40 0.60 0.80];
passbandFc=[0.35 1.65;0.40 1.60;0.50 1.50];
comboWindows=[3 3;5 5]; comboBetas=[0.50 1.00]; comboNoiseDraws=3;
rfWin=max(5,2*floor(fs/fc)+1);
fprintf('Testing RF banks and cascades over %d settings and %d seeds...\n', ...
    size(passbandFc,1)*numel(bankCounts)*(1+numel(strengths)*(1+ ...
    size(comboWindows,1)*numel(comboBetas))),nSeeds);
for ipass=1:size(passbandFc,1)
  for ibank=1:numel(bankCounts)
    nBands=bankCounts(ibank); loFc=passbandFc(ipass,1); hiFc=passbandFc(ipass,2);
    [W,nfft,pad,crop]=makeBank(fs,fc,Ns,nBands,loFc,hiFc);
    cleanPass=filterBank(rf,[],0,W,nfft,pad,crop,rfWin,0);
    imgCleanPass=picmusDAS(cleanPass,G);
    noisyPassSNR=zeros(nSeeds,1); passCorr=zeros(nSeeds,1); passAmp=zeros(nSeeds,1);
    for s=1:nSeeds
        passRF=filterBank(rf+sigma*noiseRealizations{s},[],0, ...
            W,nfft,pad,crop,rfWin,0);
        passImg=picmusDAS(passRF,G);
        noisyPassSNR(s)=errorSNR(imgCleanPass,passImg,validP);
        passCorr(s)=envCorr(imgCleanPass,passImg,validP);
        passAmp(s)=mean(abs(passImg(validP)))/(mean(abs(imgCleanPass(validP)))+eps);
    end
    caseNo=caseNo+1;
    name=sprintf('RF_BANK_%d_BANDPASS_%.2f_%.2f_FC',nBands,loFc,hiFc);
    results(caseNo)=packResult(name,'RF bandpass baseline',noisyPassSNR, ...
        noisyPassSNR,passCorr,passAmp); %#ok<SAGROW>
    for ia=1:numel(strengths)
        strength=strengths(ia); scores=zeros(nSeeds,1); ec=scores; ar=scores;
        denImages=cell(nSeeds,1);
        comboNoise=zeros(size(imgClean));
        fprintf('  %d-band bank, shrink %.2f: calibrating combo noise (%d draws)\n', ...
            nBands,strength,comboNoiseDraws);
        for q=1:comboNoiseDraws
            rng(6000+100*ibank+10*ia+q,'twister'); nn=randn(size(rf),'single');
            noiseRF=filterBank(sigma*nn,nn,sigma^2,W,nfft,pad,crop,rfWin,strength);
            noiseImg=picmusDAS(noiseRF,G);
            comboNoise=comboNoise+abs(noiseImg).^2;
        end
        comboNoise=comboNoise/comboNoiseDraws;
        for s=1:nSeeds
            denRF=filterBank(rf+sigma*noiseRealizations{s},noiseRealizations{s}, ...
                sigma^2,W,nfft,pad,crop,rfWin,strength);
            denImg=picmusDAS(denRF,G);
            denImages{s}=denImg;
            scores(s)=errorSNR(imgCleanPass,denImg,validP);
            ec(s)=envCorr(imgCleanPass,denImg,validP);
            ar(s)=mean(abs(denImg(validP)))/(mean(abs(imgCleanPass(validP)))+eps);
        end
        caseNo=caseNo+1;
        name=sprintf('RF_BANK_%d_BANDS_%.2f_%.2f_FC_S%.2f', ...
            nBands,loFc,hiFc,strength);
        results(caseNo)=packResult(name,'RF filter bank',scores,noisyPassSNR,ec,ar); %#ok<SAGROW>
        if results(caseNo).meanGain>bestOverallGain
            bestOverallGain=results(caseNo).meanGain;
            bestOverallImage=denImages{1}; bestOverallName=name;
        end

        % Test the requested cascade: adaptive RF bank, DAS, then local Wiener.
        for iw=1:size(comboWindows,1)
            win=comboWindows(iw,:); kernel=ones(win)/prod(win);
            pn=conv2(comboNoise,kernel,'same');
            for ib=1:numel(comboBetas)
                beta=comboBetas(ib); comboScores=zeros(nSeeds,1);
                comboCorr=comboScores; comboAmp=comboScores;
                comboEstimates=cell(nSeeds,1);
                for s=1:nSeeds
                    denImg=denImages{s};
                    p=conv2(abs(denImg).^2,kernel,'same');
                    g=max(0,1-beta*pn./max(p,eps)); g(~valid)=0;
                    comboImg=denImg.*g;
                    comboEstimates{s}=comboImg;
                    comboScores(s)=errorSNR(imgCleanPass,comboImg,validP);
                    comboCorr(s)=envCorr(imgCleanPass,comboImg,validP);
                    comboAmp(s)=mean(abs(comboImg(validP)))/(mean(abs(imgCleanPass(validP)))+eps);
                end
                caseNo=caseNo+1;
                comboName=sprintf('RF%d_%.2f_%.2fFC_S%.2f_THEN_WIENER_%dx%d_B%.2f', ...
                    nBands,loFc,hiFc,strength,win(1),win(2),beta);
                results(caseNo)=packResult(comboName,'RF bank + post-DAS Wiener', ...
                    comboScores,noisyPassSNR,comboCorr,comboAmp); %#ok<SAGROW>
                if results(caseNo).meanGain>bestOverallGain
                    bestOverallGain=results(caseNo).meanGain;
                    bestOverallImage=comboEstimates{1}; bestOverallName=comboName;
                end
            end
        end
    end
    fprintf('  completed %d-band bank\n',nBands);
  end
end

%% Rank by mean gain, while reporting stability and image-preservation checks
names={results.name}.'; methods={results.method}.'; settings={results.settings}.';
meanSNR=[results.meanSNR].'; stdSNR=[results.stdSNR].';
meanGain=[results.meanGain].'; meanCorr=[results.meanEnvCorr].';
meanAmp=[results.meanAmpRatio].';
summary=table(methods,names,settings,meanSNR,stdSNR,meanGain,meanCorr,meanAmp, ...
    'VariableNames',{'Method','Configuration','Settings','MeanErrorSNR_dB', ...
    'StdErrorSNR_dB','MeanGain_dB','EnvelopeCorrelation','MeanEnvelopeRatio'});
summary=sortrows(summary,'MeanGain_dB','descend');
fprintf('\nTop configurations (ranked by mean matched-reference error-SNR gain):\n');
disp(summary(1:min(12,height(summary)),:));
fprintf('Best post-DAS Wiener: %s | mean gain %+0.2f dB\n', ...
    bestWienerName,results(find(strcmp(names,bestWienerName),1)).meanGain);
bestWienerGain=results(find(strcmp(names,bestWienerName),1)).meanGain;
fprintf('Best overall method: %s | mean gain %+0.2f dB\n', ...
    bestOverallName,bestOverallGain);

%% Display clean, noisy, and best Wiener example
ref=max(abs(imgClean(validP)));
figure('Color','w','Name','PICMUS filter sweep: best Wiener','Position',[80 120 1500 560]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
images={imgClean,noisyImages{1},bestWiener,bestOverallImage};
labels={'Clean DAS reference',sprintf('Noisy DAS (%.2f dB)',noisySNR(1)), ...
    sprintf('Best post-DAS Wiener (mean gain %+.2f dB)',bestWienerGain), ...
    sprintf('Best overall: %s (%+.2f dB)',bestOverallName,bestOverallGain)};
for k=1:numel(images)
    dB=20*log10(max(abs(images{k})/ref,1e-6)); dB(~valid)=NaN;
    ax=nexttile; imagesc(ax,x*1e3,z*1e3,dB); axis(ax,'image');
    set(ax,'YDir','reverse'); clim(ax,[-55 0]); colormap(ax,gray(256));
    title(ax,labels{k}); xlabel(ax,'Lateral (mm)'); ylabel(ax,'Depth (mm)');
end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='dB relative to clean peak';

outMat=fullfile(root,'PICMUS_filterbank_sweep_results.mat');
outCsv=fullfile(root,'PICMUS_filterbank_sweep_summary.csv');
writetable(summary,outCsv);
save(outMat,'summary','results','noisySNR','evalSeeds','targetImageSNR','sigma', ...
    'imgClean','noisyImages','bestWiener','bestWienerName','bestWienerDb','bestWienerGain', ...
    'bestOverallImage','bestOverallName','bestOverallGain', ...
    'noisePower','x','z','valid','fs','fc','-v7.3');
fprintf('\nSaved full sweep: %s\nSummary table: %s\n',outMat,outCsv);

function r=packResult(name,method,scores,baseline,ec,ar)
r=struct('name',name,'method',method,'settings',name, ...
    'meanSNR',mean(scores),'stdSNR',std(scores), ...
    'meanGain',mean(scores-baseline),'meanEnvCorr',mean(ec), ...
    'meanAmpRatio',mean(ar),'scores',scores,'envCorr',ec,'ampRatio',ar);
end

function [W,nfft,pad,crop]=makeBank(fs,fc,Ns,nBands,loFc,hiFc)
pad=min(128,Ns-1); nfft=2^nextpow2(Ns+2*pad);
f=(0:nfft-1)'*fs/nfft; fa=min(f,fs-f);
lo=loFc*fc; hi=min(hiFc*fc,.98*(fs/2)); tr=.08*fc;
pass=zeros(nfft,1);
q=fa>=lo & fa<lo+tr; pass(q)=.5-.5*cos(pi*(fa(q)-lo)/tr);
pass(fa>=lo+tr & fa<=hi-tr)=1;
q=fa>hi-tr & fa<=hi; pass(q)=.5+.5*cos(pi*(fa(q)-(hi-tr))/tr);
centers=linspace(lo+.5*tr,hi-.5*tr,nBands);
sigma=(hi-lo)/(max(2,nBands)*1.15); W=zeros(nfft,nBands);
for b=1:nBands, W(:,b)=exp(-.5*((fa-centers(b))/sigma).^2); end
W=W./max(sum(W,2),eps); W=W.*pass; crop=pad+(1:Ns);
end

function out=filterBank(rf,noiseUnit,noiseVar,W,nfft,pad,crop,win,strength)
[T,~,Nrx]=size(rf); nBands=size(W,2); out=zeros(size(rf),'single');
for tx=1:T
    y=double(squeeze(rf(tx,:,:))); Y=fft(reflectPad(y,pad),nfft,1);
    row=zeros(numel(crop),Nrx);
    if strength>0
        n=double(squeeze(noiseUnit(tx,:,:))); N=fft(reflectPad(n,pad),nfft,1);
    end
    for b=1:nBands
        yb=real(ifft(Y.*W(:,b),[],1)); yb=yb(crop,:);
        if strength>0
            nb=real(ifft(N.*W(:,b),[],1)); nb=nb(crop,:);
            pn=noiseVar*mean(nb.^2,1);
            p=movmean(yb.^2,win,1);
            wiener=max(0,1-pn./(p+eps));
            gain=1-strength*(1-wiener); gain=movmean(gain,9,1);
            yb=yb.*gain;
        end
        row=row+yb;
    end
    out(tx,:,:)=single(row);
end
end

function y=reflectPad(x,pad)
y=[flipud(x(1:pad,:)); x; flipud(x(end-pad+1:end,:))];
end

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
end
img=zeros(nPix,1); q=validCount>0; img(q)=accum(q)./validCount(q);
img=reshape(img,G.imageSize);
end

function s=errorSNR(reference,estimate,mask)
s=10*log10(sum(abs(reference(mask)).^2)/(sum(abs(estimate(mask)-reference(mask)).^2)+eps));
end

function r=envCorr(reference,estimate,mask)
a=abs(reference(mask)); b=abs(estimate(mask)); a=a-mean(a); b=b-mean(b);
r=sum(a.*b)/(sqrt(sum(a.^2)*sum(b.^2))+eps);
end
