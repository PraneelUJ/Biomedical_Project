%% Adaptive coherence-factor DAS (CF-DAS)
% Input: fieldII_cyst_RF.mat from the user's Field II simulation.
% Compares conventional DAS and adaptive CF-DAS on identical RF data.
% CF-DAS is adaptive coherence weighting, not MVDR/Capon.
% Requires Signal Processing Toolbox.

clear; clc; close all;
S = load('fieldII_cyst_RF.mat');
need = {'rf_clean','noise_unit','fs','f0','c','time_axis','element_x', ...
    'tx_focus_x','tx_rmax','pulse_lag','z_focus','pitch','lambda', ...
    'cyst_x','cyst_z','cyst_radius','N_elements','N_tx'};
for k=1:numel(need)
    assert(isfield(S,need{k}),'Missing %s in fieldII_cyst_RF.mat',need{k});
end
rf_clean=S.rf_clean; noise_unit=S.noise_unit;
fs=S.fs; f0=S.f0; c=S.c; time_axis=S.time_axis; element_x=S.element_x;
tx_focus_x=S.tx_focus_x; tx_rmax=S.tx_rmax; pulse_lag=S.pulse_lag;
z_focus=S.z_focus; pitch=S.pitch; lambda=S.lambda;
cyst_x=S.cyst_x; cyst_z=S.cyst_z; cyst_radius=S.cyst_radius;
[Ns,Nrx,Ntx]=size(rf_clean);
assert(Nrx==S.N_elements && Ntx==S.N_tx);
assert(isequal(size(noise_unit),size(rf_clean)));

%% Settings
% Inject channel noise before beamforming. The saved noise_unit makes the
% experiment repeatable and allows each method to see exactly the same noise.
image_SNR_dB=15; % target mean speckle/noise level at the cyst-depth background
run_snr_sweep=false; % enable only when you want the longer numeric comparison
image_snr_sweep_dB=[5 10 15 20];
sweep_results=struct([]); % always defined so the results save works when sweep is off
F_number=2;
UP=8;
CF_gamma=1.0;       % try 0.7, 1.0, 1.5 and compare
CF_floor=0.02;      % prevents complete suppression of weak echoes
Wiener_floor=0.08;  % minimum local Wiener gain
Adaptive_floor=0.08;
band_edges=[1.35 2.15; 1.85 3.25; 2.85 4.30]*1e6; % overlapping RF bands
smooth_samples=max(5,2*floor((fs/f0)*1.5)+1); % about 1.5 pulse cycles
DR=45;
x_image=-12e-3:0.10e-3:12e-3;
z_image=30e-3:0.10e-3:100e-3;
[X,Z]=meshgrid(x_image,z_image);
[Nz,Nx]=size(X); Npix=numel(X);
Xv=X(:); Zv=Z(:);

%% Receive delays and apodization
fprintf('Precomputing geometry...\n');
trx=zeros(Npix,Nrx,'single'); rw=zeros(Npix,Nrx,'single');
for r=1:Nrx
    dx=abs(Xv-element_x(r));
    trx(:,r)=single(sqrt(dx.^2+Zv.^2)/c);
    nd=dx./(0.5*Zv/F_number);
    w=ones(Npix,1);
    q=nd>0.5 & nd<=1;
    w(q)=0.5*(1+cos(pi*(nd(q)-0.5)/0.5));
    w(nd>1)=0;
    rw(:,r)=single(w);
end
G.X=X; G.Z=Z; G.fs=fs; G.UP=UP; G.c=c;
G.t_start=time_axis(1); G.lag=pulse_lag; G.z_focus=z_focus;
G.tx_focus_x=tx_focus_x; G.tx_rmax=tx_rmax;
G.t_rx=trx; G.rx_w=rw;
G.a_min=lambda*z_focus/(S.N_elements*pitch);
G.tan_half=(S.N_elements*pitch/2)/z_focus;

%% Reconstruct clean DAS and adaptive CF-DAS
fprintf('\nConventional DAS...\n'); tic;
img_das=beamform(rf_clean,G,false,CF_gamma,CF_floor);
fprintf('DAS time %.1f s\n',toc);
fprintf('\nAdaptive CF-DAS...\n'); tic;
[img_cf,cf_map]=beamform(rf_clean,G,true,CF_gamma,CF_floor);
fprintf('CF-DAS time %.1f s\n',toc);

%% Calibrate channel noise to a meaningful image-domain target level
% Global RF SNR is dominated by shallow echoes and does not control whether
% the cyst-depth image is interpretable. Match the existing DAS_beamforming.m
% approach: set noise from mean speckle/noise envelope in the background ROI.
d=sqrt((X-cyst_x).^2+(Z-cyst_z).^2);
cyst_mask=d<=0.7*cyst_radius;
bg_mask=d>=1.4*cyst_radius & d<=2*cyst_radius;
img_noise_unit=beamform(noise_unit,G,false,CF_gamma,CF_floor);
speckle_level=mean(abs(img_das(bg_mask)));
unit_noise_level=mean(abs(img_noise_unit(bg_mask)));
noise_sigma=speckle_level/(10^(image_SNR_dB/20)*max(unit_noise_level,eps));
noise_power=noise_sigma^2;
rf_noisy=rf_clean+noise_sigma*noise_unit;
fprintf('Target cyst-background image SNR %.1f dB; measured %.2f dB\n', ...
    image_SNR_dB,20*log10(speckle_level/(mean(abs(noise_sigma*img_noise_unit(bg_mask)))+eps)));
N_sig=min(S.N_sig,Ns);
signal_power=mean(double(rf_clean(1:N_sig,:,:)).^2,'all');
fprintf('Equivalent global channel SNR %.2f dB\n', ...
    10*log10(signal_power/noise_power));

fprintf('\nNoisy channel DAS and CF-DAS...\n'); tic;
img_das_n=beamform(rf_noisy,G,false,CF_gamma,CF_floor);
[img_cf_n,cf_noisy]=beamform(rf_noisy,G,true,CF_gamma,CF_floor);
fprintf('Noisy reconstruction time %.1f s\n',toc);

%% Adaptive three-band channel Wiener denoising + coherence weighting
% Noise variance is known from the injected channel noise. Each frequency
% band gets a local Wiener gain; CF from the band-limited beamformed data then
% downweights bands whose delayed channel signals do not agree.
nBands=size(band_edges,1);
img_bank=cell(nBands,1); img_bank_clean=cell(nBands,1);
cf_bank=cell(nBands,1); cf_clean_bank=cell(nBands,1); band_noise_power=zeros(nBands,1);
rf_band_clean=cell(nBands,1); rf_band_noisy=cell(nBands,1);
fprintf('\nAdaptive filter bank (%d bands)...\n',nBands); tic;
for b=1:nBands
    [rf_band_clean{b},rf_band_noisy{b},band_noise_power(b)]= ...
        band_wiener(rf_clean,rf_noisy,noise_unit,fs,band_edges(b,:), ...
        noise_power,smooth_samples,Wiener_floor);
    [img_bank_clean{b},cf_clean]=beamform(rf_band_clean{b},G,true,CF_gamma,CF_floor);
    [img_bank{b},cf_n]=beamform(rf_band_noisy{b},G,true,CF_gamma,CF_floor);
    cf_bank{b}=cf_n;
    cf_clean_bank{b}=cf_clean;
    fprintf('  band %d: %.2f-%.2f MHz, noise variance %.3g\n', ...
        b,band_edges(b,1)/1e6,band_edges(b,2)/1e6,band_noise_power(b));
end
% Local band strength avoids selecting a band just because its CF is high.
bank_stack=cat(3,img_bank{:}); cf_stack=cat(3,cf_bank{:});
band_energy=abs(bank_stack).^2;
energy_total=sum(band_energy,3);
energy_weight=band_energy./(energy_total+eps);
coh_weight=max(cf_stack,Adaptive_floor).^2;
adaptive_weights=energy_weight.*coh_weight;
adaptive_weights=adaptive_weights./(sum(adaptive_weights,3)+eps);
% Soft image-domain gate. Channel-domain noise variance cannot be subtracted
% directly from beamformed pixel energy (the units/scaling differ). Estimate
% a conservative local floor from the lower tail of band image energy, then
% combine that SNR cue with coherence without hard-zeroing the image.
positive_energy=energy_total(isfinite(energy_total) & energy_total>0);
if isempty(positive_energy)
    noise_image_floor=0;
else
    noise_image_floor=0.05*prctile(positive_energy,20);
end
energy_gate=energy_total./(energy_total+noise_image_floor+eps);
energy_gate(~isfinite(energy_gate))=1;
coherence_gate=0.5+0.5*sqrt(max(mean(cf_stack,3),0));
pixel_gate=min(max(energy_gate.*coherence_gate,0.35),1);
img_adaptive=sum(adaptive_weights.*bank_stack,3);
img_adaptive=img_adaptive.*pixel_gate;

% Estimate weights from denoised noisy bands. Freeze those weights for the
% clean counterpart so the comparison measures reconstruction error rather
% than differences in a data-dependent weighting rule.
clean_stack=cat(3,img_bank_clean{:});
clean_power=abs(clean_stack).^2; clean_total=sum(clean_power,3);
clean_cf_stack=cat(3,cf_clean_bank{:});
clean_weights=clean_power./(clean_total+eps).*max(clean_cf_stack,Adaptive_floor).^2;
clean_weights=clean_weights./(sum(clean_weights,3)+eps);
clean_positive=clean_total(isfinite(clean_total)&clean_total>0);
if isempty(clean_positive), clean_floor=0; else, clean_floor=0.05*prctile(clean_positive,20); end
clean_gate=clean_total./(clean_total+clean_floor+eps); clean_gate(~isfinite(clean_gate))=1;
clean_gate=min(max(clean_gate.*(0.5+0.5*sqrt(max(mean(clean_cf_stack,3),0))),0.35),1);
img_adaptive_clean=sum(clean_weights.*clean_stack,3).*clean_gate;
fprintf('Filter bank processing time %.1f s\n',toc);
%% COHERENCE FACTOR DIAGNOSTICS

cf_values = cf_map(isfinite(cf_map));

fprintf('\n========== CF DIAGNOSTICS ==========\n');
fprintf('Minimum CF : %.8f\n', min(cf_values));
fprintf('Mean CF    : %.8f\n', mean(cf_values));
fprintf('Median CF  : %.8f\n', median(cf_values));
fprintf('Maximum CF : %.8f\n', max(cf_values));

fprintf('\nCF percentiles:\n');
fprintf('P10 : %.8f\n', prctile(cf_values,10));
fprintf('P25 : %.8f\n', prctile(cf_values,25));
fprintf('P50 : %.8f\n', prctile(cf_values,50));
fprintf('P75 : %.8f\n', prctile(cf_values,75));
fprintf('P90 : %.8f\n', prctile(cf_values,90));
fprintf('P99 : %.8f\n', prctile(cf_values,99));

fprintf('\nPercentage of image below CF floor: %.2f%%\n', ...
    100*mean(cf_values < CF_floor));
%% Cyst masks and metrics
%% Metrics
e_dc=abs(img_das); e_dn=abs(img_das_n);
e_cc=abs(img_cf);  e_cn=abs(img_cf_n);
e_ac=abs(img_adaptive_clean); e_an=abs(img_adaptive);
% Normalize clean/noisy pairs to a shared clean peak. Do not normalize a
% noisy image to its own peak, which makes noise appear as a white image.
refd=max(e_dc(:)); refc=max(e_cc(:)); refa=max(e_ac(:));
b_dc=20*log10(max(e_dc/refd,1e-8));
b_dn=20*log10(max(e_dn/refd,1e-8));
b_cc=20*log10(max(e_cc/refc,1e-8));
b_cn=20*log10(max(e_cn/refc,1e-8));
b_ac=20*log10(max(e_ac/refa,1e-8));
b_an=20*log10(max(e_an/refa,1e-8));
mdc=metrics(e_dc,cyst_mask,bg_mask);
mdn=metrics(e_dn,cyst_mask,bg_mask);
mcc=metrics(e_cc,cyst_mask,bg_mask);
mcn=metrics(e_cn,cyst_mask,bg_mask);
m_ac=metrics(e_ac,cyst_mask,bg_mask); m_an=metrics(e_an,cyst_mask,bg_mask);
snr_d=region_snr(img_das,img_das_n,bg_mask);
snr_cf=region_snr(img_cf,img_cf_n,bg_mask);
snr_adaptive=region_snr(img_adaptive_clean,img_adaptive,bg_mask);

fprintf('\n====== CYST METRICS, background image SNR %.1f dB ======\n',image_SNR_dB);
fprintf('Method             CNR      CR(dB)    gCNR   region SNR(dB)\n');
fprintf('DAS clean       %7.3f   %8.2f   %6.3f        --\n',mdc.CNR,mdc.CR,mdc.gCNR);
fprintf('DAS noisy       %7.3f   %8.2f   %6.3f      %7.2f\n',mdn.CNR,mdn.CR,mdn.gCNR,snr_d);
fprintf('CF-DAS clean    %7.3f   %8.2f   %6.3f        --\n',mcc.CNR,mcc.CR,mcc.gCNR);
fprintf('CF-DAS noisy    %7.3f   %8.2f   %6.3f      %7.2f\n',mcn.CNR,mcn.CR,mcn.gCNR,snr_cf);
fprintf('Adaptive clean  %7.3f   %8.2f   %6.3f        --\n',m_ac.CNR,m_ac.CR,m_ac.gCNR);
fprintf('Adaptive noisy  %7.3f   %8.2f   %6.3f      %7.2f\n',m_an.CNR,m_an.CR,m_an.gCNR,snr_adaptive);
fprintf('Region SNR is 10log10(clean complex-image energy / squared complex-image error).\n');
fprintf('Mean CF background %.4f | cyst %.4f\n', ...
    mean(cf_noisy(bg_mask),'omitnan'),mean(cf_noisy(cyst_mask),'omitnan'));
fprintf('Adaptive pixel gate: mean %.4f, max %.4f\n',mean(pixel_gate(:)),max(pixel_gate(:)));

%% Optional image-SNR sweep using the same clean RF and noise realization
% The requested SNR is defined in the cyst-depth background ROI.
if run_snr_sweep
    sweep_results=struct('image_SNR_dB',{},'DAS_SNR',{},'Adaptive_SNR',{}, ...
        'DAS_CNR',{},'Adaptive_CNR',{},'DAS_CR',{},'Adaptive_CR',{}, ...
        'DAS_gCNR',{},'Adaptive_gCNR',{});
    fprintf('\n====== CYST-BACKGROUND IMAGE SNR SWEEP ======\n');
    fprintf('Image SNR DAS err-SNR Adapt err-SNR DAS CNR Adapt CNR DAS CR Adapt CR DAS gCNR Adapt gCNR\n');
    for si=1:numel(image_snr_sweep_dB)
        this_snr=image_snr_sweep_dB(si);
        this_sigma=speckle_level/(10^(this_snr/20)*max(unit_noise_level,eps));
        this_noise_power=this_sigma^2;
        this_rf=rf_clean+this_sigma*noise_unit;
        this_das=beamform(this_rf,G,false,CF_gamma,CF_floor);
        [this_adaptive,this_adaptive_clean]=adaptive_reconstruct(this_rf,rf_clean,noise_unit,G,fs, ...
            band_edges,this_noise_power,smooth_samples,Wiener_floor, ...
            CF_gamma,CF_floor,Adaptive_floor);
        md=metrics(abs(this_das),cyst_mask,bg_mask);
        ma=metrics(abs(this_adaptive),cyst_mask,bg_mask);
        sd=region_snr(img_das,this_das,bg_mask);
        sa=region_snr(this_adaptive_clean,this_adaptive,bg_mask);
        sweep_results(si)=struct('image_SNR_dB',this_snr,'DAS_SNR',sd, ...
            'Adaptive_SNR',sa,'DAS_CNR',md.CNR,'Adaptive_CNR',ma.CNR, ...
            'DAS_CR',md.CR,'Adaptive_CR',ma.CR,'DAS_gCNR',md.gCNR, ...
            'Adaptive_gCNR',ma.gCNR);
        fprintf('%5.1f dB %10.2f dB %14.2f dB %9.3f %10.3f %8.2f %9.2f %9.3f %10.3f\n', ...
            this_snr,sd,sa,md.CNR,ma.CNR,md.CR,ma.CR,md.gCNR,ma.gCNR);
    end
end

%% Display only the phantom depth range with one shared clean-reference scale
th=linspace(0,2*pi,300);
cx=(cyst_x+cyst_radius*cos(th))*1e3;
cz=(cyst_z+cyst_radius*sin(th))*1e3;
show_z=z_image>=28e-3 & z_image<=90e-3;
figure('Color','w','Name','Cyst phantom denoising','Position',[100 100 1350 540]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
titles={'Clean reference','Noisy channel data: DAS','Adaptive filter-bank output'};
ims={b_dc(show_z,:),b_dn(show_z,:),b_an(show_z,:)};
for k=1:3
    ax=nexttile;
    imagesc(ax,x_image*1e3,z_image(show_z)*1e3,ims{k});
    axis(ax,'image'); set(ax,'YDir','reverse'); colormap(ax,gray(256)); clim(ax,[-DR 0]);
    hold(ax,'on'); plot(ax,cx,cz,'r','LineWidth',1.5);
    title(ax,titles{k}); xlabel(ax,'Lateral position (mm)'); ylabel(ax,'Depth (mm)');
end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='dB relative to clean-image peak';

save('DAS_adaptive_CF_results.mat','img_das','img_cf','img_das_n','img_cf_n', ...
    'img_adaptive','img_adaptive_clean','cf_map','cf_noisy','adaptive_weights','pixel_gate', ...
    'band_edges','band_noise_power','b_dc','b_cc','b_dn','b_cn','b_ac','b_an', ...
    'mdc','mdn','mcc','mcn','m_ac','m_an','snr_d','snr_cf','snr_adaptive', ...
    'x_image','z_image','cyst_mask','bg_mask','image_SNR_dB','noise_sigma', ...
    'CF_gamma','CF_floor','Wiener_floor','Adaptive_floor','F_number','UP', ...
    'run_snr_sweep','image_snr_sweep_dB','sweep_results','-v7.3');
fprintf('\nSaved DAS_adaptive_CF_results.mat\n');

%% Local functions
function [img,cfimg]=beamform(rf,G,usecf,gamma,cffloor)
[~,M,T]=size(rf); [nz,nx]=size(G.X); fsu=G.fs*G.UP;
x=G.X(:); z=G.Z(:); P=numel(x);
num=complex(zeros(P,1)); den=zeros(P,1);
cfnum=zeros(P,1); cfden=zeros(P,1);
ah=sqrt(G.a_min^2+((z-G.z_focus)*G.tan_half).^2);
for t=1:T
    xf=G.tx_focus_x(t); u=abs(x-xf)./ah;
    pix=find(u<1); if isempty(pix), continue; end
    wt=0.5*(1+cos(pi*u(pix)));
    dz=z(pix)-G.z_focus;
    ttx=(G.tx_rmax(t)+sign(dz).*sqrt((x(pix)-xf).^2+dz.^2))/G.c+G.lag;
    a=hilbert(resample(double(rf(:,:,t)),G.UP,1));
    L=size(a,1); acc=complex(zeros(numel(pix),1));
    wsum=zeros(numel(pix),1); energy=zeros(numel(pix),1);
    for r=1:M
        w=double(G.rx_w(pix,r));
        idx=(ttx+double(G.t_rx(pix,r))-G.t_start)*fsu+1;
        ok=idx>=1 & idx<L & w>0;
        i=floor(idx); i(~ok)=1; f=idx-i;
        col=a(:,r);
        s=(1-f).*col(i)+f.*col(i+1);
        s(~ok)=0; w(~ok)=0;
        acc=acc+w.*s; wsum=wsum+w;
        energy=energy+w.*abs(s).^2;
    end
    good=wsum>0; y=zeros(numel(pix),1);
    y(good)=acc(good)./wsum(good);
    if usecf
        q=wsum.*energy; cf=zeros(size(q));
        ok=q>eps; cf(ok)=abs(acc(ok)).^2./q(ok);
        cf=min(max(real(cf),0),1);
        adapt=max(cf,cffloor).^gamma; adapt(~good)=0;
        effective_weight = wt.*adapt.*good;

        num(pix)=num(pix)+effective_weight.*y;
        den(pix)=den(pix)+effective_weight;
        cfnum(pix)=cfnum(pix)+wt.*cf.*good;
        cfden(pix)=cfden(pix)+wt.*good;
    else
        num(pix)=num(pix)+wt.*y;
        den(pix)=den(pix)+wt.*good;
    end
    if mod(t,8)==0 || t==T, fprintf('  transmit %d/%d\n',t,T); end
end
img=zeros(P,1); q=den>0; img(q)=num(q)./den(q);
img=reshape(img,nz,nx);
if nargout>1
    cfimg=zeros(P,1); q=cfden>0; cfimg(q)=cfnum(q)./cfden(q);
    cfimg=reshape(cfimg,nz,nx);
end
end

function m=metrics(e,cm,bm)
a=e(cm); b=e(bm);
m.CNR=abs(mean(b)-mean(a))/sqrt(var(b,1)+var(a,1));
m.CR=20*log10(max(mean(a),realmin)/max(mean(b),realmin));
edges=linspace(min([a;b]),max([a;b]),101);
pa=histcounts(a,edges,'Normalization','probability');
pb=histcounts(b,edges,'Normalization','probability');
m.gCNR=1-sum(min(pa,pb));
end

function [clean_out,noisy_out,pn_band]=band_wiener(clean,noisy,noise_unit,fs,edges,noise_power,win,min_gain)
% Zero-phase FIR band isolation followed by local variance-subtraction Wiener gain.
% Filter the independent unit-noise realization to estimate its band variance.
order=64;
b=fir1(order,edges/(fs/2),'bandpass',hann(order+1));
clean_band=filtfilt(b,1,double(clean));
noise_band=sqrt(noise_power)*filtfilt(b,1,double(noise_unit));
pn_band=mean(noise_band(:).^2);
noisy_band=filtfilt(b,1,double(noisy));
local_power=movmean(noisy_band.^2,win,1);
gain=max(local_power-pn_band,0)./(local_power+eps);
gain=max(gain,min_gain);
gain=movmean(gain,5,1);
% Apply the estimated gain to both noisy and clean reference RF.
clean_out=clean_band.*gain;
noisy_out=noisy_band.*gain;
end

function value=region_snr(clean_img,noisy_img,mask)
% Complex-image signal-to-error ratio against the method-matched clean result.
signal=clean_img(mask); error_image=noisy_img(mask)-signal;
value=10*log10(sum(abs(signal).^2)/(sum(abs(error_image).^2)+eps));
end

function [img,img_clean]=adaptive_reconstruct(noisy,clean,noise_unit,G,fs,edges,noise_power,win, ...
    min_gain,gamma,cffloor,adaptive_floor)
nBands=size(edges,1); band_imgs=cell(nBands,1); clean_imgs=cell(nBands,1);
band_cfs=cell(nBands,1); clean_cfs=cell(nBands,1);
for bi=1:nBands
    [clean_band,noisy_band]=band_wiener(clean,noisy,noise_unit,fs,edges(bi,:), ...
        noise_power,win,min_gain);
    [clean_imgs{bi},clean_cfs{bi}]=beamform(clean_band,G,true,gamma,cffloor);
    [band_imgs{bi},band_cfs{bi}]=beamform(noisy_band,G,true,gamma,cffloor);
end
stack=cat(3,band_imgs{:}); cf=cat(3,band_cfs{:});
power=abs(stack).^2; total=sum(power,3);
weights=power./(total+eps).*max(cf,adaptive_floor).^2;
weights=weights./(sum(weights,3)+eps);
positive=total(isfinite(total)&total>0);
if isempty(positive), floor_img=0; else, floor_img=0.05*prctile(positive,20); end
gate=total./(total+floor_img+eps); gate(~isfinite(gate))=1;
gate=min(max(gate.*(0.5+0.5*sqrt(max(mean(cf,3),0))),0.35),1);
img=sum(weights.*stack,3).*gate;
clean_stack=cat(3,clean_imgs{:}); clean_cf=cat(3,clean_cfs{:});
clean_power=abs(clean_stack).^2; clean_total=sum(clean_power,3);
clean_weights=clean_power./(clean_total+eps).*max(clean_cf,adaptive_floor).^2;
clean_weights=clean_weights./(sum(clean_weights,3)+eps);
clean_positive=clean_total(isfinite(clean_total)&clean_total>0);
if isempty(clean_positive), clean_floor=0; else, clean_floor=0.05*prctile(clean_positive,20); end
clean_gate=clean_total./(clean_total+clean_floor+eps); clean_gate(~isfinite(clean_gate))=1;
clean_gate=min(max(clean_gate.*(0.5+0.5*sqrt(max(mean(clean_cf,3),0))),0.35),1);
img_clean=sum(clean_weights.*clean_stack,3).*clean_gate;
end
