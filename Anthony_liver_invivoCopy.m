%Pontificia universidad catolica del Peru%
%Anthony Carrera%

clc;close all; clear all;

addpath('./functions_att');
addpath('./bm3d');
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultTextFontSize',12);
font = 12;
% Data block size
db_size= [15];

blocksize= db_size;

%windows filter
winsize=0.5;
window_type = 5;
saran_layer=0;
type_inclusion = 1;
rf_class=102;
% mu_values for the RSLD estimation
mu = 10.^(1:5);
mu_red = 10.^(3:9);
% groundtruth of the reference
att_ref_dB = [0 , 0.54 , 0];

%scale of color in the results
bottom = 0;
top = 1.2;

x_inf=0;
x_sup=6;
z_inf=0;
z_sup=8.8;

%SELECT REGION ANALISIS
%rectangle 1
c1x = 3;
c1y =7.0;
%rectangle 2
c2x = 3;
c2y=7.2;
% size of the rectangles
width=4.5;
height=2.1;
width2=4.5;
height2=1;

%range of the B-mode initial image
dyn_range =65;
overall_limits=[-inf inf 0 0.12];
%freq of analysis
freq_L = 4;
freq_H = 8;
ground_truth = [0.49 0.49];
overlap_pc =0.8;

%% load the sample and the reference

out = load("Q:\smerino\REDjournalResults\rf\invivoLiver.mat");
x = out.xRf; z = out.zRf;
Im_db = out.bMode;
Im_db_ori = Im_db;
x_ori = out.xBm; z_ori = out.zBm;
fs = out.fs;
sam1 = out.rf;
ref1 = out.ref(:,:,1);
ref2 = out.ref(:,:,1);
ref3 = out.ref(:,:,1);
ref4 = out.ref(:,:,1);
clear out

dx=(x_ori(end)-x_ori(1))/(length(x_ori)-1);
dz=(z_ori(end)-z_ori(1))/(length(z_ori)-1);
c0 = 1540;


figure(1); set(1,'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
set(1,'Color',[1 1 1]);
set(gca,'FontSize',font);
set(1,'Color',[1 1 1]);
imagesc(x,z,Im_db); axis image; colormap gray; caxis([0-dyn_range 0]);
hb1 = colorbar;
ylabel(hb1,'dB','FontSize', font)
%title('B-mode image');
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
set(gca,'FontSize',font);

%% Region of interest selection
ind_x = x_inf <= x_ori*100 & x_ori*100 <= x_sup;
ind_z = z_inf <= z_ori*100 & z_ori*100 <= z_sup;

freq_L = freq_L*1e6;   % (Hz)
freq_H = freq_H*1e6;   % (Hz)

wl = c0/mean([freq_L freq_H]);   % Wavelength (m)
disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);

% Number of repetitions of datablocks within a single
% datablock
rpt = 1/(1-overlap_pc);   % r = rep = 1/(1-overlap)
rpt = round(rpt);
overlap = 1 - (1/rpt);

% Part without lateral overlap
wx  = round((blocksize*wl)/(dx*rpt));
% Number of lines laterally = repetitions * block without repetition
nx  = rpt*wx;
new_blocksize = round(nx*dx/(wl));

% RF data columns
L2   = size(sam1,2);
% Blocksize colums
ncol = floor((L2-(rpt-1)*wx)/wx);
sam1 = sam1(:,1:wx*(ncol+rpt-1));
% Actual rf data columns
L2 = size(sam1,2);
x  = x(1:L2);

xi = 1;
xf = L2;
x0 = (xi:wx:xf+1-nx);
x_ACS = x(1,x0+round(nx/2));
nx_ACS = nx*dx;

n  = length(x0);

wz = floor(nx*dx/(dz*rpt));
nz = rpt*wz;

% winsize: Percentage of window (max 0.5)
% nw: Samples of each window axially
nw = 2*floor(winsize*nx*dx/(2*dz)) - 1 ;
L = (nz - nw)*dz*100;   % (cm)

NFFT = 2^(nextpow2(nw)+2);
band = fs*linspace(0,1,NFFT)';   % [Hz] Band of frequencies
rang = (floor(freq_L/fs*NFFT)+1:round(freq_H/fs*NFFT));   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
L3 = length(rang);
p = L3;

L1   = size(sam1,1);   % RF data: rows
nrow = floor((L1-(rpt-1)*wz)/wz);        % Blocksize rows
sam1 = sam1(1:wz*(nrow+rpt-1),:);
L1   = size(sam1,1);   % RF data: rows
z    = z(1:L1);

zi = 1;
zf = L1;
z0 = (zi:wz:zf+1-nz);
m  = length(z0);

z_ACS = z(z0+round(nz/2));
nz_ACS = nz*dz;

z0p = z0 + (nw-1)/2;
z0d = z0 + (nz-1) - (nw-1)/2;

%%
ref1 = ref1(1:L1,1:L2);
ref2 = ref2(1:L1,1:L2);
ref3 = ref3(1:L1,1:L2);
ref4 = ref4(1:L1,1:L2);

disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);


%window_type = input('Windowing type (Tukey-0.25 (5), Hamming (6), Boxcar (7)): ');
switch window_type
    case 5
        windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25
    case 6
        windowing = hamming(nw);   % Hamming
    case 7
        windowing = rectwin(nw);   % Boxcar

end

% Windowing neccesary before Fourier transform
windowing = windowing*ones(1,nx);

figure(3); set(3,'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
figure(6); set(6,'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
for jj=1:n
    for ii=1:m

        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = sam1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        blockP(ii,jj) = effective_lines(sam1(zp-(nw-1)/2:2:zp+(nw-1)/2,xw:xw+nx-1) );
        blockD(ii,jj) = effective_lines( sam1(zd-(nw-1)/2:2:zd+(nw-1)/2,xw:xw+nx-1) );
        blockT(ii,jj) = effective_lines( sam1(zp-(nw-1)/2:2:zd+(nw-1)/2,xw:xw+nx-1) );

        [tempSp,temp_psnr_Sp] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd,temp_psnr_Sd] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);

        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));


        if jj==floor(3*n/6)
            if ii==ceil(m/6)
                figure(3)
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k');
                title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);

                figure(6)
                spmax = max(tempSp((1:NFFT/2+1)));
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'k');
                %title('Spectrum');
                xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);

                %                                 figure(7); set(7,'units','normalized','outerposition',[0 0 1 1]);
                %                                 s1 = sub_block_p(:,round(end/2));
                %                                 smax = max(abs(s1));
                %                                 set(gca,'FontSize',font);
                %                                 plot(s1/smax,'k');
                %                                 %title('Gated ultrasound signal');
                %                                 xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
                %                                 %axis([0 25 -70 0]);
                %                                 set(gca,'FontSize',font);
                %                                 axis tight;
                %                                 ylim([-1 1])
            elseif ii==round(m/2)
                figure(3)
                hold on;
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'r');
                title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);

                figure(6)
                hold on;
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'r');
                %title('Spectrum');
                xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);

                %                                 figure(8); set(8,'units','normalized','outerposition',[0 0 1 1]);
                %                                 s2 = sub_block_p(:,round(end/2));
                %                                 %hold on;
                %                                 set(gca,'FontSize',font);
                %                                 plot(s2/smax,'r');
                %                                 %title('Gated ultrasound signal');
                %                                 xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
                %                                 %axis([0 25 -70 0]);
                %                                 set(gca,'FontSize',font);
                %                                 axis tight;
                %                                 ylim([-1 1])

            elseif ii==floor(5*m/6)
                figure(3)
                hold on;
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'b');   %
                title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);
                %pause
                legend('Top','Half','Bottom');

                figure(6)
                hold on;
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'b');   %
                %title('Spectrum');
                xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity normalized (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);
                %pause
                legend('1 cm above focal dist.','At focal distance','1 cm below focal dist.');

                %                                 figure(9); set(9,'units','normalized','outerposition',[0 0 1 1]);
                %                                 %hold on;
                %                                 s3 = sub_block_p(:,round(end/2));
                %                                 set(gca,'FontSize',font);
                %                                 plot(s3/smax,'b');
                %                                 %title('Gated ultrasound signal');
                %                                 xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
                %                                 %axis([0 25 -70 0]);
                %                                 set(gca,'FontSize',font);
                %                                 %legend('Top','Half','Bottom');
                %                                 axis tight;
                %                                 ylim([-1 1])

            end

        end

    end

end
figure(3)
hold on;
plot(band((1:NFFT/2+1))*1e-6,-20*ones(NFFT/2+1,1),'k--');
plot(band((1:NFFT/2+1))*1e-6,-15*ones(NFFT/2+1,1),'k--');
plot(band((1:NFFT/2+1))*1e-6,-10*ones(NFFT/2+1,1),'k--');

%% Correlation in reference phantom
blockP_avg = mean2(blockP);
blockD_avg = mean2(blockD);

rxx = 5;
index_block_row = (1:rxx:size(blockT,1));
index_block_col = (1:rxx:size(blockT,2));

switch rf_class
    case {3,6,14,102,106}
        att_ref = (att_ref_dB(1)*(f.^2) + att_ref_dB(2)*f + att_ref_dB(3) )/8.686;

    case {2,5,7,4,103}
        if phantom_ref == 1
            att_ref = (0.0076*(f).^2 +0.1189*(f) -0.0319)/8.686;   % [Np/cm]
        elseif phantom_ref == 2
            att_ref = (0.0057*(f).^2 +0.4432*(f) -0.1000)/8.686;   % [Np/cm]
        end
end


%% Diffraction compensation

% Attenuation reference
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Four 'samples' of the reference phantom

for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);
        % Reference 1
        sub_block_p = ref1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp1,temp_psnr_Sp1] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd1,temp_psnr_Sd1] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 2
        sub_block_p = ref2(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref2(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp2,temp_psnr_Sp2] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd2,temp_psnr_Sd2] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 3
        sub_block_p = ref3(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref3(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp3,temp_psnr_Sp3] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd3,temp_psnr_Sd3] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 4
        sub_block_p = ref4(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref4(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp4,temp_psnr_Sp4] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd4,temp_psnr_Sd4] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);

        tempSp = 1/4*(tempSp1 + tempSp2 + tempSp3 + tempSp4);
        tempSd = 1/4*(tempSd1 + tempSd2 + tempSd3 + tempSd4);
        temp_psnr_Sp = 1/4*(temp_psnr_Sp1 + temp_psnr_Sp2 + temp_psnr_Sp3 + temp_psnr_Sp4);
        temp_psnr_Sd = 1/4*(temp_psnr_Sd1 + temp_psnr_Sd2 + temp_psnr_Sd3 + temp_psnr_Sd4);

        Sp_ref(ii,jj,:) = (tempSp(rang));
        Sd_ref(ii,jj,:) = (tempSd(rang));

        if jj==floor(2.5*n/5)
            if ii==ceil(m/10)
                figure(3)
                %subplot(311);
                %set(gca,'FontSize',sl);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k--');   %
                title('Spectrum REFERENCE'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);
            elseif ii==round(m/2)
                figure(3);
                hold on;
                %subplot(312);
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'r--');   %
                title('Spectrum REFERENCE'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);
            elseif ii==floor(9*m/10)
                figure(3);
                hold on;
                %subplot(313);
                set(gca,'FontSize',font);
                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'b--');   %
                title('Spectrum REFERENCE'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                axis([0 25 -70 0]);
                set(gca,'FontSize',font);
                %pause
                legend(' REF Top','REF Half','REF Bottom');
                hold off
            end

        end


    end
end

diffraction_compensation = ( log(Sp_ref) - log(Sd_ref) ) - 4*L*att_ref_map;

figure(3); hold off;

%% Au = b
b = (log(Sp) - log(Sd)) - (diffraction_compensation);
% watch curves
temp_ratio = 0; count_ratio = 0;
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

%% Regularization using RED  with median
disp(['Start RED median']);
for q = 1:length(mu_red)
    disp(['RED median:  ', num2str(q)]);
    mu_aux = mu_red(q) ;
    [err_fp2,nfp2 ,u2]  =  admm_red_median(A'*A,A'*b(:),mu_aux,0.001,size(A'*b(:),1),1500,4,1,5,m,n,mu_aux/1.0);
    BS2 = u2(1:end/2); %CS = u(end/2+1:end);
    BS2 = 8.686*BS2;   % [dB.cm^{-1}.MHz^{-1}]
    BRED(:,:,q) = reshape(BS2,m,n);
    BSRED_interp(:,:,q) = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BRED(:,:,q), x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 200, 1, [bottom top], '(a)', c1x ,c1y , width , height , c2x , c2y , width2 , height2 );

    pause(0.1)
end

%% Measurement ROI
for ii = 1:size( BSRED_interp,3)
    R1_REDSLD(:,:,ii) = rect_center( BSRED_interp(:,:,ii),x,z,c1x,c1y,width,height);
    R2_REDSLD(:,:,ii) = rect_center(BSRED_interp(:,:,ii),x,z,c2x,c2y,width2,height2);
end

for ii = 1:size(BSRED_interp,3)
    REDSLD.result(ii,:,1) = [mu_red(ii) nanmean(reshape(R1_REDSLD(:,:,ii),[],1)) nanstd(reshape(R1_REDSLD(:,:,ii),[],1))];
    REDSLD.result(ii,:,2) = [mu_red(ii) nanmean(reshape(R2_REDSLD(:,:,ii),[],1)) nanstd(reshape(R2_REDSLD(:,:,ii),[],1))];

    REDSLD.error{ii,1} = (reshape(R1_REDSLD(:,:,ii),[],1) - ground_truth(1))/ground_truth(1);
    mae = (abs(reshape(R1_REDSLD(:,:,ii),[],1) - ground_truth(1)))/ground_truth(1);

    disp(['Regularized SLD RED. median ',num2str( REDSLD.result(ii,1,1),'%8.3e'),'. '...
        'R1: ',num2str( REDSLD.result(ii,2,1),'%3.2f'),' +/- ',num2str( REDSLD.result(ii,3,1),'%3.2f'),' dB.cm^{-1}.MHz^{-1}. '...
        'R2: ',num2str( REDSLD.result(ii,2,2),'%3.2f'),' +/- ',num2str( REDSLD.result(ii,3,2),'%3.2f'),' dB.cm^{-1}.MHz^{-1}.'...
        'CV:' , num2str(REDSLD.result(ii,3,1)/REDSLD.result(ii,2,1),'%3.2f') ...
        ' MPE:',num2str(  abs(nanmean( REDSLD.error{ii,1})),'%3.2f' ) ...
        ' MAE:',num2str( abs(nanmean( mae)) ,'%3.2f' )

        ]);


    % inclusion mean and sd
    REDSLD.error{ii,1} = (abs(reshape(R1_REDSLD(:,:,ii),[],1) - ground_truth(1)))/ground_truth(1);
    REDSLD.error{ii,2} = [abs(nanmean( REDSLD.error{ii,1})) nanstd( REDSLD.error{ii,1}) 100*abs(nanmean(REDSLD.error{ii,1})) 100*nanstd( REDSLD.error{ii,1}) ];
    REDSLD.vrr(ii,1) = 100*( 1 - nanvar(reshape(R1_REDSLD(:,:,ii),[],1))/ nanvar(reshape(R1_REDSLD(:,:,1),[],1)) );
    % background mean and sd
    REDSLD.error{ii,3} = (reshape(R2_REDSLD(:,:,ii),[],1) - ground_truth(2))/ground_truth(2);
    REDSLD.error{ii,4} = [abs(nanmean( REDSLD.error{ii,3})) nanstd( REDSLD.error{ii,3}) 100*abs(nanmean(REDSLD.error{ii,3})) 100*nanstd( REDSLD.error{ii,3}) ];
    REDSLD.vrr(ii,2) = 100*( 1 - nanvar(reshape(R2_REDSLD(:,:,ii),[],1))/ nanvar(reshape(R2_REDSLD(:,:,1),[],1)) );

    REDSLD.bias(ii,1) =  REDSLD.error{ii,2}(3);
    REDSLD.bias(ii,2) =  REDSLD.error{ii,4}(3);

    REDSLD.cv(ii,1) =100*REDSLD.result(ii,3,1)/REDSLD.result(ii,2,1);
    REDSLD.cv(ii,2) = 100*REDSLD.result(ii,3,2)/REDSLD.result(ii,2,2);

end


%% Figures
figure(1012); set(1012,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(121)
set(gca,'FontSize',font,'xscale','log');
%set(gca,'FontSize',font,'xscale','log');
errorbar(mu_red, REDSLD.result(:,2,2), REDSLD.result(:,3,2),'r-.d','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);

%xlim([10^-5 10^11]);

hold on;
errorbar(mu_red, REDSLD.result(:,2,1), REDSLD.result(:,3,1),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);

%plot(mu1, ground_truth(1)*ones(size(mu1)),'r--s','LineWidth',3);
plot(mu_red, ground_truth(1)*ones(size(mu_red)),'b--','LineWidth',3);
plot(mu_red, ground_truth(2)*ones(size(mu_red)),'r--','LineWidth',3);
%%plot(ones(1,17)*10^1.25, linspace(bottom,top,17), 'k--','LineWidth',3);
%%plot(ones(1,17)*10^3.75, linspace(bottom,top,17), 'k--','LineWidth',3);

%xlim([results_BR(1,1,1) results_BR(end,1,1)]);
axis tight
xlim([REDSLD.result(1,1,1) REDSLD.result(end,1,1)]);
ylim([bottom top]);
%title('(a)');
xlabel('\bfRegularization parameter \mu'); ylabel('\bfMean and SD (dB.cm^{-1}.MHz^{-1})');
set(gca,'FontSize',font);
%legend('Regularized SLD', 'Through-transmission')
set(gca,'XScale','log');
%errorbarlogx();
h = legend('Background','Inclusion','Location','NorthEast');
set(h,'fontsize',font);



figure(1014); set(1014,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;
plot(mu_red, REDSLD.bias(:,2),'r-d','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);
plot(mu_red, REDSLD.bias(:,1),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);
%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfMPE (%)');
set(gca,'FontSize',font);
h = legend('Background','Inclusion','Location','NorthWest');
set(h,'fontsize',font);
plot(mu_red, 20*ones(size(mu_red)),'k--','LineWidth',3);
plot(mu_red, 10*ones(size(mu_red)),'k--','LineWidth',3);

figure(1015); set(1015,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;
plot(mu_red, REDSLD.cv(:,2),'r-d','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);
plot(mu_red, REDSLD.cv(:,1),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);
%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfCV (%)');
set(gca,'FontSize',font);
h = legend('Background','Inclusion','Location','NorthWest');
set(h,'fontsize',font);
plot(mu_red, 20*ones(size(mu_red)),'k--','LineWidth',3);
plot(mu_red, 10*ones(size(mu_red)),'k--','LineWidth',3);

%% Regularization: Au = b
tol = 1e-3;
%mu = logspace(2,3,3)';
%mu = logspace(1,2,3)';
%mu = [0 1];

clear mask
mask1 = ones(m,n);
mask4 = speye(m*n,m*n);
mask_diff = speye(2*m*n,2*m*n);

for kk=1:p
    mask(:,:,kk)=mask1(:,:);
end

for mm = 1:length(mu)

    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    BR_interp(:,:,mm) = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BR(:,:,mm), x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 100, 1, [bottom top], ' (b) ',c1x , c1y ,width , height , c2x , c2y ,width2, height2);
    % R1_RSLDg = square_centerg(BR_interp(:,:,mm),x,z,c1x,c1y,square_L);

    %figure(300+mm); imagesc(R1_RSLDg); colormap('jet') ;colorbar;
    %BR_interp(:,:,mm) = attROT(BR_interp(:,:,mm), x, z, angle, rot_lat, rot_axi, bottom, top);

end



%% Selection of regions with differeent attenuation: layers or inclusions
% Type 1: One circular inclusion
% R1: Region 1: Inclusion (Standards SLD, Regularized SLD)
% R2: Region 2: Background (Standards SLD, Regularized SLD)
% Type_inclusion = input('Type of regions (No Inc. (1), Line(2), Inc_Back(3), Inc_Back_V2 (4)): ');

for ii = 1:size(BR_interp,3)

    R1_RSLD(:,:,ii) = rect_center(BR_interp(:,:,ii),x,z,c1x,c1y,width,height);

    R2_RSLD(:,:,ii) = rect_center(BR_interp(:,:,ii),x,z,c2x,c2y,width2,height2);

end

for ii = 1:size(BR_interp,3)
    RSLD.result(ii,:,1) = [mu(ii) nanmean(reshape(R1_RSLD(:,:,ii),[],1)) nanstd(reshape(R1_RSLD(:,:,ii),[],1))];
    RSLD.result(ii,:,2) = [mu(ii) nanmean(reshape(R2_RSLD(:,:,ii),[],1)) nanstd(reshape(R2_RSLD(:,:,ii),[],1))];

    RSLD.error{ii,1} = (reshape(R1_RSLD(:,:,ii),[],1) - ground_truth(1))/ground_truth(1);
    mae = (abs(reshape(R1_RSLD(:,:,ii),[],1) - ground_truth(1)))/ground_truth(1);

    disp(['Regularized SLD. ',num2str(RSLD.result(ii,1,1),'%8.3e'),'. '...
        'R1: ',num2str(RSLD.result(ii,2,1),'%3.2f'),' +/- ',num2str(RSLD.result(ii,3,1),'%3.2f'),' dB.cm^{-1}.MHz^{-1}. '...
        'R2: ',num2str(RSLD.result(ii,2,2),'%3.2f'),' +/- ',num2str(RSLD.result(ii,3,2),'%3.2f'),' dB.cm^{-1}.MHz^{-1}.'...
        ' CV:',num2str(  RSLD.result(ii,3,1) /RSLD.result(ii,2,1),'%3.2f' )...
        ' MPE:',num2str(  abs(nanmean( RSLD.error{ii,1})),'%3.2f' ) ...
        ' MAE:',num2str( abs(nanmean( mae)) ,'%3.2f' )
        ]);
    % inclusion
    RSLD.error{ii,1} = (abs(reshape(R1_RSLD(:,:,ii),[],1) - ground_truth(1)))/ground_truth(1);
    RSLD.error{ii,2} = [abs(nanmean(RSLD.error{ii,1})) nanstd(RSLD.error{ii,1}) 100*abs(nanmean(RSLD.error{ii,1})) 100*nanstd(RSLD.error{ii,1}) ];
    RSLD.vrr(ii,1) = 100*( 1 - nanvar(reshape(R1_RSLD(:,:,ii),[],1))/ nanvar(reshape(R1_SLD(:,:,1),[],1)) );
    % background
    RSLD.error{ii,3} = (reshape(R2_RSLD(:,:,ii),[],1) - ground_truth(2))/ground_truth(2);
    RSLD.error{ii,4} = [abs(nanmean(RSLD.error{ii,3})) nanstd(RSLD.error{ii,3}) 100*abs(nanmean(RSLD.error{ii,3})) 100*nanstd(RSLD.error{ii,3}) ];
    RSLD.vrr(ii,2) = 100*( 1 - nanvar(reshape(R2_RSLD(:,:,ii),[],1))/ nanvar(reshape(R2_SLD(:,:,1),[],1)) );

    RSLD.bias(ii,1) = RSLD.error{ii,2}(3);
    RSLD.bias(ii,2) = RSLD.error{ii,4}(3);

    RSLD.cv(ii,1) = 100*RSLD.result(ii,3,1)/RSLD.result(ii,2,1);
    RSLD.cv(ii,2) = 100*RSLD.result(ii,3,2)/RSLD.result(ii,2,2);

end

% table only with mu = 10^2.5 results of MPE, SDPE and CNR
RSLD.table.inc(1,:)  = [SLD.result(:,:,1) , SLD.error{2}(1:2) , RSLD.result(2,2:3,1), RSLD.error{2,2}(1:2) ];
RSLD.table.back(1,:) = [SLD.result(:,:,2) , SLD.error{4}(1:2) , RSLD.result(2,2:3,2), RSLD.error{2,4}(1:2) ];
RSLD.table.vrr(1,1)  = RSLD.vrr(2,1);

%%
figure(12); set(12,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(121)
set(gca,'FontSize',font,'xscale','log');
%set(gca,'FontSize',font,'xscale','log');
errorbar(mu, RSLD.result(:,2,2), RSLD.result(:,3,2),'r-.d','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);

%xlim([10^-5 10^11]);

hold on;
errorbar(mu, RSLD.result(:,2,1), RSLD.result(:,3,1),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);

%plot(mu1, ground_truth(1)*ones(size(mu1)),'r--s','LineWidth',3);
plot(mu, ground_truth(1)*ones(size(mu)),'b--','LineWidth',3);
plot(mu, ground_truth(2)*ones(size(mu)),'r--','LineWidth',3);
%plot(ones(1,17)*10^1.25, linspace(bottom,top,17), 'k--','LineWidth',3);
% plot(ones(1,17)*10^3.75, linspace(bottom,top,17), 'k--','LineWidth',3);

%xlim([results_BR(1,1,1) results_BR(end,1,1)]);
axis tight
xlim([RSLD.result(1,1,1) RSLD.result(end,1,1)]);
ylim([bottom top]);
%title('(a)');
xlabel('\bfRegularization parameter \mu'); ylabel('\bfMean and SD (dB.cm^{-1}.MHz^{-1})');
set(gca,'FontSize',font);
%legend('Regularized SLD', 'Through-transmission')
set(gca,'XScale','log');
%errorbarlogx();
h = legend('Background','Inclusion','Location','NorthEast');
set(h,'fontsize',font);



figure(14); set(14,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;
plot(mu, RSLD.bias(:,2),'r-d','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);
plot(mu, RSLD.bias(:,1),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);
%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfMPE (%)');
set(gca,'FontSize',font);
h = legend('Background','Inclusion','Location','NorthWest');
set(h,'fontsize',font);
plot(mu, 20*ones(size(mu)),'k--','LineWidth',3);
plot(mu, 10*ones(size(mu)),'k--','LineWidth',3);



figure(15); set(15,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;
plot(mu, RSLD.cv(:,2),'r-d','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);
plot(mu, RSLD.cv(:,1),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);
%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfCV (%)');
set(gca,'FontSize',font);
h = legend('Background','Inclusion','Location','NorthWest');
set(h,'fontsize',font);
plot(mu_red, 20*ones(size(mu_red)),'k--','LineWidth',3);
plot(mu_red, 10*ones(size(mu_red)),'k--','LineWidth',3);


figure(70); set(70,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;
plot(mu, RSLD.cv(:,1),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);
plot(mu_red, REDSLD.cv(:,1),'r-s','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);


%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfCV (%) inclusion');
set(gca,'FontSize',font);
h = legend('TV','RED-median','RED-NLM','RED-BM3D','Location','NorthWest');
set(h,'fontsize',font);
plot(mu_red, 20*ones(size(mu_red)),'k--','LineWidth',3);
plot(mu_red, 10*ones(size(mu_red)),'k--','LineWidth',3);


figure(71); set(71,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;
plot(mu, RSLD.cv(:,2),'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);
plot(mu_red, REDSLD.cv(:,2),'r-s','LineWidth',3,'MarkerFaceColor','r','MarkerSize',20);


%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfCV (%)');
set(gca,'FontSize',font);
h = legend('TV','RED-median','RED-NLM','RED-BM3D','Location','NorthWest');
set(h,'fontsize',font);
plot(mu_red, 20*ones(size(mu_red)),'k--','LineWidth',3);
plot(mu_red, 10*ones(size(mu_red)),'k--','LineWidth',3);

color1 = '#1f77b4'; % Azul
color2 = '#ff7f0e'; % Naranja
color3 = '#2ca02c'; % Verde

rgb1 = sscanf(color1(2:end), '%2x%2x%2x', [1 3]) / 255;
rgb2 = sscanf(color2(2:end), '%2x%2x%2x', [1 3]) / 255;
rgb3 = sscanf(color3(2:end), '%2x%2x%2x', [1 3]) / 255;

figure(73); set(73,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;

plot(mu, RSLD.bias(:,1),'s-','LineWidth',3,'MarkerFaceColor',rgb1,'MarkerSize',20, 'Color', rgb1);
plot(mu_red, REDSLD.bias(:,1),'s-','LineWidth',3,'MarkerFaceColor',rgb2,'MarkerSize',20, 'Color', rgb2);


%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfMAE (%)');
set(gca,'FontSize',font);
h = legend('TV','RED-median','RED-NLM','RED-BM3D','Location','NorthWest');
set(h,'fontsize',font);
plot(mu_red, 20*ones(size(mu_red)),'k--','LineWidth',3);
plot(mu_red, 10*ones(size(mu_red)),'k--','LineWidth',3);



figure(74); set(74,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;

plot(mu, RSLD.bias(:,2),'s-','LineWidth',3,'MarkerFaceColor',rgb1,'MarkerSize',20, 'Color', rgb1);
plot(mu_red, REDSLD.bias(:,2),'s-','LineWidth',3,'MarkerFaceColor',rgb2,'MarkerSize',20, 'Color', rgb2);


%plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',3);
%plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',3);
%ylim([0.2 1.4]);
%ylim([0 2.0]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfMPE (%)');
set(gca,'FontSize',font);
h = legend('TV','RED-median','RED-NLM','Location','NorthWest');
set(h,'fontsize',font);
plot(mu_red, 20*ones(size(mu_red)),'k--','LineWidth',3);
plot(mu_red, 10*ones(size(mu_red)),'k--','LineWidth',3);


figure(75); set(75,'units','normalized','outerposition',[0 0 1 1]); box on;
%subplot(212);
set(gca,'FontSize',font,'xscale','log');
%xlim([10^0 10^5]); % semilogx;
hold on;

errorbar(mu, RSLD.result(:,2,2), RSLD.result(:,3,2)/2,'b-.d','LineWidth',3,'MarkerFaceColor','b','MarkerSize',20);
errorbar(mu_red, REDSLD.result(:,2,2), REDSLD.result(:,3,2)/2,'g-.d','LineWidth',3,'MarkerFaceColor','g','MarkerSize',20);

ylim([-0.2 1.5]);
xlabel('\bfRegularization parameter \mu');
ylabel('\bfMean and SD (%)');
set(gca,'FontSize',font);
h = legend('TV','RED-median','RED-NLM','Location','NorthWest');
set(h,'fontsize',font);



save_all_figures_to_directory("compliv20")
save(['compliv20','\stats_RSLD.mat'],'REDSLD','RSLD');

close all;

%% Auxiliary functions
function save_all_figures_to_directory(dir_name)

figlist=findobj('type','figure');

for i=1:numel(figlist)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.fig']));
    figure(figlist(i))
    set(gcf,'PaperPositionMode','auto')
    %pause(2)
    saveas(figlist(i),fullfile(dir_name,['figure_' num2str(figlist(i).Number) '.png']));
    %pause(2)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.eps']));

end

end

function Neff = effective_lines(data_block)


[~, N] =size(data_block);
RHO = corrcoef(data_block);
%v = 0:N-1;
%factor =toeplitz([v(1) fliplr(v(2:end))], v);

%val = factor.*(RHO.^2);
%val = sum(val);

rho = diagSumCalc(RHO,1);
rho = rho(1:N-1)./(1:N-1);

val = (rho.^2)*(1:N-1)';
Neff = N./( 1 + 2/N*( val ) );
%[mean(Neff) std(Neff) median(Neff) mode(Neff)]

% for ii = 1:d2-1,
%    y = data_block(:,ii);
%    bb = corrcoef(x,y);
%    rho(ii) = bb(1,2);
%    factor(ii) = ii;
%    %rho(ii) =  sum( xcorr(x,y') )/( norm(x)*norm(y) );
%
% end
%
% N = d2
% Neff = N/( 1 + 2/N*( factor*(rho').^2 ) )

end

function [diagSum] = diagSumCalc(squareMatrix, LLUR0_ULLR1)
%
% Input: squareMatrix: A square matrix.
%        LLUR0_ULLR1:  LowerLeft to UpperRight addition = 0
%                      UpperLeft to LowerRight addition = 1
%
% Output: diagSum: A vector of the sum of the diagnols of the matrix.
%
% Example:
%
% >> squareMatrix = [1 2 3;
%                    4 5 6;
%                    7 8 9];
%
% >> diagSum = diagSumCalc(squareMatrix, 0);
%
% diagSum =
%
%       1 6 15 14 9
%
% >> diagSum = diagSumCalc(squareMatrix, 1);
%
% diagSum =
%
%       7 12 15 8 3
%
% Written by M. Phillips
% Oct. 16th, 2013
% MIT Open Source Copywrite
% Contact mphillips@hmc.edu fmi.
%

if (nargin < 2)
    disp('Error on input. Needs two inputs.');
    return;
end

if (LLUR0_ULLR1 ~= 0 && LLUR0_ULLR1~= 1)
    disp('Error on input. Only accepts 0 or 1 as input for second condition.');
    return;
end

[M, N] = size(squareMatrix);

if (M ~= N)
    disp('Error on input. Only accepts a square matrix as input.');
    return;
end

diagSum = zeros(1, M+N-1);

if LLUR0_ULLR1 == 1
    squareMatrix = rot90(squareMatrix, -1);
end

for i = 1:length(diagSum)
    if i <= M
        countUp = 1;
        countDown = i;
        while countDown ~= 0
            diagSum(i) = squareMatrix(countUp, countDown) + diagSum(i);
            countUp = countUp+1;
            countDown = countDown-1;
        end
    end
    if i > M
        countUp = i-M+1;
        countDown = M;
        while countUp ~= M+1
            diagSum(i) = squareMatrix(countUp, countDown) + diagSum(i);
            countUp = countUp+1;
            countDown = countDown-1;
        end
    end
end

end

function D = diffraction_JR(fn,fl,band,c0,za)

% Calcula el factor de compensacion por difraccion en la profundidad
% de analisis za a la frecuencia fa (X.Chen)
%
% Entradas:     fn: #focal del transductor
%               fl: Longitud focal del tx (m)
%               band: Frecuencia de analisis (Hz)
%               c: Velocidad del sonido en el medio (m/s)
%               za: Profundidad de analisis (m)
%
% Salidas:      D: Factor de compensacion por difraccion

a  = (fl/(fn*2));                        % Radio del transductor.
Gp = ( (pi*band/c0)*(a^2) ) / fl;           % Preasure gain factor.

for i = 1 : length(Gp)

    if 1/(1+pi/Gp(i))<=za/fl && za/fl<= 1/(1-pi/Gp(i))
        D(i) = ((pi*a^2)/(za^2))* 0.46*exp(-(0.46/pi)*Gp(i).^2*((fl/za)-1).^2);
    else
        %D(i) = ((pi*a^2)/(za.^2))* 1.07*(Gp(i)*((fl/za)-1)).^-2;
        D(i) = ((pi*a^2)/(za.^2))* 1.00*(Gp(i)*((fl/za)-1)).^-2;
    end
    % To PLOT
    %D(i)= D(i)/((pi*a^2)/(za^2));

end

end

function [x,ite_cgs] = cgs2(A,b,tol,maxit,varargin)

if length(varargin) == 1
    x = varargin{1};
else
    x= zeros(size(A,2),1);
end
r = A*x - b;
p = -r;
ite_cgs = 0;

while norm(r,2) > tol && ite_cgs < maxit ;

    alpha = (r'*r)/(p'*(A*p));
    x = x + alpha*p;
    rn = r + alpha*(A*p);
    beta = (rn'*rn)/(r'*r);
    p = -rn + beta*p;
    r = rn;
    ite_cgs = ite_cgs + 1;
end

end

function [swe_Vi] = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BS, x_ate, z_ate, nx_ate, nz_ate, m, n, i, k, c, t, c1x , c1y , width , height , c2x , c2y , width2 , height2)

%sl = 60;
sl = 12;

Bmode = Im_db;

SWS_PD = BS;
%Properties.Width_B = x_ori*1e3;
%Properties.Depth_B = z_ori*1e3;
Properties.Width_B = x*10;
Properties.Depth_B = z*10;

Properties.Width_S = x*10;
Properties.Depth_S = z*10;

rd = c;
% BMode with SWS map superposed
% close all
BW=ones(size(Bmode));
[X1,Y1] = meshgrid(1:size(SWS_PD,2),1:size(SWS_PD,1));
[X2,Y2] = meshgrid(linspace(1,size(SWS_PD,2),size(Bmode,2)),linspace(1,size(SWS_PD,1),size(Bmode,1)));

%[XX1,YY1] = meshgrid(x_ate,z_ate);
%[XX2,YY2] = meshgrid(x,z);

swe_Vi = interp2(X1,Y1,SWS_PD,X2,Y2);


%swe_Vi = interp2(XX1,YY1,SWS_PD,XX2,YY2,'linear');
%figure; imagesc(swe_Vi); colorbar;
%figure; imagesc(SWS_PD); colorbar; caxis([0,1])

% a=Properties.Depth_S;
%a=Properties.Depth_S-Properties.Depth_S(1);
a = Properties.Depth_S;
%figure(500)
figure(i);
set(i,'units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
set(i,'color',[1 1 1]);
transparency=0.55;
alpha_mat = (1-transparency)*ones(size(Bmode));

Bmode = Im_db_ori;
Properties.Width_B = x_ori*1e3;
Properties.Depth_B = z_ori*1e3;

caxis_bmode = [-60 0];


%subimage(1e-1*Properties.Width_B, 1e-1*Properties.Depth_B, 64*mat2gray(Bmode,[-60 0]), gray(64));

subimage(1e-1*Properties.Width_B, 1e-1*Properties.Depth_B, 256*mat2gray(Bmode,caxis_bmode), gray);
%subimage(1e-1*Properties.Width_B, 1e-1*Properties.Depth_B, 64*mat2gray(Bmode), gray);

% imagesc(1e3*Properties.Width_B, 1e3*Properties.Depth_B,Bmode)
hold on
%h=subimage(1e-1*Properties.Width_S, 1e-1*a, 64*mat2gray(BW.*swe_Vi, [rd]), jet(64));
h=subimage(1e-1*Properties.Width_S, 1e-1*a, 256*mat2gray(BW.*swe_Vi,[rd]), jet);
set( h, 'AlphaData', 0.5) ; % .5 transparency
% h=subimage(1e3*Properties.Width_S, 1e3*Properties.Depth_S, swe_Vi);
%set(gca,'FontSize',sl);

rectangle('Position',[c1x - width /2 ,c1y - height/2 ,width,height],'LineWidth',3);
hold on;

% rectangle('Position',[c2x - width2 /2 ,c2y - height2/2 ,width2,height2],'EdgeColor','r','LineWidth',3,'LineStyle','--');
% hold on;



%set( h, 'AlphaData', alpha_mat) ; % .5 transparency
set(gca,'FontSize',sl);

%xlabel('Width [mm]','fontsize',16);ylabel('Depth [mm]','fontsize',sl)
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
set(gca,'FontSize',sl);
%title(['Vibration Frequency ' num2str(Properties.VibFreq) ' Hz'],'fontsize',14)

%h2 = colorbar;
%set(get(h2,'YLabel'),'String','SWS [m/s]','FontSize',sl);
h2 = colorbar;
colormap(jet);
ylabel(h2,'dB.cm^{-1}.MHz^{-1}','FontSize', sl);

%n_tick = rd(end);  % number of tick marks desired on color axis
n_tick = 2;
% get the current colorbar axis limits
cbar_lim = get(h2, 'YLim');
% label colorbar correctly
set(h2, 'YTick',cbar_lim(1):range(cbar_lim)/n_tick:cbar_lim(2));
set(h2,'YTickLabel', [rd(1) rd(end)/2 rd(end)],'FontSize',sl);
set(gca,'FontSize',sl);


hold off

end

function img_inside = rect_center(img,x,z,c1x,c1y,W,H)
%L = 0.6;   % [cm]
u = size(img);
mx = c1x;
my = c1y;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
for jj=1:u(2)
    for ii=1:u(1)
        if abs(x(1) + (jj-1)/(u(2)-1)*(x(end)-x(1)) - mx) > W/2 || abs(z(1) + (ii-1)/(u(1)-1)*(z(end)-z(1)) - my) > H/2
            %if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 > radius^2,
            img(ii,jj) = NaN;
        end
    end
end
img_inside = img;
end


%% andres functions
function [B,C] = AlterOpti_ADMM(A1,A2,b,mu1,mu2,m,n,tol,mask)

p = length(mask)/(m*n);
minimask = reshape(mask,[m n p]);
minimask = minimask(:,:,1);
minimask = minimask(:);
b = mask.*b;
b(isnan(b)) = 0;
A = [A1 A2];
[u,~] = cgs2(A'*A,A'*b,1e-6,20);
B = reshape(u(1:end/2),m,n);
C = reshape(u(end/2+1:end),m,n);
%figure(109); imagesc(8.686*reshape(B,m,n)); colormap pink; caxis([0 1.2])
B = B(:);
C = C(:);
D = 0;
v = 0;

F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);

ite  = 0;
error = 1;

while abs(error) > tol && ite < 20;
    ite = ite + 1;

    rho = 1;
    % First part of ADMM algorithm: B
    % B = IRLS_ANIS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);
    B = IRLS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);
    %F(2*ite,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);
    %error = F(2*ite,1) - F(2*ite-1,1);

    % Second part of ADMM algorithm: C
    % C = IRLS_ANIS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask);
    C = IRLS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask);
    %F(2*ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);
    %error = F(2*ite+1,1) - F(2*ite,1);

    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;

    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);

end

end


% TV Andres Leonel Coila
function [TV] = TVcalc(B,M,N,mask)

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = Dh.^2 + Dv.^2;
P = sqrt(P);
TV = norm(P(:).*mask,1);

end

% TV Andres Leonel Coila - ANISOTROPIC
function [TV] = TVcalc2(B,M,N,mask)

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = abs(Dh) + abs(Dv);
%P = sqrt(P);
TV = norm(P(:).*mask,1);

end

% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,mask,minimask)

[u,~] = cgs2(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol

    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];

    %Dx*X(:) - Dh(:);
    %Dy*X(:) - Dv(:);

    P = Dh.^2 + Dv.^2;
    eps = 0.01;
    P = 2*sqrt(P.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    W = kron(speye(2),omega);

    AtA = A'*A;
    [u] = cgs2( AtA + mu*D'*W*D, A'*b, 1e-6 , 20, u );

    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);
    error = abs(G(ite_irls+1) - G(ite_irls));

end

%figure(909); plot(1:length(G),G);

end

% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_ANIS_TV(b,A,mu,M,N,tol,mask,minimask)

[u,~] = cgs2(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc2(u,M,N,minimask);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol && ite_irls < 20

    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];

    %Dx*X(:) - Dh(:);
    %Dy*X(:) - Dv(:);

    P = Dh.^2 + Dv.^2;
    eps = 0.1;
    P = 2*sqrt(P.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    %W = kron(speye(2),omega);

    Px = abs(Dh + eps);
    Px = 1./Px;
    Px = Px(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Px,0,omega);
    Wx = kron(speye(1),omega);

    Py = abs(Dv + eps);
    Py = 1./Py;
    Py = Py(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Py,0,omega);
    Wy = kron(speye(1),omega);

    AtA = A'*A;
    %mu=5000;
    %[u] = cgs(AtA + mu*D'*W*D, A'*b,1e-6,200);
    [u] = cgs2( AtA + mu*Dx'*Wx*Dx + mu*Dy'*Wy*Dy , A'*b, 1e-6 , 20, u );

    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc2(u,M,N,minimask);
    error = abs(G(ite_irls+1) - G(ite_irls));

end

%figure(909); plot(1:length(G),G);

end

%%

function [res_admm,n_out , out ] = admm_red_median(A,y,lambda,error,N,max_iter,inner_iters,inner_iters2,median,m,ni,beta)

x_est= ones([N 1]);
v_est=ones([N 1]);
u_est=ones([N 1]);

%     x_est= 0.5*abs(x_est)+0.75;
%     v_est= 0.5*abs(v_est)+0.75;
%     u_est= 0.5*abs(u_est)+0.75;

n_out = [x_est v_est u_est];
alpha=1.5;


for k = 1:1:max_iter

    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2

    for j = 1:1:inner_iters
        %                e = A'*(A*x_est-y) + beta*(x_est-(v_est-u_est));
        %                 r = A'*A*e + beta*e;
        %                 mu_opt = (e'*e)/(e'*r);
        %                 x_est = x_est + mu_opt*e;
        %                 %x_est = max( min(x_est, 10), 0);


        b = (A*y) + beta*(v_est - u_est);
        A_x_est = (A'*A*x_est) + beta*x_est;
        res = b - A_x_est;
        %res = (1 / (input_sigma^2))*(A'*(A*x_est-y))+beta*(x_est -(v_est-u_est));
        a_res = (A'*A*res) + beta*res;
        mu_opt = mean(res(:).*res(:))/mean(res(:).*a_res(:));
        x_est = x_est + mu_opt*res;
        %x_est = max( min(x_est, 5), 0);

    end


    % relaxation
    x_hat = alpha*x_est + (1-alpha)*v_est;

    %x_hat = x_est;
    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    % using gradient descent
    for j = 1:1:inner_iters2

        %              v_est = reshape(v_est,30,100);
        %           f_v_est = imnlmfilt(v_est,'DegreeOfSmoothing',10,'SearchWindowSize',15);
        %           f_v_est = f_v_est(:);

        v_est = reshape(x_est,m,2*ni);
        v_est1=v_est((1:m),(1:ni));
        v_est2=v_est((1:m),(ni+1:2*ni));
        f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
        f_v_est2 = medfilt2(v_est2, [median median],'symmetric');
        %
        %             f_v_est1 = BM3D(v_est1, median);
        %             f_v_est2 = BM3D(v_est2, median);


        %             f_v_est1 = imnlmfilt(v_est1,'DegreeOfSmoothing',median);
        %             f_v_est2 = imnlmfilt(v_est2,'DegreeOfSmoothing',median);


        f_v_est = [f_v_est1 f_v_est2];
        f_v_est = f_v_est(:);



        v_est = (beta*(x_hat + u_est) + lambda*f_v_est)/(lambda + beta);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_hat - v_est;

    r(k) = norm(A*x_hat-y)^2/norm(y)^2;

    if r(k) < error
        out = x_hat;
        %n_out=k;
        res_admm = r;
        break
    end


end
out = x_hat;
% n_out=k;
res_admm = r;
end

function [res_admm,n_out , out ] = admm_red_nlm(A,y,lambda,error,N,max_iter,inner_iters,inner_iters2,median,size,comp,m,ni,beta)

x_est= ones([N 1]);
v_est=ones([N 1]);
u_est=ones([N 1]);

%     x_est= 0.5*abs(x_est)+0.75;
%     v_est= 0.5*abs(v_est)+0.75;
%     u_est= 0.5*abs(u_est)+0.75;

n_out = [x_est v_est u_est];
alpha=1.5;


for k = 1:1:max_iter

    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2

    for j = 1:1:inner_iters
        %                e = A'*(A*x_est-y) + beta*(x_est-(v_est-u_est));
        %                 r = A'*A*e + beta*e;
        %                 mu_opt = (e'*e)/(e'*r);
        %                 x_est = x_est + mu_opt*e;
        %                 %x_est = max( min(x_est, 10), 0);


        b = (A*y) + beta*(v_est - u_est);
        A_x_est = (A'*A*x_est) + beta*x_est;
        res = b - A_x_est;
        %res = (1 / (input_sigma^2))*(A'*(A*x_est-y))+beta*(x_est -(v_est-u_est));
        a_res = (A'*A*res) + beta*res;
        mu_opt = mean(res(:).*res(:))/mean(res(:).*a_res(:));
        x_est = x_est + mu_opt*res;
        %x_est = max( min(x_est, 5), 0);

    end


    % relaxation
    x_hat = alpha*x_est + (1-alpha)*v_est;

    %x_hat = x_est;
    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    % using gradient descent
    for j = 1:1:inner_iters2

        %              v_est = reshape(v_est,30,100);
        %           f_v_est = imnlmfilt(v_est,'DegreeOfSmoothing',10,'SearchWindowSize',15);
        %           f_v_est = f_v_est(:);

        v_est = reshape(x_est,m,2*ni);
        v_est1=v_est((1:m),(1:ni));
        v_est2=v_est((1:m),(ni+1:2*ni));
        %            f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
        %             f_v_est2 = medfilt2(v_est2, [median median],'symmetric');
        %
        %             f_v_est1 = BM3D(v_est1, median);
        %             f_v_est2 = BM3D(v_est2, median);


        f_v_est1 = imnlmfilt(v_est1,'DegreeOfSmoothing',median, 'SearchWindowSize', size,'ComparisonWindowSize', comp);
        f_v_est2 = imnlmfilt(v_est2,'DegreeOfSmoothing',median,'SearchWindowSize',size ,'ComparisonWindowSize',comp);


        f_v_est = [f_v_est1 f_v_est2];
        f_v_est = f_v_est(:);



        v_est = (beta*(x_hat + u_est) + lambda*f_v_est)/(lambda + beta);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_hat - v_est;

    r(k) = norm(A*x_hat-y)^2/norm(y)^2;

    if r(k) < error
        out = x_hat;
        %n_out=k;
        res_admm = r;
        break
    end


end
out = x_hat;
% n_out=k;
res_admm = r;
end

function [res_admm,n_out , out ] = admm_red_bm3d(A,y,lambda,error,N,max_iter,inner_iters,inner_iters2,median,m,ni,beta)

x_est= ones([N 1]);
v_est=ones([N 1]);
u_est=ones([N 1]);

%     x_est= 0.5*abs(x_est)+0.75;
%     v_est= 0.5*abs(v_est)+0.75;
%     u_est= 0.5*abs(u_est)+0.75;

n_out = [x_est v_est u_est];
alpha=1.5;


for k = 1:1:max_iter

    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2

    for j = 1:1:inner_iters
        %                e = A'*(A*x_est-y) + beta*(x_est-(v_est-u_est));
        %                 r = A'*A*e + beta*e;
        %                 mu_opt = (e'*e)/(e'*r);
        %                 x_est = x_est + mu_opt*e;
        %                 %x_est = max( min(x_est, 10), 0);


        b = (A*y) + beta*(v_est - u_est);
        A_x_est = (A'*A*x_est) + beta*x_est;
        res = b - A_x_est;
        %res = (1 / (input_sigma^2))*(A'*(A*x_est-y))+beta*(x_est -(v_est-u_est));
        a_res = (A'*A*res) + beta*res;
        mu_opt = mean(res(:).*res(:))/mean(res(:).*a_res(:));
        x_est = x_est + mu_opt*res;
        %x_est = max( min(x_est, 5), 0);

    end


    % relaxation
    x_hat = alpha*x_est + (1-alpha)*v_est;

    %x_hat = x_est;
    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    % using gradient descent
    for j = 1:1:inner_iters2

        %              v_est = reshape(v_est,30,100);
        %           f_v_est = imnlmfilt(v_est,'DegreeOfSmoothing',10,'SearchWindowSize',15);
        %           f_v_est = f_v_est(:);

        v_est = reshape(x_est,m,2*ni);
        v_est1=v_est((1:m),(1:ni));
        v_est2=v_est((1:m),(ni+1:2*ni));
        %            f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
        %             f_v_est2 = medfilt2(v_est2, [median median],'symmetric');
        %
        f_v_est1 = BM3D(v_est1, median);
        f_v_est2 = BM3D(v_est2, median);


        %             f_v_est1 = imnlmfilt(v_est1,'DegreeOfSmoothing',median);
        %             f_v_est2 = imnlmfilt(v_est2,'DegreeOfSmoothing',median);


        f_v_est = [f_v_est1 f_v_est2];
        f_v_est = f_v_est(:);



        v_est = (beta*(x_hat + u_est) + lambda*f_v_est)/(lambda + beta);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_hat - v_est;

    r(k) = norm(A*x_hat-y)^2/norm(y)^2;

    if r(k) < error
        out = x_hat;
        %n_out=k;
        res_admm = r;
        break
    end


end
out = x_hat;
% n_out=k;
res_admm = r;
end