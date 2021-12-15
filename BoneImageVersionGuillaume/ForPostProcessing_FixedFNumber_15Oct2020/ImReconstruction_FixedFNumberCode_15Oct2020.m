clear all
close all
clc

%----------
%----------


filtering_RF_data_YES = 1;
apodization_start_end_YES = 1;
TGC_YES = 1;
correction_delay_matching_layers_YES = 1;

NCORE                   = 1; % number of CPU cores for OpenMP parallel computing
ReConTo                          = 2; % 1 for soft tissue only, 2 for cutaneous tissue and cortical bone, 3 for until the marrow
NeedTravelTime                   = 0; % 0, 1, 2 for compute without saving, compute travel times/receive & transmit angles and save them, do not calculate but load travel times/receive & transmit angles.

%--------------------------------------------------%
%           P4-1 probe
%--------------------------------------------------%
% lens parameters found with semblance on single reflection on interface
% silicone/water (15 oct 2019)
correction_delay    = -70.9e-9; % [s]
offset              = 9; % time to peak of the round-trip waveform [samples]
Fs                  = 10e6; % [Hz]
tone_burst_cycles   = 4;
C_LENS              = 970; % [m/s]
LENS_THICKNESS      = 1.3e-3; % [m]
NELEMENTS           = size(SIG,2);
PITCH               = 0.295e-3; % [m]
FREQ_Transducteur   = 2.5e6; % [Hz]

% horizontal positions of the elements
XR = (0:NELEMENTS-1) *PITCH; % [m]
XR = XR - mean(XR); % place the center of coordinate system in the middle of transducer array
% vertical positions of the elements
ZR = zeros(1,length(XR)); % [m]
% coordinates of sources
XS = XR;
ZS = ZR;
NSOURCES = size(SIG,3);
% for synthetic aperture imaging, zero transmit delay added (only used for virtual sources)
add_to_delay_firing = zeros(1,NSOURCES);

% velocity model
C_TISSUE            = 1570; % [m/s]
C_AXIAL             = 3200;%3950; % [m/s]
C_RADIAL            = 3200; % [m/s]
ANISO_SHAPE_COEF    = 1.5; % unitless
C_MARROW            = 1400; % [m/s]

% a priori knownledge on anatomy
D_PROBE_PERIOS       = 20e-3;
MIN_CORTICAL         = 4e-3;
MAX_CORTICAL         = 12e-3;

% define dimension and pixel size of reconstructed image
xmin    = -8e-3; % [m]
xmax    = 8e-3; % [m]
zmin    = 1.5e-3; % [m]
zmax    = 35e-3; % [m]
resolution = C_TISSUE/(FREQ_Transducteur*2); % [m]

% create image axes/pixels coordinates
X       = xmin:resolution:xmax; % image width [m]
Z       = zmin:resolution:zmax; % image depth [m]

% Sources to be used during reconstruction
Tx_to_be_used       = 1:NSOURCES;

% parameter only for ImageTissue, ImageBone, ImageMarrow:
HalfOpeningAngInSkinDeg = 40; % half opening angle in receive only for image reconstructed with full aperture in receive [degree]
HalfOpeningAngInLensRad = asin(C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]

% Remark: no acceptance angle HalfOpeningAng used in transmit

% parameters only for Beam_I_... and Beam_Q_... images:
FnumberMin                       = 0.5; % fixed F-number, for large depth actual F-number might be larger
SubApertureApodis                = 2; % 0, 1, 2 for no windowing, hamming, hann.
% Receive angle(s)
ReceiveAngle    = deg2rad([0]); % receive angles, can be a vector for vector flow imaging [rad]
Max_err_allowed_ReceiveAngle_rad = deg2rad(5); % [rad]

% Parabolas (if needed the user can force the parabolas of periosteum and/or endosteum)
PeriParabIn = [ 0 0 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros
EndoParabIn = [ 0 0 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros


%%

% backward shift of raw echo signals by time to peak of the round-trip waveform
SIG = SIG(offset:end,:,:);

Nt = size(SIG,1); % number of temporal samples

if apodization_start_end_YES
    muting_duration                     = 3e-6;
    SIG(1:ceil(muting_duration*Fs),:,:) = 0;
end

if filtering_RF_data_YES
b = fir1(90,[0.2 1.8]*FREQ_Transducteur/(Fs/2),'bandpass');
    for ii=1:size(SIG,2)
        for jj=1:size(SIG,3)
        SIG(:,ii,jj) = filtfilt(b,1,SIG(:,ii,jj));
        end
    end
end

if correction_delay_matching_layers_YES
    nb_pts_fft = 2^nextpow2(Nt);
    frek = [[0:nb_pts_fft/2-1] [-nb_pts_fft/2:-1]]*Fs/nb_pts_fft;
    frek = repmat(frek.',1,NELEMENTS);
    delay_time = correction_delay*ones(1,NELEMENTS);
    for iTx=1:NSOURCES
    FFT_matrix = fft(squeeze(SIG(:,:,iTx)),nb_pts_fft);
    sig_delayed = pulse_delayingXFFT(FFT_matrix,delay_time,frek);
    SIG(:,:,iTx) = sig_delayed(1:Nt,:);
    end
end

if TGC_YES
    %%% Time Gain Compensation
    TGC = repmat((linspace(0,size(SIG,1)-1,size(SIG,1)).').^1.9,1,NELEMENTS,NSOURCES);
    TGC = TGC/max(TGC(:));
    SIG = SIG .* TGC;
end

Time = [0:Nt-1]/Fs;
displayed_dynamic_range_dB = 50;
figure(1)
subplot 131
TX_Element = 1;
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
subplot 132
TX_Element = round(NSOURCES/2);
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
subplot 133
TX_Element = NSOURCES;
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
pause(0.5)



%%

% I/Q separation.
I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
for Tx = 1:size(SIG,3)
    tmp =  hilbert(squeeze(SIG(:,:,Tx)));
    I_SIG(:,:,Tx) = real(tmp);
    Q_SIG(:,:,Tx) = imag(tmp);
end



%%


Setup = [ANISO_SHAPE_COEF LENS_THICKNESS PITCH...
    D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL HalfOpeningAngInLensRad FnumberMin...
    Max_err_allowed_ReceiveAngle_rad Fs C_LENS C_TISSUE C_MARROW...
    C_AXIAL C_RADIAL NCORE ReConTo NeedTravelTime SubApertureApodis];

if ReConTo == 1
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue] =...
    PostProcessingFNumber(...
            Setup, XR, ZR, XS, ZS, X, Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
            ReceiveAngle, PeriParabIn, EndoParabIn);
        
    % sum over transmissions for image with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));

elseif ReConTo == 2
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
    APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
    PostProcessingFNumber(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngle, PeriParabIn, EndoParabIn);

    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);

    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
    ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));

elseif ReConTo == 3
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone, APix_T_Bone, APix_R_Bone,...
    ImageMarrow, Time_T_Marrow, Time_R_Marrow, Angle_T_Marrow,...
    Angle_R_Marrow, Beam_I_Marrow, Beam_Q_Marrow, APix_T_Marrow, APix_R_Marrow, PeriParab, EndoParab] =...
    PostProcessingFNumber(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngle, PeriParabIn, EndoParabIn);
    
    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);
    
    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
    ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));
    ImageMarrow_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Marrow,4).^2+sum(Beam_Q_Marrow,4).^2));

end
    


%%    

%---- calculate thickness of cortical bone ----%

figure(33)
plot(X*1e3,fit_curve_Periosteum*1e3,'k','linewidth',2)
hold on
plot(X*1e3,fit_curve_Endosteum*1e3,'k','linewidth',2)
axis ij
xlabel('[mm]')
ylabel('[mm]')
axis equal

P_midLine = mean([P_Periosteum; P_Endosteum]);
fit_curve_midLine = polyval(P_midLine,X);
plot(X*1e3,fit_curve_midLine*1e3,'b--','linewidth',2)

x1 = X;
fit_Periosteum = fit_curve_Periosteum;
for ii=[1 length(x1)]
    pt_Periosteum = [x1(ii), fit_Periosteum(ii)];
    N = [2* P_Periosteum(1)*pt_Periosteum(1) + P_Periosteum(2), -1];
    N = -N / norm(N);
    
    a = P_Endosteum(1);
    b = P_Endosteum(2);
    c = P_Endosteum(3);
    thickness1 = intersection_rayon_interface(a,b,c,pt_Periosteum,N);
    pt_Endosteum = pt_Periosteum + thickness1*N;
    if ii==1
        pt_Endosteum_left = pt_Endosteum;
    elseif ii==length(x1)
        pt_Endosteum_right = pt_Endosteum;
    end
end
X_start_Endosteum = pt_Endosteum_left(1);
X_end_Endosteum = pt_Endosteum_right(1);


x1 = X;
thicknessB = [];
thicknessT = [];
for ii=1:length(x1)
    pt_midLine = [x1(ii), fit_curve_midLine(ii)];
    N = [2* P_midLine(1)*pt_midLine(1) + P_midLine(2), -1];
    N = -N / norm(N);
    
    a = P_Endosteum(1);
    b = P_Endosteum(2);
    c = P_Endosteum(3);
    thickness = intersection_rayon_interface(a,b,c,pt_midLine,N);
    pt_Endosteum = pt_midLine + thickness*N;
    
    if (pt_Endosteum(1)>X_start_Endosteum) && (pt_Endosteum(1)<X_end_Endosteum)
        thicknessB = [thicknessB thickness];
        figure(33)
        plot([pt_midLine(1) pt_Endosteum(1)]*1e3,[pt_midLine(2) pt_Endosteum(2)]*1e3,'r')

        N = -N;
        a = P_Periosteum(1);
        b = P_Periosteum(2);
        c = P_Periosteum(3);
        thickness = intersection_rayon_interface(a,b,c,pt_midLine,N);
        thicknessT = [thicknessT thickness];
        pt_Periosteum = pt_midLine + thickness*N;
        figure(33)
        plot([pt_midLine(1) pt_Periosteum(1)]*1e3,[pt_midLine(2) pt_Periosteum(2)]*1e3,'m')
    end
end

thicknessMidLine = thicknessT+thicknessB;
mean_cortical_thickness = mean(abs(thicknessMidLine));
title(['mean cortical thickness = ' num2str(mean_cortical_thickness*1e3,3) ' mm']);
disp(['mean cortical thickness (normal to mean axis) is ' num2str(mean_cortical_thickness*1e3,3) ' mm']);

%%    

% for iTx=1:NSOURCES
%     
%     figure(1)
%     imagesc(X*1e3,Z*1e3,20*log10(sqrt(squeeze(Beam_I_Tissue(:,:,1,iTx)).^2+squeeze(Beam_Q_Tissue(:,:,1,iTx)).^2)))
%     hold on
%     plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
%     hold off
%     axis image
% 
%     pause
% end

if ReConTo == 1

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

elseif ReConTo == 2

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

    figure(3)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image

    FullImage = ImageTissue;
    FullImage(ImageBone~=0) = ImageBone(ImageBone~=0);
    max_depth_index_FullImage = PeriParab(5) + round(MAX_CORTICAL/resolution);
    FullImage  = FullImage(1:max_depth_index_FullImage,:);

    FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
    FullImage_FixedReceiveAng(ImageBone_FixedReceiveAng~=0) = ImageBone_FixedReceiveAng(ImageBone_FixedReceiveAng~=0);
    FullImage_FixedReceiveAng  = FullImage_FixedReceiveAng(1:max_depth_index_FullImage,:);
    
    
    figure(10)
    subplot 121
    imagesc(X*1e3,Z(1:max_depth_index_FullImage)*1e3,20*log10(FullImage))
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    subplot 122
    imagesc(X*1e3,Z(1:max_depth_index_FullImage)*1e3,20*log10(FullImage_FixedReceiveAng))
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    
elseif ReConTo == 3

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

    figure(3)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image

    figure(4)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageMarrow))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageMarrow_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    
    FullImage = ImageTissue;
    FullImage(ImageBone~=0) = ImageBone(ImageBone~=0);
    FullImage(ImageMarrow~=0) = ImageMarrow(ImageMarrow~=0);
    
    FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
    FullImage_FixedReceiveAng(ImageBone_FixedReceiveAng~=0) = ImageBone_FixedReceiveAng(ImageBone_FixedReceiveAng~=0);
    FullImage_FixedReceiveAng(ImageMarrow_FixedReceiveAng~=0) = ImageMarrow_FixedReceiveAng(ImageMarrow_FixedReceiveAng~=0);

    figure(10)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(FullImage))
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(FullImage_FixedReceiveAng))
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    
end





