addpath /home_local/SimulBone/BQUI/ReconImage/Beamforming/
clearvars -except N k repet_numb
clc, %close all

% Adding utilities functions
addpath /home_local/Experimental/functions/
addpath /home_local/SimulBone/BQUI/ReconImage/BoneImageVersionGuillaume/
addpath /home_local/SimulBone/BQUI/ReconImage/Beamforming/functions/
addpath /home_local/SimulBone/BQUI/Simulation/utils/

% Boolean for parameters 
%
filtering_RF_data_YES = 0;
apodization_start_end_YES = 0;
TGC_YES = 0;
correction_delay_matching_layers_YES=0;
mean_th_YES=0;

FnumberMin = 1.90; % fixed F-number,

experi_dir = '/home_local/Experimental/etude_clinique_paris/';
sujet_dir = [experi_dir 'SUJET4'];



ACQ_TECHNIC = 'SA';    % 'PW' PWHQ' 'SA'
SECTION     = 'Trans'; % 'Longi' 'Trans'
BONE_TYPE = 'Tibia';   % 'Radius' 'Tibia'
HEAD_POSITION = 'U';  % 'D' 'U'
% N = 1;% sequence number [1 2 3 4]
% k = 5;% distance [5 6] for distal and  [7 8] for proximal
% repet_numb = 1;%1 2
[SIG, PROBE_PARAM] = extract_experimental_data('SUJET4', ACQ_TECHNIC, SECTION, ...
    BONE_TYPE, HEAD_POSITION, N, k, repet_numb);

%
% speed_data_name = sprintf('%s/%s/all_bone_sound_speed_%s.mat',sujet_dir,abs_dir,abs_type);
%
% all_speed = load(speed_data_name);
% all_speed = all_speed.all_sound_speed;
%
rf_data  = SIG;
NCORE                   = 24; % number of CPU cores for OpenMP parallel computing
ReConTo                 = 2; % 1 for soft tissue only, 2 for cutaneous tissue and cortical bone, 3 for until the marrow
NeedTravelTime          = 0; % 0, 1, 2 for compute without saving, compute travel times/receive & transmit angles and save them, do not calculate but load travel times/receive & transmit angles.

%--------------------------------------------------%
%           P4-1 probe
%--------------------------------------------------%
correction_delay    = -70.9e-9; % [s]


NSOURCES = size(SIG,3);
% for synthetic aperture imaging, zero transmit delay added (only used for virtual sources)
add_to_delay_firing = zeros(1,NSOURCES);

% velocity model
C_TISSUE            = 1450;%layers_o(1).cp*1e3; % [m/s]
C_AXIAL             = 3260;%all_speed(p(1:5));% 3496;% [m/s]
C_RADIAL            = C_AXIAL; % [m/s]
ANISO_SHAPE_COEF    = 0;%1.5; % unitless
C_MARROW            = C_TISSUE;%layers_o(3).cp*1e3; % [m/s]

% a priori knownledge on anatomy
D_PROBE_PERIOS       = 8e-3; % 
MIN_CORTICAL         = 2e-3; % Minimum cortical thickness (thickness lower will not be considered)
MAX_CORTICAL         = 8e-3; % Maximum cortical thickness (thickness higher will not be considered) 
% horizontal positions of the elements
XR = (0:PROBE_PARAM.NELEMENTS-1) *PROBE_PARAM.PITCH; % [m]
XR = XR - mean(XR); % place the center of coordinate system in the middle of transducer array
% define dimension and pixel size of reconstructed image
xmin    = -8e-3;%XR(1)-PROBE_PARAM.PITCH;%-10e-3;XR(1);%- PITCH;% [m]
xmax    = -xmin;% [m]
zmin    = PROBE_PARAM.LENS_THICKNESS;%C_TISSUE/FREQ_Transducteur*2; % [m]
z1 = 15e-3;
zmax    = max(z1,D_PROBE_PERIOS+MAX_CORTICAL); % [m]
resolution = C_TISSUE/(PROBE_PARAM.FREQ_Transducteur*4); % [m]
%
% create image axes/pixels coordinates
X       = xmin:resolution:xmax; % image width [m]
Z       = zmin:resolution:zmax; % image depth [m]

%
% Sources to be used during reconstruction
Tx_to_be_used       = 1:NSOURCES;
%

% parameter only for ImageTissue, ImageBone, ImageMarrow:
HalfOpeningAngInSkinDeg = 40; % half opening angle in receive only for image reconstructed with full aperture in receive [degree]
HalfOpeningAngInLensRad = asin(PROBE_PARAM.C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]

% Remark: no acceptance angle HalfOpeningAng used in transmit

% parameters only for Beam_I_... and Beam_Q_... images:
SubApertureApodis                = 2; % 0, 1, 2 for no windowing, hamming, hann.

% Receive angle(s)
ReceiveAngle    = deg2rad([0]); % receive angles, can be a vector for vector flow imaging [rad]
Max_err_allowed_ReceiveAngle_rad = deg2rad(10); % [rad]

% Parabolas (if needed the user can force the parabolas of periosteum and/or endosteum)
PeriParabIn = [ 0 0 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros
EndoParabIn = [ 0 0 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros

% vertical positions of the elements
ZR = zeros(1,length(XR));% [m]
% coordinates of sources
XS = XR;
ZS = ZR;
%
% PROBE_PARAM.offset = 8;
% backward shift of raw echo signals by time to peak of the round-trip waveform
SIG = SIG(PROBE_PARAM.offset+1:end,:,:);

Nt = size(SIG,1); % number of temporal samples

if apodization_start_end_YES
    muting_duration                     = 1.2e-6;
    SIG(1:ceil(muting_duration*PROBE_PARAM.Fs),:,:) = 0;
end
%
if filtering_RF_data_YES
    b = fir1(90,[0.2 1.8]*PROBE_PARAM.FREQ_Transducteur/(PROBE_PARAM.Fs/2),'bandpass');
    for ii=1:size(SIG,2)
        for jj=1:size(SIG,3)
        SIG(:,ii,jj) = filtfilt(b,1,SIG(:,ii,jj));
        end
    end
end


Time = (0:Nt-1)/PROBE_PARAM.Fs;

%

% I/Q separation.
I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
for Tx = 1:size(SIG,3)
    tmp =  hilbert(squeeze(SIG(:,:,Tx)));
    I_SIG(:,:,Tx) = real(tmp);
    Q_SIG(:,:,Tx) = imag(tmp);
end


Setup = [ANISO_SHAPE_COEF PROBE_PARAM.LENS_THICKNESS PROBE_PARAM.PITCH...
    D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL HalfOpeningAngInLensRad FnumberMin...
    Max_err_allowed_ReceiveAngle_rad PROBE_PARAM.Fs PROBE_PARAM.C_LENS C_TISSUE C_MARROW...
    C_AXIAL C_RADIAL NCORE ReConTo NeedTravelTime SubApertureApodis];
%
addpath /home_local/Experimental/
autofocus_all
return
%%
tic
[ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
PostProcessingFNumber_new(...
    Setup, XR, ZR, XS, ZS, X, Z,...
    I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
    ReceiveAngle, PeriParabIn, EndoParabIn);
toc
%%
display_dynamic_range = 40;

P_Periosteum = PeriParab(1:3);
P_Endosteum = EndoParab(1:3);
fit_curve_Periosteum = polyval(PeriParab(1:3),X);
fit_curve_Endosteum = polyval(EndoParab(1:3),X);

% sum over transmissions for images with fixed receive angle and F-number
% if multiple receive angles, then receive angle must be chose with
% dimension 3 in Beam_I and Beam_Q
ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));

max_depth_index_FullImage = PeriParab(5) + round(MAX_CORTICAL/resolution)-1;
% FullImage  = FullImage(1:max_depth_index_FullImage,:);
FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
FullImage_FixedReceiveAng(ImageBone_FixedReceiveAng~=0) = ImageBone_FixedReceiveAng(ImageBone_FixedReceiveAng~=0);
IQ_bf = FullImage_FixedReceiveAng;
FullImage_FixedReceiveAng  = FullImage_FixedReceiveAng(1:max_depth_index_FullImage,:)/max(FullImage_FixedReceiveAng(:));
figure
imagesc(X*1e3,Z(1:max_depth_index_FullImage)*1e3,20*log10(FullImage_FixedReceiveAng))%/max(FullImage_FixedReceiveAng(:))))%))
hold on
plot(X*1e3,fit_curve_Endosteum*1e3, 'k.')
plot(X*1e3,fit_curve_Periosteum*1e3,'k.')
colormap gray
colorbar
caxis([-display_dynamic_range 0])
if k==5 || k==6
    title('Distal Tibia')
elseif k==7 || k==8
    title('Proximal Tibia')
end


%%
% save the signal and the corresponding parameters
save_dir = [sujet_dir '/BoneImage_output_data/'];
save_name = sprintf('output_%s',...
    rf_name);


save([save_dir save_name],'rf_data','XR','XS','X','Z','ZR','ZS','FS', ...
    'Time_R_Bone', 'Time_R_Tissue', 'Time_T_Bone', 'Time_T_Tissue',...
    'C_AXIAL', 'C_RADIAL', 'C_TISSUE','FREQ_Transducteur',...
    'fit_curve_Periosteum','fit_curve_Endosteum');

return