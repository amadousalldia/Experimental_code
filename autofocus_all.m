%
%---- Input parameters -----%
% LENS_THICKNESS: Lens thickness (meter)
% C_LENS: known wavespeed in the lens of the probe (m/s)
% FS: Sampling frequency (Hz)
% PITCH: array pitch (meter)
% HalfOpeningAngInLensRad: max half opening angle at receive element (radian)
% C_TISSUE: a priori wavespeed in soft tissue between probe and bone (m/s)
% NCORE: number of CPU cores to use for OMP parallel computing
% (X_El, Z_El): coordinates of array elements (meter)
% (XS, ZS): coordinates of virtual sources (meter)
% add_to_delay_firing: travel time from virtual sources to nearest array element (second)
% (X, Z): coordinates of image pixels (meter)
% (I_SIG, Q_SIG): in-phase and quadrature sampled RF data, 3D array dimension is (Time,Receiver index,Transmission index)
% Tx_to_be_used: indices of transmissions to be used
C_TISSUE = 1500;
C_TISSUE_TEST_VALUES = C_TISSUE-200:15:C_TISSUE+200; % define here the range of values you want to test

% Fill in the Setup array.
Setup_soft = [PROBE_PARAM.LENS_THICKNESS PROBE_PARAM.PITCH ...
    HalfOpeningAngInLensRad PROBE_PARAM.Fs PROBE_PARAM.C_LENS C_TISSUE NCORE];
tic
[ImageTissue, spatial_energy_2D, ...
    sharpness_NormVariance, sharpness_Brenner] = AutofocusTissue(...
    Setup_soft, XR, ZR, XS, ZS, X, Z,...
    I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, C_TISSUE_TEST_VALUES);
toc
%
focus_qual = struct('sharpness_NV',sharpness_NormVariance, 'sharpness_B',sharpness_Brenner, ...
    'intensity', spatial_energy_2D, 'TEST_VAL', C_TISSUE_TEST_VALUES);
file_name = sprintf('%s_%s_Sujet4_%s_N%d_H%s_%d_%d',...
    ACQ_TECHNIC, SECTION, ...
    BONE_TYPE, N, HEAD_POSITION, k, repet_numb);

MAP_TISS = containers.Map(file_name,focus_qual);
% save([file_name '_tiss_speed.mat'],'MAP_TISS');
C_TISSUE = C_TISSUE_TEST_VALUES(sharpness_Brenner==max(sharpness_Brenner));
figure(1)
imagesc(X*1e3,Z*1e3,ImageTissue)
axis image


figure(111)
plot(C_TISSUE_TEST_VALUES,spatial_energy_2D,'kv-')
hold on
plot(C_TISSUE_TEST_VALUES,sharpness_Brenner,'g.-')
plot(C_TISSUE_TEST_VALUES,sharpness_NormVariance,'m.-')
ylabel('focus quality','fontsize',16)
xlabel('wavespeed','fontsize',16)
legend('intensity','sharpness Brenner','sharpness NormVariance')
set(gca,'fontsize',16)
%%

%---- Input parameters -----%
% LENS_THICKNESS: Lens thickness (meter)
% C_LENS: known wavespeed in the lens of the probe (m/s)
% FS: Sampling frequency (Hz)
% PITCH: array pitch (meter)
% D_PROBE_PERIOS: a priori max depth of outer surface of bone (meter)
% MIN_CORTICAL: a priori min thickness of bone
% MAX_CORTICAL: a priori max thickness of bone
% HalfOpeningAngInLensRad: max half opening angle at receive element (radian)
% C_TISSUE: known wavespeed in soft tissue between probe and bone (m/s)
% C_AXIAL: wavespeed in bone in direction of max wavespeed (m/s)
% C_RADIAL: wavespeed in bone in direction of min wavespeed (m/s)
% ANISO_SHAPE_COEF: anisotropy shape factor (unitless)
% NCORE: number of CPU cores to use for OMP parallel computing
% WhichModel: 1, 2, 3 or 4 (see below)
% (X_El, Z_El): coordinates of array elements (meter)
% (XS, ZS): coordinates of virtual sources (meter)
% add_to_delay_firing: travel time from virtual sources to nearest array element (second)
% (X, Z): coordinates of image pixels (meter)
% (I_SIG, Q_SIG): in-phase and quadrature sampled RF data, 3D array dimension is (Time,Receiver index,Transmission index)
% Tx_to_be_used: indices of transmissions to be used


    
%---- Select which parameter to test ----%
%     1 = CAxial     % for anisotropic model
%     2 = CRadial    % for anisotropic model
%     3 = AnisoShape % for anisotropic model
%     4 = Isotropic velocity model
%----------------------------------------%
    C_RADIAL = 3000;
    WhichModel = 4; 
    C_RADIAL_TEST_VALUES = 2900:25:3500; % define here the range of values you want to test

    % Fill in the Setup array.
    Setup_bone= [PROBE_PARAM.LENS_THICKNESS PROBE_PARAM.PITCH...
        D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
        HalfOpeningAngInLensRad PROBE_PARAM.Fs PROBE_PARAM.C_LENS ...
        C_TISSUE C_AXIAL C_RADIAL ANISO_SHAPE_COEF NCORE, WhichModel];
tic
[ImageTissue, ImageBone, spatial_energy_2D, ...
    sharpness_NormVariance, sharpness_Brenner, PeriParab] = ...
    AutofocusBone(Setup_bone, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, C_RADIAL_TEST_VALUES);
toc
figure(2)
subplot 121
imagesc(X*1e3,Z*1e3,ImageTissue)
axis image
subplot 122
imagesc(X*1e3,Z*1e3,ImageBone)
axis image

figure(222), hold on
plot(C_RADIAL_TEST_VALUES,spatial_energy_2D, 'vk-','DisplayName','Intensity')
plot(C_RADIAL_TEST_VALUES,sharpness_Brenner, 'g.-','DisplayName','Brenner')
plot(C_RADIAL_TEST_VALUES, sharpness_NormVariance, 'm.-','DisplayName','NormVar')
ylabel('focus quality','fontsize',16)
xlabel('wavespeed','fontsize',16)
legend
axis tight

best=[1450, 3260];
focus_qual = struct('sharpness_NV',sharpness_NormVariance, 'sharpness_B',sharpness_Brenner, ...
    'intensity', spatial_energy_2D, 'TEST_VAL', C_RADIAL_TEST_VALUES);
file_name = sprintf('%s_%s_Sujet4_%s_N%d_H%s_%d_%d',...
    ACQ_TECHNIC, SECTION, ...
    BONE_TYPE, N, HEAD_POSITION, k, repet_numb);

MAP_BONE = containers.Map(file_name,focus_qual);

save([file_name '_foc_qual_speed.mat'],'MAP_BONE', 'MAP_TISS'); 
        