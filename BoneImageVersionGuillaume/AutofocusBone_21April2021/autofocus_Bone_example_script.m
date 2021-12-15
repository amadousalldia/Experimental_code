%--- IMPORTANT ---%
% ALL PARAMETERS MUST BE OF CLASS REAL DOUBLE
% VECTORS MUST BE LINE VECTORS


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
%     1 = CAxial
%     2 = CRadial
%     3 = AnisoShape
%     4 = Isotropic velocity model
%----------------------------------------%
    WhichModel = 4; 
    model_Values = C_RADIAL-150:25:C_RADIAL+150; % define here the range of values you want to test
    
    % Fill in the Setup array.
    Setup = [LENS_THICKNESS PITCH...
        D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
        HalfOpeningAngInLensRad FS C_LENS C_TISSUE C_AXIAL C_RADIAL ANISO_SHAPE_COEF...
        NCORE, WhichModel];

    [ImageTissue, ImageBone, spatial_energy_2D, sharpness_NormVariance, sharpness_Brenner, PeriParab] = AutofocusBone(Setup, X_El, Z_El, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, model_Values);
        

    
%%

figure(1)
subplot 121
imagesc(X*1e3,Z*1e3,ImageTissue)
hold on
interf = PeriParab(1)*X.^2+PeriParab(2)*X+PeriParab(3);
plot(X*1e3,interf*1e3,'r')
hold off
axis image
subplot 122
imagesc(X*1e3,Z*1e3,ImageBone)
hold on
plot(X*1e3,interf*1e3,'r')
hold off
axis image

figure(111)
subplot 121
plot(model_Values,spatial_energy_2D,'kv-')
hold on
plot(model_Values,sharpness_Brenner,'g.-')
plot(model_Values,sharpness_NormVariance,'m.-')
ylabel('focus quality','fontsize',16)
xlabel('wavespeed','fontsize',16)
legend('intensity','sharpness Brenner','sharpness NormVariance')
set(gca,'fontsize',16)
subplot 122
plot(model_Values,spatial_energy_2D+sharpness_Brenner+sharpness_NormVariance,'bo-')
ylabel('sum of metrics','fontsize',16)
xlabel('wavespeed','fontsize',16)
set(gca,'fontsize',16)



