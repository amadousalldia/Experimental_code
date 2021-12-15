%--- IMPORTANT ---%
% ALL PARAMETERS MUST BE OF CLASS REAL DOUBLE
% VECTORS MUST BE LINE VECTORS


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


C_TISSUE_TEST_VALUES = c_tissue-30:5:c_tissue+30; % define here the range of values you want to test

% Fill in the Setup array.
Setup = [LENS_THICKNESS PITCH ...
    HalfOpeningAngInLensRad FS C_LENS C_TISSUE NCORE];

[ImageTissue, spatial_energy_2D, sharpness_NormVariance, sharpness_Brenner] = AutofocusTissue(Setup, X_El, Z_El, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, C_TISSUE_TEST_VALUES);

%%

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



