% This code outputs two types of images:
% 1) envelope images using the full aperture in receive, provided the receive ray angle at the element is smaller than an acceptance angle
% 2) I/Q images for vector flow imaging (or B-mode imaging), with fixed F-number and receive angle at the pixel

% If multiple transmissions and/or multiple receive angles are used, I/Q images are four-dimensional arrays (Z,X,Rx,Tx)


% The code also outputs:
% - the transmit and receive travel times between sources and pixels and, pixels and array elements (Time_T_..., Time_R_...)
% - the transmit and receive ray angles at sources and at the array elements (Angle_T_...,Angle_R_...)
% - the transmit and receive ray angles at the pixel (APix_T_..., APix_R_...)
% - the coefficients of the parabola (z = a*x^2+b*x+c) for segmentation of the outer surface of the bone layer: PeriParab = [a b c IdxMinDepthSegmentation IdxMaxDepthSegmentation]
% - the coefficients of the parabola (z = a*x^2+b*x+c) for segmentation of the inner surface of the bone layer: EndoParab = [a b c IdxMinDepthSegmentation IdxMaxDepthSegmentation]


%--- IMPORTANT ---%
% ALL PARAMETERS MUST BE OF CLASS REAL DOUBLE
% VECTORS MUST BE LINE VECTORS


%---- Input parameters -----%
% LENS_THICKNESS: Lens thickness (meter)
% C_LENS: wavespeed in the lens of the probe (m/s)
% FS: Sampling frequency (Hz)
% PITCH: array pitch (meter)
% D_PROBE_PERIOS: a priori max depth of outer surface of bone (meter)
% MIN_CORTICAL: a priori min thickness of bone
% MAX_CORTICAL: a priori max thickness of bone
% C_TISSUE: wavespeed in soft tissue between probe and bone (m/s)
% C_AXIAL: wavespeed in bone in direction of max wavespeed (m/s)
% C_RADIAL: wavespeed in bone in direction of min wavespeed (m/s)
% ANISO_SHAPE_COEF: anisotropy shape factor (unitless)
% C_MARROW: wavespeed in the marrow/brain (m/s)
% NCORE: number of CPU cores to use for OpenMP parallel computing
% (X_El, Z_El): coordinates of array elements (meter)
% (XS, ZS): coordinates of virtual sources (meter)
% add_to_delay_firing: travel time from virtual sources to nearest array element (second)
% (X, Z): coordinates of image pixels (meter)
% (I_SIG, Q_SIG): in-phase and quadrature sampled RF data, 3D array dimension is (Time,Receiver index,Transmission index)
% Tx_to_be_used: indices of transmissions to be used
% ReConTo: flag that indicates how many layers should be modeled, 1 for uniform soft tissue only, 2 for cutaneous tissue and cortical bone, 3 for until the marrow/brain
% NeedTravelTime: flag = 0, 1, 2 for compute without saving, compute travel times/receive & transmit angles and save them, do not calculate but load travel times/receive & transmit angles.

% Parabolas (if needed the user can force the parabolas of periosteum and/or endosteum)
% PeriParabIn = [ 0 0 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros
% EndoParabIn = [ 0 0 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros



%---- Parameter only for envelope images reconstructed with full aperture in receive (ImageTissue, ImageBone, ImageMarrow):
% HalfOpeningAngInLensRad: max half opening angle at receive element (radian)
% HalfOpeningAngInSkinDeg = 45; % half opening angle in receive only for image reconstructed with full aperture in receive [degree]
% HalfOpeningAngInLensRad = asin(C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]

% Remark: no acceptance angle HalfOpeningAng used in transmit



%---- Parameters only for I/Q images (Beam_I_... and Beam_Q_... images) used for flow analysis:
% FnumberMin: fixed receive F-number, for large depth F-number might be larger if aperutre is not large enough to ensure the nominal F-number
% SubApertureApodis: receive apodization --> 0, 1, 2 for no apodization, hamming, hann.

% ReceiveAngleRad: receive angles at the pixel, typically 0 but can be a vector for vector flow imaging [rad]
% MaxErrorReceiveAngleRad: tolerance on receive angle when selecting the center element of the receive aperture (sub-group of elements) [rad]



%%


% Fill in the Setup array.
Setup = [ANISO_SHAPE_COEF LENS_THICKNESS PITCH...
    D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL HalfOpeningAngInLensRad FnumberMin...
    MaxErrorReceiveAngleRad FS C_LENS C_TISSUE C_MARROW...
    C_AXIAL C_RADIAL NCORE ReConTo NeedTravelTime SubApertureApodis];

if ReConTo == 1
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue] =...
    PostProcessingFNumber(...
            Setup, X_El, Z_El, XS, ZS, X, Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
            ReceiveAngleRad, PeriParabIn, EndoParabIn);
        
    % Create a B-mode image from I/Q images:
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

    % Create a B-mode image from I/Q images:
    % sum over transmissions for image with fixed receive angle and F-number
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
    
    % Create a B-mode image from I/Q images:
    % sum over transmissions for image with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
    ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));
    ImageMarrow_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Marrow,4).^2+sum(Beam_Q_Marrow,4).^2));

end
    


%%    

% Display reconstructed image

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
    max_depth_index_FullImage = PeriParab(5) + round(MAX_CORTICAL/pixel_size);
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





