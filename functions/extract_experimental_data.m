function [SIG, PROBE_PARAM] = extract_experimental_data(subject_name, ACQ_TECHNIC, SECTION, ...
    BONE_TYPE, HEAD_POSITION, N, k, repet_numb)
    % Inputs
%     ACQ_TECHNIC = 'SA';    % 'PW' 'PWHQ' 'SA'
%     SECTION     = 'Trans'; % 'Longi' 'Trans'
%     BONE_TYPE = 'Tibia';   % 'Radius' 'Tibia'
%     HEAD_POSITION = 'U';  % 'D' 'U'
%     N = 1;% sequence number [1 2 3 4]
%     k = 5;% distance [5 6] for distal and  [7 8] for proximal
%     repet_numb = 1;%1 2

    experi_dir = '/home_local/Experimental/etude_clinique_paris/';
    sujet_dir = [experi_dir subject_name];
    subject_name = upper(subject_name);
    subject_name1 = lower(subject_name);
    subject_name1 = [subject_name(1) subject_name1(2:end)];
    rf_name = sprintf('%s_%s_%s_%s_N%d_H%s_%d_%d.mat',...
                    ACQ_TECHNIC,SECTION,subject_name1,BONE_TYPE,...
                    N,HEAD_POSITION,k,repet_numb);

    data_path = sprintf('%s/%s/%s',sujet_dir,'Matlab RF',rf_name);

    experimental_data = load(data_path);

    %
    SIG = experimental_data.rf;
    rf_data  = SIG;
    Fs = 4*experimental_data.param_val.FREQ_Transducteur*1e6;
    %
    %
    %--------------------------------------------------%
    %           P4-1 probe
    %--------------------------------------------------%
    % lens parameters found with semblance on single reflection on interface
    % silicone/water (15 oct 2019)
%     tone_burst_cycles   = 3;
    %offset=1;
    C_LENS              = experimental_data.param_val.C_LENS; % [m/s]
    LENS_THICKNESS      = experimental_data.param_val.LENS_THICKNESS;%1.3e-3; % [m]
    NELEMENTS           = experimental_data.param_val.NELEMENTS;
    PITCH               = experimental_data.param_val.PITCH;%0.295e-3; % [m]
    FREQ_Transducteur   = experimental_data.param_val.FREQ_Transducteur*1e6; % [Hz]
    offset              = experimental_data.param_val.offset; % time to peak of the round-trip waveform [samples]

    PROBE_PARAM = struct('C_LENS',C_LENS, 'LENS_THICKNESS',LENS_THICKNESS,...
        'NELEMENTS', NELEMENTS, 'PITCH', PITCH,'FREQ_Transducteur', FREQ_Transducteur,...
        'offset', offset, 'Fs', Fs);
end

