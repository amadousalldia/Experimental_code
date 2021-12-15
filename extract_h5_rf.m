clc, clearvars, close all
show=0;save_bool=1;
sujet_dir = '/home_local/Experimental/etude_clinique_paris/SUJET4/';

technicAcq = 'SA'; % 'PW' PWHQ'
sectionAcq = 'Trans';% 'Longi'
boneType = 'Tibia'; % 'Radius'
% acq_seq = 1; % [1 2 3 4]
% HEAD_POSITION = 'D';% [D U]
% k=5;% [5 6 7 8] % 5,6 - distal   and 7,8 proximal
% formatstr = '%s_%s_Sujet4_%s_N%d_H%s_%d';
for acq_seq=1:4
    for HEAD_POSITION=['U' 'D']
        for k=5:8
            file_name = sprintf('%s_%s_Sujet4_%s_N%d_H%s_%d.h5',...
                technicAcq,sectionAcq,boneType,acq_seq,HEAD_POSITION,k);
            save_name = sprintf('%s/%s/%s_%s_Sujet4_%s_N%d_H%s_%d',...
                sujet_dir,'Matlab RF', technicAcq,sectionAcq,boneType,acq_seq,HEAD_POSITION,k);
            
            % file_name =  'SA_Trans_Sujet4_Tibia_N2_HD_5.h5';
            % h5disp([sujet_dir file_name]);
            rfdataval = h5read([sujet_dir file_name], '/RF');
            trans_val = h5read([sujet_dir file_name], '/Trans');
            param_val = h5read([sujet_dir file_name], '/param');
            lambda=trans_val.spacingMm/trans_val.spacing;
            XS = (0:param_val.NELEMENTS-1)*param_val.PITCH;
            ZS = zeros(size(XS));
            if strcmp(technicAcq, 'SA')            
                param_val.XS = (0:param_val.NELEMENTS-1)*param_val.PITCH;
                param_val.ZS = zeros(size(XS));
                param_val.NSOURCES = param_val.NELEMENTS;
                param_val.offset = 8;
            end
            connector = trans_val.Connector;
            rf1=reshape(rfdataval(:,:,1),[384, 96, 128]);
            rf1 = rf1(:,:,connector);
            rf2=reshape(rfdataval(:,:,2),[384, 96, 128]);
            rf2 = rf2(:,:,connector);
            if show
                tx=48;
                figure, imagesc(20*log10(abs(hilbert(rf2(:,:,tx)))))
                colorbar,colormap gray
                figure, imagesc(20*log10(abs(hilbert(rf1(:,:,tx)))))
                colorbar,colormap gray
            end
            if save_bool
                rf=rf1;
                save([save_name '_1.mat'],'rf', 'param_val')
                rf=rf2;
                save([save_name '_2.mat'],'rf', 'param_val')
            end
            
        end
    end
end
%%
