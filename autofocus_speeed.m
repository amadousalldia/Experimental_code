%focus_qual = struct('sharpness_NV',sharpness_NormVariance, 'sharpness_B',sharpness_Brenner, ...
 %   'intensity', spatial_energy_2D, 'TEST_VAL', C_TISSUE_TEST_VALUES);
ACQ_TECHNIC = 'SA';    % 'PW' PWHQ' 'SA'
SECTION     = 'Trans'; % 'Longi' 'Trans'
BONE_TYPE = 'Tibia';   % 'Radius' 'Tibia'
HEAD_POSITION = 'U';  % 'D' 'U'
%
j=1;
L1T= zeros([1 32]);
L2T= zeros([1 32]);
L3T= zeros([1 32]);
L1B= zeros([1 32]);
L2B= zeros([1 32]);
L3B= zeros([1 32]);
for k=5:8
    for acq_seq=1:4
        for repet_numb=1:2
        file_name = sprintf('%s_%s_Sujet4_%s_N%d_H%s_%d_%d',...
            ACQ_TECHNIC, SECTION, ...
            BONE_TYPE, acq_seq, HEAD_POSITION, k, repet_numb);
            FOC_QUAL = load([file_name '_foc_qual_speed.mat']);
            MAP_TISS = FOC_QUAL.MAP_TISS;
            MAP_TISS=MAP_TISS(file_name);
            L1T(j) = MAP_TISS.TEST_VAL(MAP_TISS.intensity==max(MAP_TISS.intensity));
            L2T(j) = MAP_TISS.TEST_VAL(MAP_TISS.sharpness_NV==max(MAP_TISS.sharpness_NV));
            L3T(j) = MAP_TISS.TEST_VAL(MAP_TISS.sharpness_B==max(MAP_TISS.sharpness_B));
            
            MAP_BONE = FOC_QUAL.MAP_BONE;
            MAP_BONE=MAP_BONE(file_name);
            L1B(j) = MAP_BONE.TEST_VAL(MAP_BONE.intensity==max(MAP_BONE.intensity));
            L2B(j) = MAP_BONE.TEST_VAL(MAP_BONE.sharpness_NV==max(MAP_BONE.sharpness_NV));
            L3B(j) = MAP_BONE.TEST_VAL(MAP_BONE.sharpness_B==max(MAP_BONE.sharpness_B));            
            j=j+1;
        end
    end
end

figure, 
subplot 211,
plot(L1T,'kd--'), hold on, plot(L2T, 'b+-.'), 
plot(L3T,'r*:'), hold off
axis tight
xline(8);
xline(16);
xline(24);
ylim([min(MAP_TISS.TEST_VAL) max(MAP_TISS.TEST_VAL)])
xticks([4 12 20 28])
xticklabels([{'Distal(5)'} {'Distal(6)'} {'Proximal(7)'} {'Proximal(8)'}])
title('Optimal speed of sound by diffferent image quality metrics in soft tissues',...
    'Interpreter', 'latex')
subplot 212,
plot(L1B,'kd--'), hold on, plot(L2B, 'b+-.'), 
plot(L3B,'r*:'), hold off
axis tight
xline(8);
xline(16);
xline(24);
ylim([min(MAP_BONE.TEST_VAL) max(MAP_BONE.TEST_VAL)])
xticks([4 12 20 28])
xticklabels([{'Distal(5)'} {'Distal(6)'} {'Proximal(7)'} {'Proximal(8)'}])
title('Optimal speed of sound by diffferent image quality metrics in bone tissues',...
    'Interpreter', 'latex')
