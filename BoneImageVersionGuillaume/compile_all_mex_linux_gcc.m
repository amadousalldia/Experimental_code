
%% function compiling a given mex with 03 optimization level and OpemMP lib linked.
% has to be called where the source code is: 
% comprising main.cpp + functions/*cpp + include/*hpp
% will output a bin.mexa64 that needs to be renamed in your main script
base_dir = '/home_local/SimulBone/BQUI/ReconImage/BoneImageVersionGuillaume/';
cd([base_dir 'AutofocusSoftTissue_21April2021/']);
mex('-setup', 'CPP');
mex -DMEX -Iinclude CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS="$LDFLAGS -fopenmp"...
functions/*.cpp main.cpp -output ../AutofocusTissue

cd([base_dir 'AutofocusBone_21April2021/']);
mex('-setup', 'CPP');
mex -DMEX -Iinclude CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS="$LDFLAGS -fopenmp"...
functions/*.cpp main.cpp -output ../AutofocusBone

cd([base_dir 'ForPostProcessing_FixedFNumber_15Oct2020/']);
mex('-setup', 'CPP');
mex -DMEX -Iinclude CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS="$LDFLAGS -fopenmp"...
functions/*.cpp main.cpp -output ../PostProcessingFNumber_old

cd([base_dir 'ForPostProcessing_FixedFNumber_3May2021/']);
mex('-setup', 'CPP');
mex -DMEX -Iinclude CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS="$LDFLAGS -fopenmp"...
functions/*.cpp main.cpp -output ../PostProcessingFNumber_new

cd(base_dir)