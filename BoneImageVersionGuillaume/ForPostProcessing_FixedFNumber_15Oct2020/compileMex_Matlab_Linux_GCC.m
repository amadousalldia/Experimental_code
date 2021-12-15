
%% function compiling a given mex with 03 optimization level and OpemMP lib linked.
% has to be called where the source code is: 
% comprising main.cpp + functions/*cpp + include/*hpp
% will output a bin.mexa64 that needs to be renamed in your main script

% Compilation with a gcc compiler (MinGW for Windows)

mex('-setup', 'CPP');
mex -DMEX -Iinclude CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS="$LDFLAGS -fopenmp"...
functions/*.cpp main.cpp -output PostProcessingFNumber

