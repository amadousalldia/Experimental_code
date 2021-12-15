
%% function compiling a given mex with 03 optimization level and OpemMP lib linked.
% has to be called where the source code is: 
% comprising main.cpp + functions/*cpp + include/*hpp
% will output a bin.mexa64 that needs to be renamed in your main script

% Compilation on Windows with a Visual C++ compiler

mex('-setup', 'CPP');
mex -Iinclude -DMEX COMPFLAGS="$COMPFLAGS /openmp" ...
functions/*.cpp main.cpp -output PostProcessingFNumber



