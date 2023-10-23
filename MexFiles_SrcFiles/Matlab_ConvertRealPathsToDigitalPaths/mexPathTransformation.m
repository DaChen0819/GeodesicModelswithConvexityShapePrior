
if verLessThan('matlab','8.1')
    cxxFlags = ['CXXFLAGS="-std=c++17" ' ...
        'CXXLIBS="\$CXXLIBS -lc++" ' ]; % This flag is required on some platforms, but must be commented on others...
    outputFlag = '-o ';
else
    cxxFlags = 'CXXFLAGS="-std=c++17" ';
    outputFlag = '-output ';
end





% compile ConvertRealPathsToDigitalPaths.
compilePathTransformation= @(binary_Dir) eval(['mex ' ...
    outputFlag 'ConvertRealPathsToDigitalPaths_JMMLIB' ' ConvertRealPathsToDigitalPaths.cxx' ... % define the filter name.
    ' COMPFLAGS="/std:c++17"' ... % needed on windows platform.
    ' MACOSX_DEPLOYMENT_TARGET=14.0'...% needed on macOS to remove warnings. You have to modify this sentence for your macOS version.
    ' -outdir ' binary_Dir ...
    ' ' cxxFlags...
    ' -I' '../JMM_CPPLibs' ...
    ' -I' 'PathTransformation' ...
    ' -I' 'AlgebraDataSturctures']);

compilePathTransformation(binary_Dir);



