if verLessThan('matlab','8.1')
    cxxFlags = ['CXXFLAGS="-std=c++17" ' ...
        'CXXLIBS="\$CXXLIBS -lc++" ' ]; % This flag is required on some platforms, but must be commented on others...
    outputFlag = '-o ';
else
    cxxFlags = 'CXXFLAGS="-std=c++17" ';
    outputFlag = '-output ';
end


CompileHamiltonSymmetricTube = @(binary_Dir) eval(['mex ' ...
    outputFlag 'HamiltonSymmetricTubularNeighbourhood' ' HamiltonRawGeometricTubularNeighbourhood.cxx' ...
    ' COMPFLAGS="/std:c++17"' ... % needed on windows
    ' MACOSX_DEPLOYMENT_TARGET=14.0'...% needed on macOS to remove warnings. 
    ' -outdir ' binary_Dir ...
    ' ' cxxFlags ...
    ' -I' '../GeoMetrics' ...
    ' -I' '../JMM_CPPLibs' ...
    ' -I' 'FastMarchingBase' ...
    ]);


CompileHamiltonSymmetricTube(binary_Dir);