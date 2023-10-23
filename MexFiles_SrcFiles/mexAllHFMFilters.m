
if verLessThan('matlab','8.1')
    cxxFlags = ['CXXFLAGS="-std=c++17" ' ...
        'CXXLIBS="\$CXXLIBS -lc++" ' ]; % This flag is required on some platforms, but must be commented on others...
    outputFlag = '-o ';
else
    cxxFlags = 'CXXFLAGS="-std=c++17" ';
    outputFlag = '-output ';
end


binary_Dir='../MexFiles_Binary'; %note that you can change this dir to anywhere you want.

% The variable "binary_Dir" will pass to the next four mex functions. 
cd('Matlab_ConvertRealPathsToDigitalPaths');
mexPathTransformation;


cd('../Matlab_ConvexityShapedCurvaturePenalizedMinimalPaths');
mexHamiltonianConvexityShapedMinimalPaths;
cd(binary_Dir);

cd('../Matlab_GeodesicDistanceComputation');
mexHamiltonianDistanceComputation;
cd(binary_Dir);

cd('../Matlab_RawGeometricTubularNeighbourhood');
mexRawGeometricTubularNeighbourhood;


