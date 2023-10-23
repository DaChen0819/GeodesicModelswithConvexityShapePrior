function [mPhysicalGeos,turningAnkles]=RunCircularConvexityShapedGeodesicModel(options,rawImg,maxIters)
%% The interface which which provides the initial curves and the region-based homogeneity features.

%define the maximal number of iterations for evoluting the convexity shaped
%elastica geodeic paths.
if nargin==2
    maxIters=30;
end

%define which type of region-based homogeneity criterion is used. Typical
%examples include the piecewise-constant model (PConst), Gaussian mixture
%model (GMM) and Bhattacharyya coefficient (BhaCoeff) between two histograms.
if isfield(options,'shapeHomogeneity')
    shapeHomogeneity=options.shapeHomogeneity;
else
    shapeHomogeneity='GMM';
end

%define the discretization scale of the image domain. See The first
%paragraph of Section II.
if isfield(options,'gridScale')
    gridScale=options.gridScale;
else
    gridScale=1.0;
end

%define the origin for the physcal domain \Omega. In C++, the origin is 0,
%while for MATLAB, the origin is [1,1]. Note that For both C++ and MATLAB,
%the origin of the angular dimension is set as 0.
if isfield(options,'origin')
    origin=options.origin;
else
    origin=[0.0;0.0];
end

% a parameter for relaxization in the Hamiltonian fast marching method.
if ~isfield(options,'eps')
    options.eps=0.1;
end


% define the size (i.e. columns and rows) of the image to process.
if isfield(options,'imageSize')
    imageSize=options.imageSize;
    numRows=imageSize(1);
    numCols=imageSize(2);
else
    error('imageSize must be provided.');
end

% angularCost_ImageGrad is a cost function defined in Eq. (71).
angularCost_ImageGrad=OrienCostFromImageFeatures(rawImg,options);

% contains the initial curve.
initializations=GenerateInitialConvexityShapedCurve(options,rawImg,angularCost_ImageGrad);
clear orienCost;


if strcmp(shapeHomogeneity,'None') || strcmp(shapeHomogeneity,'none')
    mPhysicalGeos=initializations.mPhysicalGeodesics;
    turningAnkles=initializations.turningAnkles;
    return;
end
if isempty(initializations.mPhysicalGeodesics)
    mPhysicalGeos=[];
    turningAnkles=[];
    return;
end

%% tubular neighbourhood parameters.
tubeOptions=options;
if ~isfield(tubeOptions,'maximalTubeWidth')
    tubeOptions.maximalTubeWidth=10;
end
if ~isfield(tubeOptions,'tubeDistanceModel')
    tubeOptions.tubeDistanceModel='Isotropic2';
end
if ~isfield(tubeOptions,'minRadiusTolerance')
    tubeOptions.minRadiusTolerance=2;
end
if ~isfield(tubeOptions,'tubeAppearance')
    tubeOptions.tubeAppearance='symmetric';
end
if ~isfield(tubeOptions,'acuteAnglePreservedShape')
    tubeOptions.acuteAnglePreservedShape=0;
end
if ~isfield(tubeOptions,'angleScalarProductTolerance')
    tubeOptions.angleScalarProductTolerance=1;%0.05
end

%% extract the information from the initialization.
% mDigPath:  a digital path derived from the intial real curve.
mDigPath=initializations.mDigPath;

% digitalPathImage: a binary image involves the digital path "mDigPath".
digitalPathImage=initializations.digitalPathImage; 
mPhysicalGeos=initializations.mPhysicalGeodesics;
turningAnkles=initializations.turningAnkles;
for iter=1:maxIters
    disp(['iteration=',num2str(iter)]);
    enclosedShape=imfill(digitalPathImage>0.5,'holes'); 
    if strcmp(shapeHomogeneity,'PConst')
        if mod(iter,2)==1
            shapeGrad=Gradient_PiecewiseConstant(rawImg,enclosedShape);
        end
    elseif strcmp(shapeHomogeneity,'GMM') 
        if mod(iter,2)==1
            shapeGrad=Gradient_GaussianMixtureModel(rawImg,enclosedShape,options);
        end
    elseif strcmp(shapeHomogeneity,'BhaCoeff')
        shapeGrad=Gradient_BhattacharyyaCoefficient(rawImg,enclosedShape,options);
    end

    %construct the tubular neighbourhood surrounding the evolving curve,
    %see Section VI-B (2) Active Geodesic Path Evolution Procedure.
    % "mDigPath" is the digital path in the domain of origin [1,1] (this is the matlab convention).
    % "cDigPath" is the corresponding digital path in the domain whose origin is [0,0] (this is the C++ convention).
    cDigPath=RescaledCoords(mDigPath,origin,[1;1],false);
    tubeData=InterfaceHamiltonRawGeometricTube(tubeOptions,cDigPath);
    tubeImage=tubeData.symTubularShape;

    %construct orientation cost.
    orienCost_Region=OrienCostFromImageFeatures(rawImg,options,shapeGrad,tubeImage);
    dataOut_Convex=CircularConvexityShapedMinimalPaths(options,orienCost_Region,tubeImage);
    if isfield(dataOut_Convex,'cGeodesics')
        cPhysicalGeo=dataOut_Convex.cGeodesics{1}(1:2,:);
        mPhysicalGeo=RescaledCoords(cPhysicalGeo,origin,[gridScale;gridScale]);
        turningAnkle=dataOut_Convex.cGeodesics{1}(3,:);
        [digitalPathImage,mDigPath]=ConvertRealPathsToDigitalPaths(mPhysicalGeo,[numRows;numCols]);
    else
        mPhysicalGeo=[];
        turningAnkle=[];
    end
    % save the computed geodesic paths.
    mPhysicalGeos=cat(2,mPhysicalGeos,{cat(2,mPhysicalGeo(:,end),mPhysicalGeo)});
    turningAnkles=cat(2,turningAnkles,{cat(2,turningAnkle(:,end),turningAnkle)});
end

end


