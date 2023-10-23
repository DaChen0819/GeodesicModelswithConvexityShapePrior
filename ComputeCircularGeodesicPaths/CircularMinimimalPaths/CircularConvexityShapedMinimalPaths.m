function dataOut=CircularConvexityShapedMinimalPaths(options,cost,mask)
if isfield(options,'imageSize') && length(options.imageSize)==3
    numRows=options.imageSize(1);
    numCols=options.imageSize(2);
    numThetas=options.imageSize(3);
else
    [numRows,numCols,numThetas]=size(cost);
end

if nargin==2
    mask=ones(numRows,numCols);
end

if isfield(options,'eps')
    input_ConvexGeo.eps= options.eps;
else
    input_ConvexGeo.eps=0.1;
end

if isfield(options,'relaxation')
    input_ConvexGeo.relaxation= options.relaxation;
else
    input_ConvexGeo.relaxation=0.0;
end

if ~isfield(options,'mAdaptiveCut')
     error('Error: mCutPoint should be provided.');
end

if ~isfield(options,'geodesicSolver')
    geodesicSolver='Discrete';
else 
    geodesicSolver=options.geodesicSolver;
end

if isfield(options,'gridScale')
    gridScale=options.gridScale;
else
    gridScale=1.0;
end

if isfield(options,'origin')
    origin=options.origin;
else
    origin=[0.0;0.0];
end

if ~isfield(options,'parameterizationOrder')
    paramOrder='AntiClockwise';
else 
    paramOrder=options.parameterizationOrder;
end

if ~isfield(options,'maxActiveRegionWidth')
    maxActiveRegionWidth=10.0;
else
    maxActiveRegionWidth=options.maxActiveRegionWidth;
end

if isfield(options,'mCutPoint')
    mCutPt=options.mCutPoint;
else
    error('Error: mCutPoint should be provided.');
end

if isfield(options,'srcOrien')
    srcOrien=options.srcOrien;
else
    error('an intial orientation must be provided.')
end

if isfield(options,'mScribbles')
    mScribbles = options.mScribbles;
    cScribbles = RescaledCoords(mScribbles,origin,[gridScale;gridScale],false);
    input_ConvexGeo.scribbles = cScribbles;
    if isfield(options,'scribbleRegion')
        input_ConvexGeo.scribbleRegion=options.scribbleRegion;
    else
        geoLengths=ComputeSpatialGeodesicDistance(options,mScribbles);
        input_ConvexGeo.scribbleRegion=double(geoLengths<=maxActiveRegionWidth);
        clear geoLengths;
    end
    clear mScribbles;
end

if isfield(options,'Simplicity_TotalCurvature')
    input_ConvexGeo.Simplicity_TotalCurvature=options.Simplicity_TotalCurvature;
else 
    input_ConvexGeo.Simplicity_TotalCurvature=1.0;
end

if isfield(options,'Simplicity_AngularObstacle')
    input_ConvexGeo.Simplicity_AngularObstacle=options.Simplicity_AngularObstacle;
else 
    input_ConvexGeo.Simplicity_AngularObstacle=0.0;
end

mAdaptiveCut=options.mAdaptiveCut;
cAdaptiveCut=RescaledCoords(mAdaptiveCut,origin,[gridScale;gridScale],false);
input_ConvexGeo.adaptiveCut=cAdaptiveCut;
if isfield(options,'adaptiveCutRegionTagged')
    input_ConvexGeo.adaptiveCutRegionTagged=options.adaptiveCutRegionTagged;
    options=rmfield(options,'adaptiveCutRegionTagged');
else
    input_ConvexGeo.adaptiveCutRegionTagged=ConstructTaggedAdaptiveCutRegion(options);
end
clear mAdaptiveCut;

cCutPt=RescaledCoords(mCutPt,origin,[gridScale;gridScale],false);
input_ConvexGeo.seeds= [cCutPt+[gridScale*0.1;gridScale*0.1];srcOrien];

%set endpoints.
if isfield(options,'mStopWhenAnyAccepted')
    mPhysicalTips=options.mStopWhenAnyAccepted;
elseif isfield(options,'mStopWhenAllAccepted')
    mPhysicalTips=options.mStopWhenAllAccepted;
elseif isfield(options,'mTips')
    mPhysicalTips=options.mTips;
else
    error('ERROR: something wrong to the provided tips');
end
numTips=size(mPhysicalTips,2);
cPhysicalTips=RescaledCoords(mPhysicalTips,origin,[gridScale;gridScale],false);
for j=1:numTips
    cPhysicalTip=cPhysicalTips(:,j);
    cPhysicalTips(:,j)=cPhysicalTip+[gridScale*0.1;gridScale*0.1];
end
input_ConvexGeo.stopWhenAnyAccepted=[cPhysicalTips;repmat(srcOrien,[1 numTips])];
input_ConvexGeo.parameterizationOrder=paramOrder;
input_ConvexGeo.xi = options.xi; %Model parameter, typical radius of curvature.

input_ConvexGeo.physicalActiveRegion=mask;
input_ConvexGeo.model=options.geoModel;
input_ConvexGeo.cost=cost;
input_ConvexGeo.dims=[numCols;numRows;numThetas];
input_ConvexGeo.exportValues=1; % distance table, of size [n,n,numberOfDirections]

input_ConvexGeo.geodesicStep=0.5;
input_ConvexGeo.verbosity=0;
input_ConvexGeo.exportActiveNeighs=0;
input_ConvexGeo.exportGeodesicFlow=0;
if strcmp(geodesicSolver,'Discrete')
    input_ConvexGeo.geodesicWeightThreshold=0.001;
    input_ConvexGeo.geodesicVolumeBound=10.985;
elseif strcmp(geodesicSolver,'ODE')
    input_ConvexGeo.geodesicCausalityTolerance=1;
    input_ConvexGeo.geodesicTargetTolerance=1;
end
input_ConvexGeo.geodesicSolver=geodesicSolver;
input_ConvexGeo.order=1;
input_ConvexGeo.seedRadius=0;

input_ConvexGeo.origin=origin;  % Physical origin
input_ConvexGeo.gridScale=gridScale; % Physical gridScale
input_ConvexGeo.arrayOrdering='YXZ_ColumnMajor';

% if size(cOriginPt,1)==2
%     input_ConvexGeo.originPoint=cat(1,cOriginPt,0);
% end


if isfield(options,'checkReachTip')
    input_ConvexGeo.checkReachTip=options.checkReachTip;
end

output=HamiltonConvexityShapedCurvaturePenalizedMinimalPaths(input_ConvexGeo);

if isfield(output,"geodesicPoints")
    cGeodesicPaths = mat2cell(output.geodesicPoints,3,output.geodesicLengths); 
    cGeodesicPath=cGeodesicPaths{1};
    cGeodesicPath=cat(2,cGeodesicPath,input_ConvexGeo.seeds);
    dataOut.cGeodesics={cGeodesicPath};
end

if isfield(output,"values")
    dataOut.values=output.values;
end


end



