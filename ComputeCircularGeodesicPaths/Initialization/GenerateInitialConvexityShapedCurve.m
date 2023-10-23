function initializations=GenerateInitialConvexityShapedCurve(options,rawImg,orienCost)

[numRows,numCols,~]=size(rawImg);

if ~isfield(options,'mCutPoint')
    error('Error: mCutPoint should be provided.');
end

if ~isfield(options,'mOriginPoint')
    error('Error: mOriginPoint should be provided');
end

if ~isfield(options,'srcOrien')
    error('an intial orientation must be provided.')
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

dataOutput=CircularConvexityShapedMinimalPaths(options,orienCost);

digitalPathImage=zeros(numRows,numCols);
if isfield(dataOutput,'cGeodesics')
    cPhysicalGeodesic=dataOutput.cGeodesics{1}(1:2,:);
    mPhysicalGeodesic=RescaledCoords(cPhysicalGeodesic,origin,[gridScale;gridScale]);
    turningAnkle=dataOutput.cGeodesics{1}(3,:);
else
    initializations.mPhysicalGeodesics=[];
    initializations.turningAnkles=[];
    initializations.digitalPathImage=[];
    initializations.mDigPath=[];
    initializations.values=[];
    return;
end
initializations.values=dataOutput.values;


[~,mDigitalPath]=ConvertRealPathsToDigitalPaths(mPhysicalGeodesic,[numRows;numCols]);


indDigPath = sub2ind([numRows numCols],mDigitalPath(2,:)',mDigitalPath(1,:)');
digitalPathImage(indDigPath)=1.0;
initializations.mPhysicalGeodesics={cat(2,mPhysicalGeodesic(:,end),mPhysicalGeodesic)};
initializations.turningAnkles={cat(2,turningAnkle(:,end),turningAnkle)};
initializations.digitalPathImage=digitalPathImage;
initializations.mDigPath=mDigitalPath;


end

