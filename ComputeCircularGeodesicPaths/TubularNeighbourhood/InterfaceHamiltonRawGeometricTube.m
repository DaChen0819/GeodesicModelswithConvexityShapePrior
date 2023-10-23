function output=InterfaceHamiltonRawGeometricTube(options,cCenterline,posFlags)
assert(nargin==2 || nargin==3);

if nargin==3
    assert(length(cCenterline)==length(posFlags));
end

input_TubularNeigh.seeds=round(cCenterline);

if nargin==2
    posFlags=[];
end


if isfield(options,'imageSize')
    numRows=options.imageSize(1);
    numCols=options.imageSize(2);
else
    error('Image Size must be provided.');
end

gridScale=1.0;
origin=[0;0];


if isfield(options,'tubeDistanceModel')
    tubeDistanceModel=options.tubeDistanceModel;
else
    tubeDistanceModel='Isotropic2';
    disp('Warning: using default setting such that tubeDistanceModel=Isotropic2.');
end


if isfield(options,'acuteAnglePreservedShape')
    input_TubularNeigh.acuteAnglePreservedShape = options.acuteAnglePreservedShape;
else
    input_TubularNeigh.acuteAnglePreservedShape=1;
    disp('Warning: using default setting such that acuteAnglePreservedShape=1.');
end

if isfield(options,'minRadiusTolerance')
    input_TubularNeigh.minRadiusTolerance = options.minRadiusTolerance;
else
    if input_TubularNeigh.acuteAnglePreservedShape>0.5
        input_TubularNeigh.minRadiusTolerance=2.0;
        disp('warning: using default setting that minRadiusTolerance=2.');
    end
end

if isfield(options,'tubeAppearance')
    tubeAppearance = options.tubeAppearance;
else
    tubeAppearance='symmetric';
    disp('warning: using default setting that tubeAppearance=symmetric.');
end

if isfield(options,'angleScalarProductTolerance')
    input_TubularNeigh.angleScalarProductTolerance = options.angleScalarProductTolerance;
else
    if input_TubularNeigh.acuteAnglePreservedShape>0.5
        input_TubularNeigh.angleScalarProductTolerance=0.05;
        disp('warning: using default setting that angleScalarProductTolerance=0.05.');
    end

end

if isfield(options,'maximalTubeWidth')
    maximalTubeWidth  =  options.maximalTubeWidth;
else
    maximalTubeWidth = 10.0;
    disp('Warning: maximalTubeWidth is set to 10 for deault.');
end

if ~isfield(options,'exportSymmetricTube')
   exportSymmetricTube=0.;
else
   exportSymmetricTube=options.exportSymmetricTube;
end

% default setting.
if ~isempty(posFlags)
    if size(posFlags,1)>size(posFlags,2)
        input_TubularNeigh.seedFlags=posFlags;
    else
        input_TubularNeigh.seedFlags=posFlags';
    end
end

input_TubularNeigh.verbosity=0.0;
input_TubularNeigh.arrayOrdering='YXZ_ColumnMajor';
input_TubularNeigh.factoringMethod='None';
input_TubularNeigh.dims=[numCols;numRows];
input_TubularNeigh.exportValues=1.0;
input_TubularNeigh.exportActiveNeighs=0.0;
input_TubularNeigh.exportGeodesicFlow=0.0;
input_TubularNeigh.order=1.0;
input_TubularNeigh.exportActiveOffsets=0.0;
input_TubularNeigh.gridScale=gridScale;
input_TubularNeigh.origin=origin;

input_TubularNeigh.stopAtDistance=maximalTubeWidth;

input_TubularNeigh.tubularAppearance='symmetric';
input_TubularNeigh.model='Isotropic2';
input_TubularNeigh.cost=1.0;
output.symTubularShape=SymmetricNeighbourhood(input_TubularNeigh);


if strcmp(tubeAppearance,'asymmetric')
    
    input_TubularNeigh=rmfield(input_TubularNeigh,'minRadiusTolerance');
    input_TubularNeigh=rmfield(input_TubularNeigh,'angleScalarProductTolerance');
    input_TubularNeigh=rmfield(input_TubularNeigh,'model');
    input_TubularNeigh=rmfield(input_TubularNeigh,'tubularAppearance');
    input_TubularNeigh=rmfield(input_TubularNeigh,'cost');
    input_TubularNeigh.mask=output.symTubularShape;
    input_TubularNeigh.tubularAppearance='asymmetric';
    input_TubularNeigh.acuteAnglePreservedShape=0.0;
    input_TubularNeigh.model=tubeDistanceModel;
    if strcmp(tubeDistanceModel,'Isotropic2') ||strcmp(tubeDistanceModel,'Isotropic3')
        input_TubularNeigh.cost=options.tubularMetric;
    else
        input_TubularNeigh.metric=options.tubularMetric;
    end

    output.asymTubularShape=AntiSymmetricNeighbourhood(input_TubularNeigh);
    
    if exportSymmetricTube<0.5
        output=rmfield(output,'symTubularShape');
    end
end

end


function tubularShape=SymmetricNeighbourhood(input_TubularNeigh)
    output=HamiltonSymmetricTubularNeighbourhood(input_TubularNeigh);
    tubularShape=output.tubularNeighbourhood;
end


function tubularShape=AntiSymmetricNeighbourhood(input_TubularNeigh)
    output=HamiltonAntiSymmetricTubularNeighbourhood(input_TubularNeigh);
    tubularShape=output.tubularNeighbourhood;
end
