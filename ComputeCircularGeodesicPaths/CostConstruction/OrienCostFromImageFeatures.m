function orienCost=OrienCostFromImageFeatures(rawImg,options,shapeGradient,tube)
assert(nargin==2 ||nargin==4);

if nargin==2
    useRegionBasedFeatures=false;
elseif nargin==4
    useRegionBasedFeatures=true;
end


if ~isfield(options,'imageSize')
    error('imageSize must be provided.');
end

if ~isfield(options,'gaussianSigma')
    options.gaussianSigma=1.0;
end

if isfield(options,'gridScale')
    options.gridScale=1.0;
end

alpha=options.alpha;
if useRegionBasedFeatures
    mu=options.mu;
    [normalizedCostE,~]=OrienCostFromPerpImageGradients(rawImg,options);
    [~,costR]=OrienCostFromRegionalHomogeneity(shapeGradient,tube,options);
    orienCost=exp(alpha.*(normalizedCostE+mu.*costR));
else
    [normalizedCostE,~]=OrienCostFromPerpImageGradients(rawImg,options);
    orienCost=exp(alpha.*normalizedCostE);
end

end

function [normalizedCost, cost]=OrienCostFromPerpImageGradients(rawImg,options)

[numRows,numCols,~]=size(rawImg);
imageSize=options.imageSize;

numOriens=imageSize(3);
orienScale=2*pi./numOriens;
gauSigma=options.gaussianSigma;

if exist("ComputePerpImageGradientTensorField","file")==2
    [Mxx,Mxy,Myy]=ComputePerpImageGradientTensorField(rawImg,gauSigma);
elseif exist("ComputePerpImageGradientTensorField_MatlabEigen","file")==2
    [Mxx,Mxy,Myy]=ComputePerpImageGradientTensorField_MatlabEigen(rawImg,gauSigma);
end

cost=zeros(numRows,numCols,numOriens);
for i=1:numOriens
    a=cos(orienScale*(i-1)+pi/2);
    b=sin(orienScale*(i-1)+pi/2);
    cost(:,:,i)=sqrt(a*a*Mxx+2*a*b*Mxy+b*b*Myy);
end

clear Mxx; clear Mxy; clear Myy; clear a; clear b; 

normalizedCost=(cost-min(cost(:)))./(max(cost(:))-min(cost(:)));

end

function [normalizedCost,cost]=OrienCostFromRegionalHomogeneity(shapeGradient,tube,options)

assert(nargin==3);

gridScale=options.gridScale;

if isfield(options,'PDESolver')
    PDESolver=options.PDESolver;
else
    PDESolver='div';
end

if isfield(options,'imageSize')
    imageSize=options.imageSize;
else
    error('imageSize must be provided.');
end

numRows=imageSize(1);
numCols=imageSize(2);

nTheta=imageSize(3);
orienScale=2*pi/nTheta;

%estimate the vector field
if strcmp(PDESolver,'div')
    u = DivSolve(shapeGradient,double(tube>0.5),gridScale);
elseif strcmp(PDESolver,'conv')
    u = ConvolutionalVectorField(shapeGradient,double(tube>0.5), greenKernel);
end


if isfield(options,'parameterizationOrder')
    paraOrder=options.parameterizationOrder;
else
    paraOrder='Clockwise';
end

u_x=u(:,:,1);
u_y=u(:,:,2);
uNorm=(u_x.^2+u_y.^2).^0.5+eps;
uNorm_Max=max(uNorm(:));clear uNorm;

if strcmp(paraOrder,'AntiClockwise')
    u_y=-1.0.*u_y;
    u_x=-1.0.*u_x;
end

uNormalized_x=u_x./uNorm_Max;
uNormalized_y=u_y./uNorm_Max;

cost=zeros(numRows,numCols,nTheta);
for i=1:nTheta
    a=cos(orienScale*(i-1));
    b=sin(orienScale*(i-1));
    cost(:,:,i)=a.*uNormalized_x+b.*uNormalized_y;
end
clear Txx Tyy Txy uNormalized_x uNormalized_y;
tube3D=repmat(tube,[1 1 nTheta]);
normalizedCost=(cost-min(cost(:)))./(max(cost(:))-min(cost(:)));
normalizedCost(tube3D<0.5)=1e10;
clear tube3D;

end
