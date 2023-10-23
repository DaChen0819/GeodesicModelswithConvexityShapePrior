function [adaptiveCutRegionTagged,mAdaptiveCutNew]=ConstructTaggedAdaptiveCutRegion(options)

if isfield(options,'parameterizationOrder')
    paramOrder=options.parameterizationOrder;
else
    paramOrder='Clockwise';
end

if isfield(options,'imageSize')
    numRows=options.imageSize(1);
    numCols=options.imageSize(2);
else
    error('Image Size must be provided.');
end

if ~isfield(options,'maxActiveRegionWidth')
    options.maxActiveRegionWidth=10.0;
end

if ~isfield(options,'mAdaptiveCut')
    error('mAdaptiveCut must be provided.');
end

if ~isfield(options,'mOriginPoint')
    error('mOriginPoint must be provided.');
end

mOriginPoint=options.mOriginPoint;
mAdaptiveCut=options.mAdaptiveCut;

if mOriginPoint(1)~=mAdaptiveCut(1,1) || mOriginPoint(2)~=mAdaptiveCut(2,1)
    error('The first point of mAdaptiveCut must be equal to mOriginPoint.');
end

geoLengths_AdaptiveCut=ComputeSpatialGeodesicDistance(options,mAdaptiveCut);  
mask_AdaptiveCut=double(geoLengths_AdaptiveCut<=options.maxActiveRegionWidth); 
clear geoLengths_AdaptiveCut;

%             figure(2);imagesc(mask_AdaptiveCut);axis off;axis image;hold on;
%             plot(options.mAdaptiveCut(:,1),options.mAdaptiveCut(:,2),'r');
%             pause;

options.maxActiveRegionWidth=2*options.maxActiveRegionWidth;
% mOriginPoint=options.mOriginPoint;
options.mask=mask_AdaptiveCut;
geoLengths_OriginPt=ComputeSpatialGeodesicDistance(options,mOriginPoint);
mask_OriginPt=double(geoLengths_OriginPt<=options.maxActiveRegionWidth);

if isfield(options,'mCutPoint')
    mCutPt=options.mCutPoint;
else
    mCutPt=mAdaptiveCut(:,round(size(mAdaptiveCut,2)/2));
end

mMiddlePoint=findMiddlePoint(mAdaptiveCut,mCutPt,3);
geoLengths_CutPt=ComputeSpatialGeodesicDistance(options,mMiddlePoint);


% 2 means anti-clockwise part; 3 means clockwise part; 4 means overlap part.
mask_Overlapped = double(geoLengths_OriginPt<=geoLengths_CutPt).*mask_OriginPt;
mask_Seperated  = double((mask_AdaptiveCut-mask_Overlapped)>0.5);
clear mask_OriginPt;clear geoLengths_CutPt;clear geoLengths_OriginPt;

adaptiveCutRegionTagged=zeros(numRows,numCols);
ind = sub2ind([numRows,numCols],mAdaptiveCut(2,:)',mAdaptiveCut(1,:)');
mask_Seperated(ind)=0;
CC = bwconncomp(mask_Seperated>0.5,8); 

%             figure(1);imagesc(mask_Overlapped);axis off;axis image;hold on;
%             plot(mAdaptiveCut(:,1),mAdaptiveCut(:,2),'r');
%             pause;

if strcmp(paramOrder,'AntiClockwise')
    adaptiveCutRegionTagged(mask_AdaptiveCut>0.5)=2.0; 
    clear mask_AdaptiveCut;
    mask_Tagged=DetermineClockWisePart_TraceBoundary(mAdaptiveCut,mask_Overlapped,CC); 
    
    adaptiveCutRegionTagged(round(mask_Tagged)==3)=3.0;
    adaptiveCutRegionTagged(mask_Overlapped>0.5)=4.0;clear mask_Overlapped;
    adaptiveCutRegionTagged(mCutPt(2),mCutPt(1))=2.0;
elseif strcmp(paramOrder,'Clockwise') 
    adaptiveCutRegionTagged(mask_AdaptiveCut>0.5)=3.0; clear mask_AdaptiveCut;
    mask_Tagged=DetermineClockWisePart_TraceBoundary(mAdaptiveCut,mask_Overlapped,CC);
    adaptiveCutRegionTagged(round(mask_Tagged)==2)=2.0;
    adaptiveCutRegionTagged(mask_Overlapped>0.5)=4.0; clear mask_Overlapped;
    adaptiveCutRegionTagged(mCutPt(2),mCutPt(1))=3.0;
else
    error('Parameterization order is wrong.')
end
if isfield(options,'mVirtualTip')
    [adaptiveCutRegionTagged,mAdaptiveCutNew]=RefineAdaptiveCut(mCutPt,adaptiveCutRegionTagged,options);
end

end
 

function mMiddlePoint=findMiddlePoint(mAdaptiveCut,mCutPt,r)
if nargin==2
    r=3;
end
count=1;
for k=1:size(mAdaptiveCut,2)
    mCurrentPt=mAdaptiveCut(:,k);
    if mCurrentPt(1)==mCutPt(1) && mCurrentPt(2)==mCutPt(2)
        break;
    end
    count=count+1;
end
indexMiddle=min(r,count);
if indexMiddle==1 
    mMiddlePoint=mCutPt;
else
    mMiddlePoint=mAdaptiveCut(:,indexMiddle);
end

end


function maskTagged=DetermineClockWisePart_TraceBoundary(mAdaptiveCut,maskOverlapped,CC)
[nRows,nCols]=size(maskOverlapped);

count=0;
for i=1:size(mAdaptiveCut,2)-1
    pt=mAdaptiveCut(:,i);
    if maskOverlapped(pt(2),pt(1))<0.5
        count=i;
        break;
    end
end
lengthCut=size(mAdaptiveCut,2)-count;
mCentrePt=mAdaptiveCut(:,round(lengthCut/2)+count);
mNextPt=mAdaptiveCut(:,round(lengthCut/2)+count+1);

group1 = CC.PixelIdxList{1};
group2 = CC.PixelIdxList{2};
indCut = sub2ind([nRows,nCols],mAdaptiveCut(2,:)',mAdaptiveCut(1,:)');

x=maskOverlapped(indCut)>0.5;
indCut(x>0.5)=[];


maskTaggedTem=zeros(nRows,nCols);
maskTaggedTem(group1)=1;
maskTaggedTem(indCut)=2;

% trace the boundary 
contourCW = bwtraceboundary(maskTaggedTem>0.5,[mCentrePt(2) mCentrePt(1)],'W',4,Inf,'clockwise');
% imagesc(maskTaggedTem);hold on;
% plot(contourCW(1,2),contourCW(1,1),'xr')
% plot(contourCW(2,2),contourCW(2,1),'gx')
% plot(mCentrePt(1),mCentrePt(2),'gs','Markersize',10)
% plot(mNextPt(1),mNextPt(2),'rs','Markersize',10)


clear maskTaggedTem;
maskTagged=zeros(nRows,nCols);
mNextContourPt=contourCW(2,:);
if mNextPt(1)==mNextContourPt(2) && mNextPt(2)==mNextContourPt(1)
    maskTagged(group1)=3;
    maskTagged(group2)=2;
else
    maskTagged(group1)=2;
    maskTagged(group2)=3;
end


end


function [cutRegionNew,mAdaptiveCutNew]=RefineAdaptiveCut(mCutPoint,cutRegion,options)
numRows=options.imageSize(1);
numCols=options.imageSize(2);
mAdaptiveCut=options.mAdaptiveCut;
refLabel=cutRegion(mCutPoint(2),mCutPoint(1));
flag=CheckNeighPoints(mCutPoint,refLabel,cutRegion,options);

cutRegionNew=cutRegion;
if ~flag
    mAdaptiveCutNew=mAdaptiveCut;
end

if flag
    linearAdaptiveCut = sub2ind([numRows,numCols],mAdaptiveCut(2,:)',mAdaptiveCut(1,:)');
    linearCutPoint    = sub2ind([numRows,numCols],mCutPoint(2,:),mCutPoint(1,:));
    indPos=find(linearAdaptiveCut==linearCutPoint,1);
    currentIndPos=1;
    for j=1:length(linearAdaptiveCut)
        if indPos+j==length(linearAdaptiveCut)-1
            break;
        end
        mPrePoint=mAdaptiveCut(:,indPos+j);
        flag2 = CheckNeighPoints(mPrePoint,refLabel,cutRegion,options);
        if ~flag2
            currentIndPos=j;
            break;
        end
    end
    if refLabel==2
        cutRegionNew(linearAdaptiveCut(indPos+1:indPos+currentIndPos))=3;
    elseif refLabel==3
        cutRegionNew(linearAdaptiveCut(indPos+1:indPos+currentIndPos))=2;
    end
    mTip=mAdaptiveCut(:,indPos+currentIndPos+1);
    mSeed=mCutPoint;
    mPath= ShortestPathLinking(mTip,mSeed,double(round(cutRegionNew)==refLabel),options);
    mAdaptiveCutNew=mAdaptiveCut(:,1:indPos);
    mAdaptiveCutNew=cat(2,mAdaptiveCutNew,mPath);
    if length(linearAdaptiveCut)>=indPos+currentIndPos+2
        mAdaptiveCutNew=cat(2,mAdaptiveCutNew,mAdaptiveCut(:,indPos+currentIndPos+2:end));
    end
end


end

function flag = CheckNeighPoints(mCutPoint,refLabel,adaptiveCutRegionTagged,options)
if isfield(options,'imageSize')
    numRows=options.imageSize(1);
    numCols=options.imageSize(2);
else
    error('Image Size must be provided.');
end

offsets=[1 -1 0 0; 0 0 1 -1];
flag=true;
for k=1:4
    offset=offsets(:,k);
    mNeighPoint=mCutPoint+offset;
    if mNeighPoint(1)<=0 || mNeighPoint(1)>numCols || mNeighPoint(2)<=0 ||mNeighPoint(2)>numRows
        continue;
    end
    currentLabel=adaptiveCutRegionTagged(mNeighPoint(2),mNeighPoint(1));
    if currentLabel~=refLabel
        flag=false;
        break;
    end
end

end

function path= ShortestPathLinking(mTip,mSeed,mask,options)

numRows=options.imageSize(1);
numCols=options.imageSize(2);

input.dims=[numCols;numRows];
input.order=1;
input.seedRadius=0;
input.origin=[0;0];  % Physical origin
input.gridScale=1.0; % Physical gridScale
input.arrayOrdering='YXZ_ColumnMajor';
input.activeRegion=mask;
input.cost=ones(numRows,numCols);

input.seeds=mSeed-1.0;
input.stopWhenAnyAccepted=mTip-1.0;
input.model='Dijkstra2';
input.verbosity=0;
input.exportValues=0.0;
output=DijkstraShortestPathFilter(input);
cGeodesicPaths = mat2cell(output.geodesicPoints,2,output.geodesicLengths);  
mGeoCurve=cGeodesicPaths{1}+1.0;
% imshow(mask);hold on
% plot(mTip(1),mTip(2),'xr')
% plot(mSeed(1),mSeed(2),'xg')
if mGeoCurve(1,end)==mSeed(1) && mGeoCurve(2,end)==mSeed(2)
    mGeoCurve(:,end)=[];
end
path=fliplr(mGeoCurve);
end

