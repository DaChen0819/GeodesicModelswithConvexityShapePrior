function [mScribblesIn,mScribblesOut,scribbleImage]=GenerateScribblesFromLandmarkPoints(mLandmarkPoints,options,r)

if nargin==2
    r=2;
end

if isfield(options,'imageSize')
    imageSize=options.imageSize;
    numRows=imageSize(1);
    numCols=imageSize(2);
else
    error('imageSize must be provided.');
end

if isfield(options,'mOriginPoint')
    mOriginPt=options.mOriginPoint;
else
    error('Error: mOriginPoint should be provided');
end

mScribblesIn=[];
mScribblesOut=[];
for i=1:size(mLandmarkPoints,2)
    mLandmarkPt=mLandmarkPoints(:,i);
    [~,mRayLine]=SeekRayLineCut([mOriginPt,mLandmarkPt],[numRows;numCols]);
    [mScribbleIn,mScribbleOut]=GenerateScribbles(mRayLine,mLandmarkPt,[numRows;numCols],r);
    mScribblesIn=cat(2,mScribblesIn,mScribbleIn);
    mScribblesOut=cat(2,mScribblesOut,mScribbleOut);
end

scribbleImage=zeros(numRows,numCols);
mLinearScribblesIn  = sub2ind(imageSize,mScribblesIn(2,:),mScribblesIn(1,:));
mLinearScribblesOut = sub2ind(imageSize,mScribblesOut(2,:),mScribblesOut(1,:));
scribbleImage(mLinearScribblesIn)  = 1;
scribbleImage(mLinearScribblesOut) =-1;

end


function [mScribblesIn,mScribblesOut]=GenerateScribbles(mRayline,mSrcPts,imageSize,r)
numRows=imageSize(1);
numCols=imageSize(2);

mLinearSrc = sub2ind(imageSize,mSrcPts(2,1), mSrcPts(1,1));
mLinearCut = sub2ind(imageSize,mRayline(2,:),mRayline(1,:));
mLinearSrc=mLinearSrc';
mLinearCut=mLinearCut';

[pos,~]=ismember(mLinearCut,mLinearSrc,'rows');
loc=find(pos==1);

mLinearCutIn=mLinearCut(1:max(loc-r,1));
mLinearCutOut=mLinearCut(min(length(pos),loc+r):end);
[mRowsIn,mColsIn]   = ind2sub([numRows;numCols],mLinearCutIn);
[mRowsOut,mColsOut] = ind2sub([numRows;numCols],mLinearCutOut);

mScribblesIn(1,:)=mColsIn;
mScribblesIn(2,:)=mRowsIn;

mScribblesOut(1,:)=mColsOut;
mScribblesOut(2,:)=mRowsOut;


end
