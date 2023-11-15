function [mScribbles,mGeoCentroid,polyin,mObstacles]=GenerateRayLinesFromScribbles(scribbleLabelImage,options)

mCutPt=options.mCutPoint;
%scribbles: a matrix where inside: 1, and outside -1.
if isfield(options,'imageSize')
    imageSize=options.imageSize;
    numRows=imageSize(1);
    numCols=imageSize(2);
else
    error('imageSize must be provided.');
end

linearFG=find(round(scribbleLabelImage) == 1);

[mRowFG,mColFG] = ind2sub([numRows;numCols],linearFG);
mRowFG=cat(1,mRowFG,mCutPt(2));
mColFG=cat(1,mColFG,mCutPt(1));
mCoordinatesFG=[mRowFG,mColFG];
indexVertices = convhull(mCoordinatesFG);
polyin = polyshape(mCoordinatesFG(indexVertices,2)',mCoordinatesFG(indexVertices,1)','KeepCollinearPoints',true);
[mGeoCentroid_x, mGeoCentroid_y]= centroid(polyin);
mGeoCentroid_x=round(mGeoCentroid_x);
mGeoCentroid_y=round(mGeoCentroid_y);
mGeoCentroid=[mGeoCentroid_x;mGeoCentroid_y];

potential=zeros(numRows,numCols);
potential(mGeoCentroid(2),mGeoCentroid(1))=1;
distMap = bwdist(potential>0.5);

FGRegion=double(round(scribbleLabelImage) == 1);
BGRegion=double(round(scribbleLabelImage) == -1);

CCFG = bwconncomp(FGRegion>0.5,4);
CCBG = bwconncomp(BGRegion>0.5,4);

mObstacles=[];
mScribbles=[];
for i=1:CCFG.NumObjects
    linearLabels_FG=CCFG.PixelIdxList{i};
    [~,indMin]=min(distMap(linearLabels_FG));
    [mRow_FG_MIN,mCol_FG_MIN] =ind2sub([numRows;numCols],linearLabels_FG(indMin));
    [mRows_FG,mCols_FG] =ind2sub([numRows;numCols],linearLabels_FG);
    mTargetPt=[mCol_FG_MIN;mRow_FG_MIN];
    mDigPathsIn=LinkTwoPointsByRayLineInside(options,mGeoCentroid,mTargetPt);
    mScribbles=cat(2,mScribbles,mDigPathsIn);
    mScribbles=cat(2,mScribbles,[mCols_FG';mRows_FG']);
    mObstacles=cat(2,mObstacles,{mDigPathsIn});
end

for i=1:CCBG.NumObjects
    linearLabels_BG=CCBG.PixelIdxList{i};
    [~,indMin]=min(distMap(linearLabels_BG));
    [mRow_BG_MIN,mCol_BG_MIN] =ind2sub([numRows;numCols],linearLabels_BG(indMin));
    mTargetPt=[mCol_BG_MIN;mRow_BG_MIN];
    mDigPathsOut=LinkTwoPointsByRayLineOutside(options,mGeoCentroid,mTargetPt);
    mScribbles=cat(2,mScribbles,mDigPathsOut);
    [mRows_BG,mCols_BG] =ind2sub([numRows;numCols],linearLabels_BG);
    mScribbles=cat(2,mScribbles,[mCols_BG';mRows_BG']);
    mObstacles=cat(2,mObstacles,{mDigPathsOut});
end


end


function mDigPath=LinkTwoPointsByRayLineInside(options,mOriginPt,mTargetPt)
mDigPath=[];
if mOriginPt(1)==mTargetPt(1) && mOriginPt(2)==mTargetPt(2)
    return;
end

if isfield(options,'imageSize')
    imageSize=options.imageSize;
    numRows=imageSize(1);
    numCols=imageSize(2);
else
    error('imageSize must be provided.');
end

tangentRayLine=mTargetPt-mOriginPt;
tNorm=norm(tangentRayLine);
tangentRayLine(1)=tangentRayLine(1)./tNorm;
tangentRayLine(2)=tangentRayLine(2)./tNorm;
mRealPath=mOriginPt;


for h=0:0.5:10000000
    mCurrentPt=mOriginPt+h.*tangentRayLine;    
    if round(mCurrentPt(1))==mOriginPt(1) && round(mCurrentPt(2))==mOriginPt(2)
        continue;
    end
    mRealPath=cat(2,mRealPath,mCurrentPt);
    if round(mCurrentPt(1))==mTargetPt(1) && round(mCurrentPt(2))==mTargetPt(2)
        break;
    end
end
[~,mDigPath]=ConvertRealPathsToDigitalPaths(mRealPath,[numRows;numCols]);

end

function mDigPath=LinkTwoPointsByRayLineOutside(options,mOriginPt,mTargetPt)

if isfield(options,'imageSize')
    imageSize=options.imageSize;
    numRows=imageSize(1);
    numCols=imageSize(2);
else
    error('imageSize must be provided.');
end

tangentRayLine=mTargetPt-mOriginPt;
tNorm=norm(tangentRayLine);
tangentRayLine(1)=tangentRayLine(1)./tNorm;
tangentRayLine(2)=tangentRayLine(2)./tNorm;
mRealPath=mTargetPt;
for h=0:0.5:10000000
    mCurrentPt=mTargetPt+h.*tangentRayLine;    
    if round(mCurrentPt(1))==mTargetPt(1) && round(mCurrentPt(2))==mTargetPt(2)
        continue;
    end
    mRealPath=cat(2,mRealPath,mCurrentPt);
    if round(mCurrentPt(1))==1 || round(mCurrentPt(2))==1
        break;
    end

    if round(mCurrentPt(1))==numRows || round(mCurrentPt(2))==numRows
        break;
    end

    if round(mCurrentPt(1))==numCols || round(mCurrentPt(2))==numCols
        break;
    end
end
[~,mDigPath]=ConvertRealPathsToDigitalPaths(mRealPath,[numRows;numCols]);

end


