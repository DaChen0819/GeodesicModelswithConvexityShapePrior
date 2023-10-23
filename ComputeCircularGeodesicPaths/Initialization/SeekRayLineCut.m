function [mRealSegment,mDigitalSegment]=SeekRayLineCut(mPoints,imageSize)

mOriginPt=mPoints(:,1);
mCutPt=mPoints(:,2);

mRealSegment=TrackRayLine(mOriginPt,mCutPt,imageSize);
[~,mDigitalSegment]=ConvertRealPathsToDigitalPaths(mRealSegment,imageSize);

end

function mRealSegment=TrackRayLine(mOriginPt,mCutPt,imageSize)
nRows=imageSize(1);
nCols=imageSize(2);
scale=0.4;
tangentVec=mCutPt-mOriginPt;
leng=sqrt(tangentVec(1).^2+tangentVec(2).^2);
tangentVec=tangentVec./leng;
mRealSegment=mOriginPt;
for i=1:2*max(nRows./scale,nCols./scale)
    mPt=mOriginPt+tangentVec.*scale.*i;
    if mPt(1)<1 || mPt(2)<1 || mPt(1)>nCols || mPt(2)>nRows
        break;
    end
    mRealSegment=cat(2,mRealSegment,mPt);
end
end
