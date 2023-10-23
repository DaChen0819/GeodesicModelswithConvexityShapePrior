function mSplittedCut=SplitAdaptiveCut(mAdaptiveCut,mSrcPts,imgeSize,r)

if nargin==3
    r=1;
end
idxSrc = sub2ind(imgeSize,mSrcPts(2,1),mSrcPts(1,1));
idxCut= sub2ind(imgeSize,mAdaptiveCut(2,:),mAdaptiveCut(1,:));
idxSrc=idxSrc';
idxCut=idxCut';

[pos,~]=ismember(idxCut,idxSrc,'rows');
loc=find(pos==1,1);

if loc>2 && loc < length(idxCut)
    loc=loc-r:loc+r;
    loc(loc>length(idxCut))=[];
    %     loc=[loc-r,loc,loc+r];
elseif loc<=2
    loc=loc+0:r;
%     loc=[loc,loc+r,loc+r+1];
end
idxCut(loc)=[];
[mRows,mCols] = ind2sub(imgeSize,idxCut);

mSplittedCut(1,:)=mCols;
mSplittedCut(2,:)=mRows;


end