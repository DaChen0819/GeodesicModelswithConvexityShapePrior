function [listSeeds,taggedBoundaries,mNewExternalBoundaries]=SampleBoundaryPoints(GTImage,numSeeds,ratio, timesInitialization,timesSample)
    
    [mBoundaries,~,~,~] = bwboundaries(GTImage>0.5,'noholes');
    mExternalBoundary=mBoundaries{1};
    
    lengTotal=size(mExternalBoundary,1);
    lengSegment=floor(ratio*lengTotal/numSeeds);
    lengSafty=floor((1.0-ratio)*lengTotal/numSeeds);
    taggedBoundaries=cell(1,timesInitialization);
    mNewExternalBoundaries=cell(1,timesInitialization);
    listSeeds=cell(timesInitialization,timesSample);
    for k=1:timesInitialization
        firstIndex=floor(rand(1)*lengTotal/2);
        mNewExternalBoundary=circshift(mExternalBoundary',[0 firstIndex]);
        taggedBoundary=TagContourPoints(mNewExternalBoundary,numSeeds,lengSegment,lengSafty);
        taggedBoundaries{k}=taggedBoundary;
        mNewExternalBoundaries{k}=mNewExternalBoundary;
        for i=1:timesSample
            mSeeds=SamplePointsOnce(taggedBoundary,mNewExternalBoundary,numSeeds,lengSegment);
            listSeeds{k,i}=mSeeds;
        end
    end
end

function taggedBoundary=TagContourPoints(mExternalBoundary,numSeeds,lengSegment,lengSafty)

    taggedBoundary=-1*ones(1,length(mExternalBoundary));
    for l=1:2*numSeeds-1
        if rem(l,2)==0
            passingLength=lengSegment*l/2+lengSafty*(l/2-1);
            taggedBoundary((passingLength+1):(passingLength+lengSafty))=-1;
        else
            passingLength=lengSegment*(l-1)/2+lengSafty*(l-1)/2;
            taggedBoundary((passingLength+1):(passingLength+lengSegment))=1+(l-1)/2;
        end
    end
end

function mSeeds=SamplePointsOnce(taggedBoundary,mExternalBoundary,numSeeds,lengSegment)
   mSeeds=zeros(2,numSeeds);
   idxSeeds=max(floor(lengSegment*rand(1,numSeeds)),1);
   for i=1:numSeeds
       cands=find(taggedBoundary==i);
       mSeeds(1,i)=mExternalBoundary(2,cands(idxSeeds(i)));
       mSeeds(2,i)=mExternalBoundary(1,cands(idxSeeds(i)));
   end
end