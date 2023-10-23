function [binaryImage,mDigitalPath]=ConvertRealPathsToDigitalPaths(mRealPath,imgSize,loopRemoval)
assert(nargin==2 ||nargin==3);

if nargin==2
    loopRemoval=0;
end

gridScale=1.;
origin=[0;0];


if isempty(mRealPath)
    mDigitalPath=[];
    binaryImage=[];
else
    mRealPath=Interpolate2DRealPath_InLine(mRealPath,1.0,0.5);
    input.loopRemoval=loopRemoval;
    input.imageSize=imgSize;
    input.realPath=RescaledCoords(mRealPath,origin,[gridScale;gridScale],false);
    output=ConvertRealPathsToDigitalPaths_JMMLIB(input);
    mDigitalPath=RescaledCoords(output.digitalPath,origin,[gridScale;gridScale]);
    binaryImage=output.binaryImage;
end

end
function mNewRealPath=Interpolate2DRealPath_InLine(mRawPath,admScale,step)
% interpolate a 2D real path by adding vertices between two path points with
% large Euclidean distance.
   assert(nargin==1 ||nargin==3);
   assert(size(mRawPath,1)==2||size(mRawPath,1)==3);
   assert(admScale>step);
   mNewRealPath=mRawPath(:,1);
   lastVertice=mRawPath(:,1);
   if nargin==1
       admScale=1;
       step=0.5;
   end
   
   for i=2:size(mRawPath,2)
       currentVertice=mRawPath(:,i);
       tangentVec=currentVertice-lastVertice;
       euclidLength=sqrt(tangentVec(1).^2+tangentVec(2).^2);
       if euclidLength<=admScale
           mNewRealPath=cat(2,mNewRealPath,currentVertice);
       else
           segment=Interpolate2DSegment_InLine(lastVertice, currentVertice,step,euclidLength);
           mNewRealPath=cat(2,mNewRealPath,segment);
       end
       lastVertice=currentVertice;
   end

end


function segment=Interpolate2DSegment_InLine(lastVertice, currentVertice,step,euclidLength)
   segment=lastVertice;
   Num=floor(euclidLength/step);
   tangentVec=(currentVertice-lastVertice)./euclidLength;
   
   for i=1:Num
       segment=cat(2,segment,lastVertice+tangentVec.*i.*step);
   end
   segment=cat(2,segment,currentVertice);
end