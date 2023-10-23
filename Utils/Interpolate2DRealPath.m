function newRealPath=Interpolate2DRealPath(rawPath,admScale,step)
% interpolate a 2D real path by adding vertices between two path points with
% large Euclidean distance.
   assert(size(rawPath,1)==2||size(rawPath,1)==3);
   assert(admScale>step);
   newRealPath=rawPath(:,1);
   lastVertice=rawPath(:,1);
   
   for i=2:size(rawPath,2)
       currentVertice=rawPath(:,i);
       tangentVec=currentVertice-lastVertice;
       euclidLength=sqrt(tangentVec(1).^2+tangentVec(2).^2);
       if euclidLength<=admScale
           newRealPath=cat(2,newRealPath,currentVertice);
       else
           segment=Interpolate2DSegment(lastVertice, currentVertice,step,euclidLength);
           newRealPath=cat(2,newRealPath,segment);
       end
       lastVertice=currentVertice;
   end
   
end

function segment=Interpolate2DSegment(lastVertice, currentVertice,step,euclidLength)
   segment=lastVertice;
   Num=floor(euclidLength/step);
   tangentVec=(currentVertice-lastVertice)./euclidLength;
   
   for i=1:Num
       segment=cat(2,segment,lastVertice+tangentVec.*i.*step);
   end
   segment=cat(2,segment,currentVertice);
end

