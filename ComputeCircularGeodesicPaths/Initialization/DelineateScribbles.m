function [mScribbles,scribbleImage]=DelineateScribbles(rawImg,radius,stringWords)
%DelineateScribbles is used to draw scribbles in the images.
%When the picture appears, left click means 'start to draw'. At the same
%time, your mouse will become '+'. After drawing one scribble, you need one
%more left click to start to draw a new scribble.
% Click 'Space' key means 'stop drawing'.

if nargin==1
    radius=3;
    stringWords='Draw Scribbles';
elseif nargin==2
    stringWords='Draw Scribbles';
end

[numRows,numCols,~]=size(rawImg);
scribbleImage=zeros(numRows,numCols);
figure;imagesc(rawImg); axis image;axis off;
title(stringWords);
mScribbles=[];
while true
    w=waitforbuttonpress;
    if w==1
        close all;
        break;
    else
        h = drawfreehand('Closed',false,'Multiclick',false);
        ctrPoints=round(h.Position);
        ctrPoints=ctrPoints';
        mLines=InterpolateBetweenPoints(ctrPoints);
        mScribbles=cat(2,mScribbles,mLines);
    end
end
if isempty(mScribbles)
    return;
end
ind = sub2ind([numRows,numCols],mScribbles(2,:),mScribbles(1,:));
BW=zeros(numRows,numCols);
BW(ind)=1;
D = bwdist(BW>0.5);
clear BW;

ind=find(D<radius);
scribbleImage=zeros(numRows,numCols);
scribbleImage(ind)=1;
[mRows,mCols] = ind2sub([numRows,numCols],ind);
mRows=mRows';
mCols=mCols';
mScribbles=[mCols;mRows];

% strcmpi(get(gcf,'CurrentCharacter'),'e')
              

end


function mLines=InterpolateBetweenPoints(ctrPoints)

numPoints=size(ctrPoints,2);
mLines=[];
if numPoints==1
    return;
end
for i=1:numPoints-1
   pCurrent= ctrPoints(:,i);
   pNext= ctrPoints(:,i+1);
   tangent=pNext-pCurrent;
   tangent=tangent./sqrt(tangent(1)^2+tangent(2)^2);
   if isnan(tangent(1))|| isnan(tangent(2))
       continue;
   end
   mLines=cat(2,pCurrent,mLines);
   for j=1:10000000
       pTem=pCurrent+0.2.*tangent.*j;
       if round(pTem(1))~=mLines(1,end) || round(pTem(2))~=mLines(2,end) 
           mLines=cat(2,round(pTem),mLines);
       end

       if round(pTem(1))==pNext(1) && round(pTem(2))==pNext(2) 
           break;
       end
   end
end



end