function visualizedImg=VisualizeScribbles(scribbleImageFG,scribbleImageBG,rawImg)

[numRows,numCols,~]=size(rawImg);
if size(rawImg,3)==3
    RC=rawImg(:,:,1);
    GC=rawImg(:,:,2);
    BC=rawImg(:,:,3);
elseif size(rawImg,3)==1
    RC=rawImg(:,:,1);
    GC=rawImg(:,:,1);
    BC=rawImg(:,:,1);
end

for i=1:numRows
    for j=1:numCols
        if scribbleImageFG(i,j)>0.5
            RC(i,j)=0;
            GC(i,j)=0;
            BC(i,j)=255;
        end
        if scribbleImageBG(i,j)>0.5
            RC(i,j)=255;
            GC(i,j)=0;
            BC(i,j)=0;
        end
    end
end
visualizedImg=cat(3,RC,GC,BC);

end