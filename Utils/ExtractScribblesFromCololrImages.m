function scribbleLabel=ExtractScribblesFromCololrImages(scribbleImage,FGColor,BGColor)
assert(nargin==1 || nargin==3);

if nargin==1
    FGColor=[0;0;255];
    BGColor=[255;0;0];
end

[numRows,numCols,nz]=size(scribbleImage);
if nz~=3
    error("ERROR:color image must be provided.");
end
scribbleLabel=zeros(numRows,numCols);
for i=1:numRows
    for j=1:numCols
        if scribbleImage(i,j,1)==FGColor(1) && scribbleImage(i,j,2)==FGColor(2) && scribbleImage(i,j,3)==FGColor(3)
            scribbleLabel(i,j)=1;
        end
        
        if scribbleImage(i,j,1)==BGColor(1) && scribbleImage(i,j,2)==BGColor(2) && scribbleImage(i,j,3)==BGColor(3)
            scribbleLabel(i,j)=-1;
        end
    end
end

end