function [Txx,Txy,Tyy]=ComputePerpImageGradientTensorField_MatlabEigen(rawImg,gauSigma)

assert(nargin==1 || nargin==2 || nargin==3);
if nargin==2
    gauSigma=1.0;
end

[numRows,numCols,~]=size(rawImg);

[~,~,numChannels]=size(rawImg);
rawImg=im2double(rawImg);

if numChannels==3
    [rGrad_x,rGrad_y]=imgradientxy(imgaussfilt(rawImg(:,:,1),gauSigma),'central');
    [gGrad_x,gGrad_y]=imgradientxy(imgaussfilt(rawImg(:,:,2),gauSigma),'central');
    [bGrad_x,bGrad_y]=imgradientxy(imgaussfilt(rawImg(:,:,3),gauSigma),'central');
    
    Mxx=1.0+rGrad_x.^2+gGrad_x.^2+bGrad_x.^2;
    Myy=1.0+rGrad_y.^2+gGrad_y.^2+bGrad_y.^2;
    Mxy=rGrad_y.*rGrad_x+gGrad_y.*gGrad_x+bGrad_y.*bGrad_x;
    
else
    [Gradx,Grady] = grayLevelGradient(rawImg,gauSigma);
    Mxx=1.0+Gradx.^2;
    Myy=1.0+Grady.^2;
    Mxy=Grady.*Gradx;
end

Txx=zeros(numRows,numCols);
Txy=zeros(numRows,numCols);
Tyy=zeros(numRows,numCols);
for i=1:numRows
    for j=1:numCols
        a=Mxx(i,j);
        b=Mxy(i,j);
        c=Myy(i,j);
        matA=[a b;b c];
        matInvA=inv(matA);
        Txx(i,j)=matInvA(1,1);
        Txy(i,j)=matInvA(1,2);
        Tyy(i,j)=matInvA(2,2);
    end
end


end

function [Grad_x,Grad_y] = grayLevelGradient(I,gaussianSigma)
[~,~,numChannels]=size(I);
if numChannels>1
    error('Input Image should be gray level image');
end
smoothImg=imgaussfilt(I,gaussianSigma);  %Guassian filter to smooth the image.
[Grad_x,Grad_y]=imgradientxy(smoothImg,'central');
end

