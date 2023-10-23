function [shapeGradient,probIn, probOut]=Gradient_BhattacharyyaCoefficient(Img,shape,options,toNormalized)

assert(nargin==3||nargin==4);
[~,~,nz]=size(Img);
numBins=options.numBins;
kernelScale=options.kernelScale;
Img=double(Img);

if nargin==3
    toNormalized=false;
end


if nz==1
    imageData=(Img-min(Img(:)))/(max(Img(:))-min(Img(:)));
    rescaledImage=RescaleImageIntensity(imageData,numBins);
    vecParameters=[2;numBins;kernelScale];
elseif nz==3
    imageDataR=(Img(:,:,1)-min(Img(:)))/(max(Img(:))-min(Img(:)));
    imageDataG=(Img(:,:,2)-min(Img(:)))/(max(Img(:))-min(Img(:)));
    imageDataB=(Img(:,:,3)-min(Img(:)))/(max(Img(:))-min(Img(:)));
    rescaledImageR=RescaleImageIntensity(imageDataR,numBins);
    rescaledImageG=RescaleImageIntensity(imageDataG,numBins);
    rescaledImageB=RescaleImageIntensity(imageDataB,numBins);
    
    rescaledImage=floor(rescaledImageG*numBins+rescaledImageB*numBins.^2+rescaledImageR);
    vecParameters=[2;numBins^3;kernelScale];
end

[shapeGradient,probIn,probOut]=BhattacharyyaCoefficient(double(rescaledImage), double(shape),vecParameters);

if toNormalized
% //    disp('Use Normalized Procedure');
% //    shapeGradient(shapeGradient<=0)=shapeGradient(shapeGradient<=0)./max(-shapeGradient(:));
% //    shapeGradient(shapeGradient>0)=shapeGradient(shapeGradient>0)./max(shapeGradient(:));
    disp('DONOT Use Normalized Procedure!!!');
end

end

function rescaledImage=RescaleImageIntensity(imageData,numBins)
imageData=double(imageData);
normalizedImage=imageData./max(imageData(:));
rescaledImage=floor(1.0+normalizedImage.*numBins);
end

