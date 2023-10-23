function [normalizedCost, cost]=OrienCostFromPerpImageGradients_PureEdge(rawImg,options)

[numRows,numCols,~]=size(rawImg);
imageSize=options.imageSize;

numOriens=imageSize(3);
orienScale=2*pi./numOriens;
gauSigma=options.gaussianSigma;


[Mxx,Mxy,Myy]=ComputePerpImageGradientTensorField(rawImg,gauSigma);
cost=zeros(numRows,numCols,numOriens);
for i=1:numOriens
    a=cos(orienScale*(i-1)+pi/2);
    b=sin(orienScale*(i-1)+pi/2);
    cost(:,:,i)=sqrt(a*a*Mxx+2*a*b*Mxy+b*b*Myy);
end

clear Mxx Mxy Myy;

normalizedCost=(cost-min(cost(:)))./(max(cost(:))-min(cost(:)));

end