function potential=PiecewiseConstantTubePotential(shape,ShapeGradient,v,tau)

if nargin<=1 
    error("the number of input should be 3 or 4.");
end

if nargin==2
    v=2;
    tau=0.01;
end

if nargin==3
    tau=0.01;
end

[nx,ny]=size(shape);

cShape=1.0*shape;
normalizedShapeGradient=ShapeGradient./max(max(abs(ShapeGradient)));
maxSG=max(max(abs(normalizedShapeGradient) ));

D1=double(normalizedShapeGradient>tau*maxSG).*double(cShape);
D2=double(normalizedShapeGradient<-tau*maxSG).*double(1-cShape);
D3=abs(normalizedShapeGradient)<tau*maxSG;

potential=ones(nx,ny)/(eps+v);
potential(D1>0.5)=v;
potential(D2>0.5)=v;
potential(D3>0.5)=1.0;

end
