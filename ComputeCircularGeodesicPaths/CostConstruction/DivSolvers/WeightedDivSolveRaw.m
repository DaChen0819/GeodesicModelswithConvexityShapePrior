function [uh,uv] = WeightedDivSolveRaw(shapeGradient,weightDistance)
    assert(ismatrix(shapeGradient));
    assert(all(size(shapeGradient)==size(weightDistance)));
    m=size(shapeGradient,1); n=size(shapeGradient,2);

    weightDistanceh=imfilter(weightDistance,[1;1]/2.,'replicate','full');
    weightDistancev=imfilter(weightDistance,[1,1]/2.,'replicate','full');

    
    %Zero entries in rho are for mask
    mask = weightDistance>eps;
    cstrLoc = find(mask>0); %Indices in array
    cstr = (1:numel(cstrLoc))';
    
    
    uhLoc=find(weightDistanceh>0);
    uvLoc=find(weightDistancev>0);
    
    uhUnk=zeros([m+1,n]);
    uvUnk=zeros([m,n+1]);
    uhUnk(uhLoc) = cstr(end)+(1:numel(uhLoc));
    uvUnk(uvLoc) = uhUnk(uhLoc(end))+(1:numel(uvLoc));
        
    % Creating the constraint matrix
    column=zeros(8*numel(cstrLoc)+numel(uhLoc)+numel(uvLoc),1);
    row=column;
    coeff=column;
    
    pos=1:(4*numel(cstr));
    unit=ones([numel(cstr),1]);
    
    uhUnk1=uhUnk(1:m,       :); uhUnk1=uhUnk1(cstrLoc);
    uhUnk2=uhUnk(2:(m+1),   :); uhUnk2=uhUnk2(cstrLoc);
    uvUnk1=uvUnk(:,       1:n); uvUnk1=uvUnk1(cstrLoc);
    uvUnk2=uvUnk(:,   2:(n+1)); uvUnk2=uvUnk2(cstrLoc);
    
    column(pos) = [uhUnk1;uhUnk2;uvUnk1;uvUnk2];
    row(pos)=[cstr;cstr;cstr;cstr];
    coeff(pos)=[unit;-unit;-unit;unit];
    

    % Creating the Lagrange multiplier matrix, transpose of constraint
    pos2=pos(end)+pos;
    column(pos2)=row(pos);
    row(pos2)=column(pos);
    coeff(pos2)=coeff(pos);
    
    %Weights
    pos3=(pos2(end)+1):numel(column);
    column(pos3) = [uhUnk(uhLoc);uvUnk(uvLoc)];
    row(pos3) = column(pos3);
    coeff(pos3)=[weightDistanceh(uhLoc);weightDistancev(uvLoc)];

    % Rhs and matrix solve
    unkn = numel(cstr)+numel(uhLoc)+numel(uvLoc);
    rhs = [shapeGradient(cstrLoc);zeros(numel(uhLoc)+numel(uvLoc),1)];
    A=sparse(column, row, coeff, unkn, unkn);
    sol = A\rhs;
    
    % Output
    uh=zeros([m+1,n]);
    uv=zeros([m,n+1]);
    uh(uhLoc)=sol(numel(cstrLoc)+(1:numel(uhLoc)));
    uv(uvLoc)=sol(numel(cstrLoc)+numel(uhLoc)+(1:numel(uvLoc)));    
end