function K = ConstructGreenKernel(Rmax, Dim, spacing)
% Compute the Green vector field kernel.
%     K = ConstructGreenKernel(R_max, type);
%     Inputs
%     R_max       the maximal  radius of the kernel.
%     type        output dimension: '2D' for 2D and '3D' for 3D kernels respectively.
%     spacing     spacing vector.

%     Outputs
%     K           the Green kernel
%
% Note that this function is adapted from the AM_VFK function obtained  from  http://viva.ee.virginia.edu/.
%%
r0 = floor(Rmax./spacing);
if Dim == 2
    [x,y] = meshgrid(spacing(1)*(r0(1):-1:-r0(1)),spacing(2)*(r0(2):-1:-r0(2)));
    dist = sqrt(x.*x+y.*y);
    SetZero = (dist>Rmax);
    x(SetZero) = 1e-2;
    y(SetZero) = 1e-2;
    
    k1=x./(eps+dist.^2);
    k2=y./(eps+dist.^2);

    K = double(cat(3,k1,k2));
elseif Dim == 3
    [x,y,z] = meshgrid(spacing(1)*(r0(1):-1:-r0(1)),spacing(2)*(r0(2):-1:-r0(2)),spacing(3)*(r0(3):-1:-r0(3)));
    dist = sqrt(x.*x+y.*y+z.*z);
    SetZero = (dist>Rmax);
    x(SetZero) = 0;
    y(SetZero) = 0;
    z(SetZero) = 0;
        
    x=x./(eps+dist.^2);
    y=y./(eps+dist.^2);
    z=z./(eps+dist.^2);
    
    K = double(cat(4,x,y,z));
else
    error('Dim must be 2 or 3!');   
end
