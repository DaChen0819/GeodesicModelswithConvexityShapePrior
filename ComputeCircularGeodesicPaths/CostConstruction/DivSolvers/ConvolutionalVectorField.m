function u = ConvolutionalVectorField(shapeGradient,tube, greenKernel, normal)
% ConvolutionalVectorField       
%     Inputs
%     xi          shape gradient map, d1-by-d2 matrix for 2D, and d1-by-d2-by-d3 matrix for 3D.   
%     K           Greem kernel
%     spacing     a 2 or 3-element vector. In 3D case, spacing=[Rx Ry Rz] define the spacing between 
%                 voxels for x, y and z dimension, default value is [1 1 1].
%               
%     Outputs
%     u           the vector field. For 2D, d1-by-d2-by-2 matrix, 
%                 For 3D, d1-by-d2-by-d3-by-3 matrix.
% 

%% - convlove xi with K using FFT to speed up the process
assert(nargin==3 || nargin==4);

if (nargin==3)
    normal=false;
end

shapeGradient=shapeGradient.*double(tube>0.5);

if ismatrix(shapeGradient) % 2D
    FFTsize = size(shapeGradient) + size(greenKernel(:,:,1)) - 1;
    k = greenKernel(:,:,1) + 1i*greenKernel(:,:,2);      % consider the vectors as complex numbers
    temp = ifft2(fft2(shapeGradient, FFTsize(1), FFTsize(2)) .* fft2(k, FFTsize(1), FFTsize(2)));

    % remove padded points
    rmv = (size(greenKernel(:,:,1)) - 1)/2;
    temp = temp(rmv(1)+1:end-rmv(1), rmv(2)+1:end-rmv(2), :);
    u(:,:,2) = imag(temp);
    u(:,:,1) = real(temp);
    u=double(u);
    
    temu=u(:,:,1);
    u(:,:,1)=u(:,:,2);
    u(:,:,2)=-temu;
    
    if normal
        disp('Convolutional Kernel-based Method WITH Normalization.');
        hKernel=ones(size(greenKernel,1),size(greenKernel,2))/(size(greenKernel,1)*size(greenKernel,2));
        normalizedFactor=conv2(double(tube>0.5),hKernel,'same');
        u(:,:,1)=u(:,:,1)./(eps+normalizedFactor);
        u(:,:,2)=u(:,:,2)./(eps+normalizedFactor);
    else
        disp('Convolutional Kernel-based Method WITHOUT Normalization.');
    end
    
else % 3D
    FFTsize = size(shapeGradient) + size(greenKernel(:,:,:,1)) - 1;
    fftf = fftn(shapeGradient, FFTsize);    
    u(:,:,:,3) = real(ifftn(fftf .* fftn(greenKernel(:,:,:,3), FFTsize)));
    temp = ifftn(fftf .* fftn(greenKernel(:,:,:,1) + 1i*greenKernel(:,:,:,2), FFTsize));
    u(:,:,:,2) = imag(temp);
    u(:,:,:,1) = real(temp);
    clear temp fftf;
    
    % remove padded points
    rmv = (size(greenKernel(:,:,:,1)) - 1)/2;
    u = u(rmv(1)+1:end-rmv(1), rmv(2)+1:end-rmv(2), rmv(3)+1:end-rmv(3), :);
    
    u=double(u);
    
end
