function shapeGradient=Gradient_PiecewiseConstant(rawImg,shape)
  assert(nargin==2);
  [nx,ny,nz]=size(rawImg);
  shapeGradient=zeros(nx,ny);
  rawImg=double(rawImg);
  shape=double(shape>0.5);
  shape_Out=1.0-shape; % region outside the shape.
        
  if nz==1 %gray level image
      rawImg=rawImg./max(rawImg(:));
      fitting_In=sum(sum(rawImg.*shape))/sum(shape(:));
      fitting_Out=sum(sum(rawImg.*(1-shape)))/sum(shape_Out(:));
      
      shapeGradient=((rawImg-fitting_In).^2-(rawImg-fitting_Out).^2);
 
  elseif nz==3 %color image     
      redChannel=rawImg(:,:,1);     
      greenChannel=rawImg(:,:,2);     
      blueChannel=rawImg(:,:,3);    
      
      redChannel=redChannel./max(redChannel(:));
      greenChannel=greenChannel./max(greenChannel(:));
      blueChannel=blueChannel./max(blueChannel(:));
      
      redFitting_In=sum(sum(redChannel.*shape))/sum(shape(:));     
      redFitting_Out=sum(sum(redChannel.*(1-shape)))/sum(shape_Out(:));
      
      greenFitting_In=sum(sum(greenChannel.*shape))/sum(shape(:));  
      greenFitting_Out=sum(sum(greenChannel.*(1-shape)))/sum(shape_Out(:));
      
      blueFitting_In=sum(sum(blueChannel.*shape))/sum(shape(:));    
      blueFitting_Out=sum(sum(blueChannel.*(1-shape)))/sum(shape_Out(:));
      
      shapeGradient=(redChannel-redFitting_In).^2-(redChannel-redFitting_Out).^2+(greenChannel-greenFitting_In).^2-(greenChannel-greenFitting_Out).^2+(blueChannel-blueFitting_In).^2-(blueChannel-blueFitting_Out).^2;    
      shapeGradient=shapeGradient/3.0;  
  end
end

