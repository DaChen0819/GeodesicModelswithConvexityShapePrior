function [shapeGradient,foreground_Prob,background_Prob]=Gradient_GaussianMixtureModel(rawImg,shape,options)

if isfield(options,'numComponents')
    numComponents=options.numComponents;
else
    numComponents=5;
end


[numRows,numCols,nz]=size(rawImg);

if nz==3
    data=im2double(reshape(rawImg,[numRows*numCols,3]));
    foregroundData = data(shape>0.5, :);
    backgroundData = data(shape<=0.5,:);
elseif nz==1
    data=im2double(reshape(rawImg,[numRows*numCols,1]));
    foregroundData = data(shape>0.5);
    backgroundData = data(shape<=0.5);
end
rng(5);
options_GMM = statset('MaxIter',1200);
foreground_GMM = fitgmdist(foregroundData,numComponents,'Options',options_GMM,...
    'CovarianceType','diagonal','RegularizationValue',0.01);

background_GMM = fitgmdist(backgroundData,numComponents,'Options',options_GMM,...
    'CovarianceType','diagonal','RegularizationValue',0.01);
foreground_Prob = foreground_GMM.pdf(data); 
background_Prob = background_GMM.pdf(data);

foreground_Prob = reshape(foreground_Prob,[numRows,numCols]);
background_Prob = reshape(background_Prob,[numRows,numCols]);

shapeGradient   = log(background_Prob./foreground_Prob);
end









