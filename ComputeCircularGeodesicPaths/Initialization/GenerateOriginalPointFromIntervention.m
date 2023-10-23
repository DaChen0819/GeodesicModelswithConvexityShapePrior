function [mGeoCentroid,polyin]=GenerateOriginalPointFromIntervention(interactions,options,model)

if isfield(options,'imageSize')
    imageSize=options.imageSize;
    numRows=imageSize(1);
    numCols=imageSize(2);
else
    error('imageSize must be provided.');
end
if strcmp(model,'LandmarkPoints')
    mCoordinatesFG=zeros(2,size(interactions,2));
    mCoordinatesFG(1,:)=interactions(2,:);
    mCoordinatesFG(2,:)=interactions(1,:);
    mCoordinatesFG=mCoordinatesFG';
    indexVertices = convhull(mCoordinatesFG);
elseif strcmp(model,'Scribbles')
    linearFG=find(round(interactions) == 1);
    [mRowFG,mColFG] = ind2sub([numRows;numCols],linearFG);
    mCoordinatesFG=[mRowFG,mColFG];
    indexVertices = convhull(mCoordinatesFG);
end

polyin = polyshape(mCoordinatesFG(indexVertices,2)',mCoordinatesFG(indexVertices,1)','KeepCollinearPoints',true);
[mGeoCentroid_x, mGeoCentroid_y]= centroid(polyin);
mGeoCentroid_x=round(mGeoCentroid_x);
mGeoCentroid_y=round(mGeoCentroid_y);
mGeoCentroid=[mGeoCentroid_x;mGeoCentroid_y];



end