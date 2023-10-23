clc; clear variables;close all;
%% This file implements the methods presented in "geodesic models with convexity shape prior, TPAMI 2023".
% The main goal is to compute a 3D geodesic path (\gamma, \eta) whose physical projection curve \gamma 
% is simple closed and convex, as well as the respective application for interactive image segmentation. 
% In this file, the user intervention is set as landmark points lying in
% the target boundary.  See Section VI-B (1)

%% In these Matlab files,  variables representing points and curves are named using the following criteria:
% if a point  is called mXXXX, it means that  this point lies at a domain
% whose origin is [1,1]. This is the MATLAB convention.
% if a point  is called cXXXX, it means that  this point lies at a domain
% whose origin is [0,0]. This is the C++ convention.

% "imageIndex" is the index of the image to process in the
% ConvexityDataSet.
imageIndex=2;

%define how many points at the target boundary should be provided.
numLandmarkPts=2;
%define how many discrete orientations will be used in the S^1:=[0,2\pi] dimention.
numOriens=60;

%read the image img_2.png
imagePath=strcat('Images/img_',num2str(imageIndex),'.png');
scribblePath=strcat('Images/scribbleImg_',num2str(imageIndex),'.png');
rawImg=imread(imagePath);

%define the image size.
[numRows,numCols,numChannels]=size(rawImg);
%Note that options is a structure data. It includes many parameters and datas.
options.imageSize=[numRows;numCols;numOriens];


%% --------- parameters for tubular neighbourhood --------
%these parameters defined a tubular neighbourhood of the physical projection of the  evolving geodesic paths. 
%Note that if only the image gradients-based cost is applied, there is no
%tubular neighbourhood is required. See Section VI-B (2) Active Geodesic Path Evolution Procedure.
options.maximalTubeWidth=10; % the width of the tubular neighbourhood.
options.tubeDistanceModel='Isotropic2';
options.acuteAnglePreservedShape=0;     
options.tubeAppearance='symmetric';     

%% --------- basic settings for the Hamiltonian Fast Marching Method ----------
%define which kind of numerical scheme to solve the gradient descent ODE for tracking geodesic paths.
options.geodesicSolver='Discrete'; 

%define the direction of the circular geodesic path.
options.parameterizationOrder='AntiClockwise';

%define the physical origin of the discrete domain. The origin of the angular dimension S^1 is 0.
options.origin=[0.0;0.0]; 

%define the grid scale h of the discrete domain. Note that the grid scale is fixed as 2*pi/numOriens.
options.gridScale=2*pi/numOriens;


%% compute the cost function 
%options.shapeHomogeneity indicates that which region-based homogeity
%criterion is used. Specifically: 
%GMM: Gaussian mixture model.
%BhaCoeff: Bhattacharyya Coefficient between two histograms.
%PConst: piecewise constants reduction of the Mumford-Shah model (or the
%region competition model).
options.shapeHomogeneity='GMM';
options.numComponents=5; % GMM components

%% curve simplicity setting.
%The folloiwing two sentences indicates which method is used for ensouring the curve simplicity. 
% Note that both methods are equvalent from the theoretical
%viewpoint, but differently in numerical implementation. Anyway, both of them have
%almostly identical experimental results. 
options.Simplicity_TotalCurvature=1.0;
options.Simplicity_AngularObstacle=0.0;

%% construct the initial information for staring the computation of geodesic paths.
% "mLandmarkPoints" are the physical positions at the target boundary.
% "mCutPt" is the physical position of the source point.
mLandmarkPoints=[391 222;144 324];
mCutPt=[171;273];
% "srcOrien" is the angle of the source point, i.e. [mCutPt,srcOrien] is the source point. 
%  srcOrien defines the direction of the gedoesic path should have at the point "mCutPt".
srcOrien=1.6234; 

%% -- you can also define the "mLandmarkPoints" and "mCutPt" by clicking the images.
% mClickPts=FindClickedPoints(rawImg,2);
% mCutPt=mClickPts(:,1);
% mTangentPt=mClickPts(:,2);
% mLandmarkPoints=FindClickedPoints(rawImg,numLandmarkPts);

%% ---if you do not want to provide the variable "srcOrien", you can use "mTangentPt" instead, as follows:
% options.srcOrien=TransferTangentToOrientation(mTangentPt-mCutPt);

%% store all these points in the structure "options".
options.mLandmarkPoints=mLandmarkPoints;
options.srcOrien=srcOrien;
options.mCutPoint=mCutPt;

%% ---- mOriginPt is the physical position that is inside the target region.
% It is detected as the center of the convex hull of the union of the "mLandmarkpoints" and "mCutPt".
% "mOriginPt" is used to define the "mAdaptiveCut" and the underlying obstacles.
if size(mLandmarkPoints,2)==1
    mOriginPt=round(0.5*(mLandmarkPoints(:,1)+mCutPt));
elseif size(mLandmarkPoints,2)>=2
    [mOriginPt,ployShape]=GenerateOriginalPointFromIntervention(cat(2,mLandmarkPoints,mCutPt),options,'LandmarkPoints');
end
options.mOriginPoint=mOriginPt;
% convexHullShape = poly2mask(ployShape.Vertices(:,1)',ployShape.Vertices(:,2)',numRows,numCols);


%% --- build the cut for computing circular geodesic paths.------
mRayLinePts=cat(2,mOriginPt,mCutPt);
[~,mAdaptiveCut]=SeekRayLineCut(mRayLinePts,[numRows;numCols]);
options.mAdaptiveCut=mAdaptiveCut;

%% ----- build the obstacles using the landmark points at the boundary of the target region.
% For each point p stored in "mLandmarkpoints", we compute a rayline segment emanating from the center of the convex hull (i.e. mOriginalPt) 
% and traveling along the direction "p - mOriginalPt".
% The we remove the point p from this rayline segment such that its remaining part is taken as obstacles such that no paths are 
% allowed to passed through. 
if ~isempty(mLandmarkPoints)
    mScribbles=[];
    for i=1:size(mLandmarkPoints,2)
        mLandmarkPt=mLandmarkPoints(:,i);
        [~,mScribbleCut]=SeekRayLineCut([mOriginPt,mLandmarkPt],[numRows;numCols]);
        mSplittedCut=SplitAdaptiveCut(mScribbleCut,mLandmarkPt,[numRows;numCols],1);
        mScribbles=cat(2,mScribbles,mSplittedCut);
    end
    options.mScribbles=mScribbles; % this is the obstacle.
end


%% -- generate the physical position of the endpoint of gedeosic paths.----
mPoints.mCutPt=mCutPt;
mPoints.mOriginPt=mOriginPt;
mPoints.srcTangent=[cos(srcOrien);sin(srcOrien)];
%% -- if you do not want to use "srcTangent", you can use "mTangentPt" instead.
% mPoints.mTangentPt=mTangentPt;

[~,mChosenTip]=DetectEndPointFromRaylineCut(mPoints,mAdaptiveCut,options.parameterizationOrder,1);
%  "mStopWhenAnyAccepted" means that once the any of the given endpoints is
%  reached, the fast marching method will be ternimated immediately to
%  reduce the CPU.
options.mStopWhenAnyAccepted=mChosenTip; 

%% ------ Note that "mChosenTip" is the physical position of the endpoint.-----

%% -----parameters setting.
%The parameter "xi" is the weight of the curvature penalization term.
% In our TPAMI paper (Geodesic models with convexity shape prior), this parameter gets to be $\beta$, see Eq. (4).
options.xi=1;

%The parameter "alpha" is the weight of the cost term, see Eq. (69).
options.alpha=4;

% the parameter "mu" is the weight for the region-based homogeneity term,
% see Eq. (70).
options.mu=0.5;

% define the maximal number of the curve iterations.
maxIter_Convex=8; % you can use a large number or give automatic stopping criterion.

% choose the model you want to use for computing gedoesic paths.
options.geoModel='ElasticaConvex2'; %DubinsConvex2 %ReedsSheppForwardConvex2 
[mGeoPaths_ConvexityShapedModel,~]=RunCircularConvexityShapedGeodesicModel(options,rawImg,maxIter_Convex);
mGeoPath_ConvexityShapedElastica=mGeoPaths_ConvexityShapedModel{end};


%show the results.
figure(1);imagesc(rawImg); axis image;axis off;hold on;
plot(mGeoPath_ConvexityShapedElastica(1,:),mGeoPath_ConvexityShapedElastica(2,:),'r','Linewidth',3);
h=quiver(mCutPt(1,:),mCutPt(2,:),25*cos(options.srcOrien),25*sin(options.srcOrien),'b');h.LineWidth=4;h.MaxHeadSize=2;
plot(mOriginPt(1,:),mOriginPt(2,:),'o','MarkerSize',12,'MarkerFaceColor','r','MarkerEdgeColor','r');
plot(mCutPt(1,:),mCutPt(2,:),'o','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','b');
plot(mLandmarkPoints(1,:),mLandmarkPoints(2,:),'o','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','b');








