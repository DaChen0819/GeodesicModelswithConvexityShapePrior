clc; clear variables;close all;


imageIndex=40; 


%define how many discrete orientations will be used in the S^1:=[0,2\pi] dimention.
numOriens=60;

%read the image img_40.png
imagePath=strcat('Images/img_',num2str(imageIndex),'.png');
rawImg=imread(imagePath);
[numRows,numCols,nz]=size(rawImg);
options.imageSize=[numRows;numCols;numOriens];

%read the scribble data provided by the dataset.
scribblePath=strcat('Images/scribbleImg_',num2str(imageIndex),'.png');

% "scribbleImage_Tem" contains the scribble data. 
% An example for the scribbles can be seen in Fig.~ 5, where the sciibles
% inside the target region are tagged as blue, and otherwise as red (i.e.
% outside the target region.
scribbleImage_Tem=imread(scribblePath);
% "varScribbleImage" is a binary image such that points of the scribble inside the target region is tagged as 1, otherwise -1; 
varScribbleImage=ExtractScribblesFromCololrImages(scribbleImage_Tem,[0;0;255],[255;0;0]);
%Then thin the scribbles. 
erodedData = imerode(255*(varScribbleImage>0),strel('disk',2));
varScribbleImage(varScribbleImage>0)=0;
varScribbleImage(erodedData>0)=1;
clear erodedData; clear scribbleImage_Tem;


%% curve simplicity setting.
%The folloiwing two sentences indicates which method is used for ensouring the curve simplicity. 
% Note that both methods are equvalent from the theoretical
%viewpoint, but differently in numerical implementation. Anyway, both of them have
%almostly identical experimental results. 
options.Simplicity_TotalCurvature=1.0;
options.Simplicity_AngularObstacle=0.0;


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
options.geodesicSolver='Discrete';% another choice is 'ODE'.

% define the direction of the simple closed physical projection curves of
% the computed geodesic paths.
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


%% construct the initial information for staring the computation of geodesic paths.
% "visualizedScribbleImg" is used to visualize the scribbles overlapped in the original image.
visualizedScribbleImg=VisualizeScribbles(varScribbleImage>0,varScribbleImage<0,rawImg);
% The details for how to use scribbles and source point are shown  Fig. 4. 
% Specifically, the blue and red regions indicate the scribbles, and the cyan dot is the "mCutPt". 

%% click points from the image. 
numPts_Clicked=2; %click two points.
% when you click the points, the first point (reffred to as p) should be exactly located at
% the target bounday. The second point (reffered to as q) should be close to the first clicked
% point, AND the direction (q-p) should be positively proportional to the
% direction of the target boundary (assumping that the direction of the
% target boundary is along the AntiClockwise order (defined by
% options.parameterizationOrder='AntiClockwise').
mClickPts=FindClickedPoints(visualizedScribbleImg,numPts_Clicked);

% In column 1 of Fig. 5, the cyan dot represents the "mCutPt", and the
% cyan arrow indicates the direction (q-p), where p is the "mCutPt" and q
% is the "mTangentPt", see bellow.
% The first point tis the "mCutPt", which is the physical position of the source point.
mCutPt=mClickPts(:,1);
% The second point is used to define the angular position for the source
% point.
mTangentPt=mClickPts(:,2);
options.srcOrien=TransferTangentToOrientation(mTangentPt-mCutPt);
options.mCutPoint=mCutPt;

%% Generate the cut and the obstacles from the scriibles and the "mCutPt".
% Details can be found in columns 2 and 3 of Fig. 5. 
% In this figure, "mOriginPt" is indicated by red dot.
% obstacles are indicated by red and blue regions.
[mScribbles,mOriginPt,polyin,mObstacles]=GenerateRayLinesFromScribbles(varScribbleImage,options);
options.mScribbles=mScribbles;
options.mOriginPoint=mOriginPt;
mRayLinePts=cat(2,mOriginPt,mCutPt);
[~,mAdaptiveCut]=SeekRayLineCut(mRayLinePts,[numRows;numCols]);
options.mAdaptiveCut=mAdaptiveCut;

%% -- generate the physical position of the endpoint of gedeosic paths.----
mPoints.mCutPt=mCutPt;
mPoints.mOriginPt=mOriginPt;
mPoints.mTangentPt=mTangentPt;
parameterizationOrder=options.parameterizationOrder;
[~,mChosenTip]=DetectEndPointFromRaylineCut(mPoints,mAdaptiveCut,parameterizationOrder,1);
%  "mStopWhenAnyAccepted" means that once the any of the given endpoints is
%  reached, the fast marching method will be ternimated immediately to
%  reduce the CPU.
%   Note that "mChosenTip" is the physical position of the endpoint.
options.mStopWhenAnyAccepted=mChosenTip;

%% -----parameters setting.
%The parameter "xi" is the weight of the curvature penalization term.
% In our TPAMI paper (Geodesic models with convexity shape prior), this parameter gets to be $\beta$, see Eq. (4).
options.xi=1;
%The parameter "alpha" is the weight of the cost term, see Eq. (69).
options.alpha=4;
% the parameter "mu" is the weight for the region-based homogeneity term,
% see Eq. (70).
options.mu=0.1;
% define the maximal number of the curve iterations.
maxIter_Convex=15;% you can use a large number or give automatic stopping criterion.
% choose the model you want to use for computing gedoesic paths.
options.geoModel='ElasticaConvex2';%DubinsConvex2 %ReedsSheppForwardConvex2 

[mGeoPaths_ConvexityShapedModel,~]=RunCircularConvexityShapedGeodesicModel(options,rawImg,maxIter_Convex);
mConvexGeo_ConvexityShapedModel=mGeoPaths_ConvexityShapedModel{end};


figure(1);imagesc(visualizedScribbleImg); axis image;axis off;hold on;
plot(mConvexGeo_ConvexityShapedModel(1,:),mConvexGeo_ConvexityShapedModel(2,:),'r','Linewidth',3);
plot(mOriginPt(1,:),mOriginPt(2,:),'o','MarkerSize',12,'MarkerFaceColor','r','MarkerEdgeColor','r');
h=quiver(mCutPt(1,:),mCutPt(2,:),20*cos(options.srcOrien),20*sin(options.srcOrien),'b');h.LineWidth=4;h.MaxHeadSize=2;
plot(mCutPt(1,:),mCutPt(2,:),'o','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','b');




return;
r=25;
visualizedScribbleImg=VisualizeScribbles(scribbleMask>1.5,and(scribbleMask<1.5,scribbleMask>0.5),rawImg);
figure(6);imagesc(visualizedScribbleImg);axis off;axis image;hold on;
h=quiver(mCutPt(1,:),mCutPt(2,:),r*cos(options.srcOrien),r*sin(options.srcOrien),'c');h.LineWidth=4;h.MaxHeadSize=2;
plot(mCutPt(1,:),mCutPt(2,:),'o','MarkerSize',12,'MarkerFaceColor','c','MarkerEdgeColor','c');
plot(mOriginPt(1,:),mOriginPt(2,:),'+','MarkerSize',12,'MarkerFaceColor','r','MarkerEdgeColor','r');






