function [mTips,mChosenTip]=DetectEndPointFromRaylineCut(points,adaptiveCut,ParameterizationOrder,maxRadius)

% find the end point for the closed geodesic curve. Let us 't' be the end
% point. Then the direction 't-mOriginPt' is at the right side of the
% direction 'mCutPt-mOriginPt'. We use the cross product between
% 't-mOriginPt' and 'mCutPt-mOriginPt' to detect the point 't'.
% IMPORTANT: (1) point 't' is a 8-neighbour point of 'mCutPt';
%            (2) is NOT involved in the adaptiveCut;
%            (3) 't-mCutPt' aligns with the source point tangent as much as possible.

if nargin==3
    maxRadius=1;
end

% points involve the original point and the source point located at the
% boundary. 
mOriginPt=points.mOriginPt;
mCutPt=points.mCutPt;

if strcmp(ParameterizationOrder,'None')
    ParameterizationOrder='Clockwise';
end

%% find all the admissible neighbour points of the 'mCutPt'.
tangentCut=mCutPt-mOriginPt;
mTips=[];
for i=-maxRadius:maxRadius  % maxRadius is generally set to 1.
    for j=-maxRadius:maxRadius
        if abs(i)+abs(j)~=0
            offset=[i;j];
            mNeigh=mCutPt+offset; % neighbnour point of 'mCutPt'
            tangent1=mNeigh-mOriginPt;
            crossProduct=tangentCut(1).*tangent1(2)-tangent1(1).*tangentCut(2); % cross product. 
            if sum(double(ismember(mNeigh',adaptiveCut','row')))<0.5  % check if the candidates are NOT located at the adaptiveCut.
                if strcmp(ParameterizationOrder,'Clockwise') && crossProduct<0 
                    mTips=cat(2,mTips,mNeigh);
                elseif strcmp(ParameterizationOrder,'AntiClockwise') && crossProduct>0
                    mTips=cat(2,mTips,mNeigh);
                end
            end
        end
    end
end

%% there are two ways to give the source point tangent. 
% 1. by using the 'mTangentPt'
% 2. by using the 'srcTangent' directly. 
if isfield(points,'mTangentPt')
    mTangentPt=points.mTangentPt;
    srcTangent=mTangentPt-mCutPt;
elseif isfield(points,'srcTangent')
    srcTangent=points.srcTangent;
end

% find the end point 't' such that 't-mCutPt' aligns with 'srcTangent' as
% much as possible.
srcTangent=srcTangent./sqrt(srcTangent(1)^2+srcTangent(2)^2);
minScalarProd=-1e10;
minIndex=1;
for j=1:size(mTips,2)
    mTip=mTips(:,j);
    tangentCurve=mCutPt-mTip;
    tangentCurve=tangentCurve./sqrt(tangentCurve(1)^2+tangentCurve(2)^2);
    scalarProdTem=tangentCurve(1).*srcTangent(1)+tangentCurve(2).*srcTangent(2);
    if scalarProdTem>=minScalarProd
        minIndex=j;
        minScalarProd=scalarProdTem;
    end
end
mChosenTip=mTips(:,minIndex);
    
end
