function oriens=FindAdmissableAnkles(mTip,options)

parameterizationOrder=options.parameterizationOrder;
mOriginPoint=options.mOriginPoint;
V1=[cos(options.srcOrien);sin(options.srcOrien)];
V2=mTip-mOriginPoint;
V2=V2./sqrt(V2(1)^2+V2(2)^2);
numThetas=options.imageSize(3);
orienScale=2*pi/numThetas;

oriens=[];
bVec=[];
if strcmp(parameterizationOrder,'AntiClockwise')
    for j=1:numThetas
        Vj=zeros(2,1);
        Vj(1)= V1(1).*cos(orienScale.*j)-V1(2).*sin(orienScale.*j);
        Vj(2)= V1(1).*sin(orienScale.*j)+V1(2).*cos(orienScale.*j);
        b=V2(1).*Vj(2)-V2(2).*Vj(1);
        orien=rem(2*pi+atan2(Vj(2),Vj(1)),2*pi);
        bVec=cat(2,bVec, b);
        if j>1 && (bVec(j-1)*b)<=0
            break;
        end
        oriens=cat(2,oriens,orien);
    end
    return;
elseif strcmp(parameterizationOrder,'Clockwise')
    for j=1:numThetas
        Vj=zeros(2,1);
        Vj(1)= V1(1).*cos(orienScale.*j)+V1(2).*sin(orienScale.*j);
        Vj(2)=-V1(1).*sin(orienScale.*j)+V1(2).*cos(orienScale.*j);
        b=V2(1).*Vj(2)-V2(2).*Vj(1);   
        orien=rem(2*pi+atan2(Vj(2),Vj(1)),2*pi);
        oriens=cat(2,oriens,Vj);
        bVec=cat(2,bVec, b);
        if j>1 && (bVec(j-1)*b)<=0
            break;
        end
        oriens=cat(2,oriens,orien);
    end
    return;
end

