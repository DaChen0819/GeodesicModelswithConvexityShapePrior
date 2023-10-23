function srcOrien=TransferTangentToOrientation(srcTangent)
assert(norm(srcTangent)~=0);
srcOrien=atan2(srcTangent(2),srcTangent(1));
if srcOrien<=0 
    srcOrien=2*pi+srcOrien;
end
srcOrien=rem(srcOrien,2*pi);
end