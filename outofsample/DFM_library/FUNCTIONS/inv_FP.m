function [ xInv ] = inv_FP( x )
    
    xInv = x\eye(size(x));

end

