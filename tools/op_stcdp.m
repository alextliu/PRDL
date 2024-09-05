function Y = op_stcdp(X,mask,forward)
%OP_STCDP performs short-time coded diffraction pattern operation on rows X

[n1,n2] = size(X);
[lenMask,nMask] = size(mask);

if forward
    X_reshape = reshape(X.',lenMask,[]);
else
    X_reshape = reshape(X.',lenMask*nMask,[]);
end
Y = op_cdp(X_reshape,mask,forward);
if forward
    Y = reshape(Y,[n2*nMask,n1]).';
else
    Y = reshape(Y,[n2/nMask,n1]).';
end

end

