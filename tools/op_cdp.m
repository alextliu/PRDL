function Y = op_cdp(X,mask,forward)
%OP_CDP performs coded diffraction pattern operation on columns of matrix X

[n1,n2] = size(X);
nM = size(mask,2);

if isdiag(X) && all(abs(diag(X)-1) < eps) % X is an identity matrix
    if forward % forward map
        Y = reshape( permute( repmat(dftmtx(n1)/sqrt(n1), [1,1,nM]) .* permute(conj(mask),[3,1,2]), [1,3,2] ), [n1*nM,n1] );
    else % adjoint
        n3 = n1/nM;
        Y = reshape( repmat(dftmtx(n3)'*sqrt(n3),[1,1,nM]) .* permute(mask,[1,3,2]), [n3,n1] );
    end
else % X is a general matrix
    if forward % forward map
        Y = reshape( myfft( conj(mask) .* permute( repmat(X,[1,1,nM]), [1,3,2] ), 1 ), [n1*nM,n2] );
    else % adjoint
        Y = squeeze( sum( mask .* myifft( reshape(X,[n1/nM,nM,n2]) , 1 ), 2 ) );
    end
end
end

