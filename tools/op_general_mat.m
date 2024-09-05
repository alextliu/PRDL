function B = op_general_mat(X,A,forward,left)
%OP_GENERAL_MAT performs linear transformation A on rows or columns of X
%(and its adjoint operations)
%
% Input:
%   - forward: logical; indicates whether performs forward map or adjoint
%   - left: logical; indicates whether performs A on rows or columns

if isdiag(X) && all( abs(diag(X)-1) < eps ) % X is an identity matrix
        if forward % forward map
            B = A;
        else % adjoint
            B = A';
        end
else % X is a general matrix
    if isdiag(A) % matrix A is diagonal
        if left % row transformation
            if forward % forward map
                B = diag(A) .* X;
            else % adjoint
                B = conj(diag(A)) .* X;
            end
        else % column transformation
            if forward % forward map
                B = X .* diag(A).';
            else % adjoint
                B = X .* diag(A)';
            end
        end
    else % matrix A is general
        if left % row transformation
            if forward % forward map
                B = A*X;
            else % adjoint
                B = A'*X;
            end
        else % column transformation
            if forward % forward map
                B = X*A;
            else % adjoint
                B = X*A';
            end
        end
    end
end

