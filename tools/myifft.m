function Y = myifft(X,dim)
%MYIFFT peforms normalized inverse DFT on dim of X

if ~exist('dim','var')
    dim = 1;
end

n = size(X,dim);
Y = ifft(X,n,dim) .* sqrt(n);
end

