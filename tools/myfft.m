function Y = myfft(X,dim)
%MYFFT peforms normalized DFT on dimension dim of X

if ~exist('dim','var')
    dim = 1;
end

n = size(X,dim);
Y = fft(X,n,dim) ./ sqrt(n);
end

