function [output] = soft_threshold(a,b)
% SOFT_THRESHOLD Perform soft-thresholding operation.
% b and output are complex vectors/matrices with same dimension;
% a is vectors/matrices with same dimension or scalar;
% elements of a > 0
output = max(abs(b)-a,0).*sign(b);
end

