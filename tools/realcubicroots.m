function r = realcubicroots(coef)
%REALCUBICROOTS Find real roots of cubic polynomial.
% If the equation has only 1 real root or 3 equal real roots,return a 
% single root;
%
% If the equation has 3 different real roots (2 of them may be equal), 
% return 3 roots.
%
% Output r is a column vector.

if ~isreal(coef)
    coef = real(coef);
end

tmp = coef(1) * coef(3);
Delta0 = coef(2)^2 - 3 * tmp;
Delta1 = 2 * coef(2)^3 - 9 * tmp * coef(2) + 27 * coef(1)^2 * coef(4);
Delta = Delta1^2 - 4*Delta0^3;

if Delta >= 0
    % Only 1 real root or 3 equal real roots
    C = nthroot( ( Delta1 + sqrt(Delta) ) / 2, 3 );
    if abs(C)<eps
        C = nthroot( ( Delta1 - sqrt(Delta) ) / 2, 3 );
    end
    r = -1/3/coef(1) .* ( coef(2) + C + Delta0 / C );
else
    % 3 different real roots, 2 of them may be equal
    C = ( ( Delta1 + sqrt(Delta) ) / 2 )^(1/3);
    if abs(C)<eps
        C = ( ( Delta1 - sqrt(Delta) ) / 2 )^(1/3);
    end
    ksi = [ 1; -0.5 + 1i*(sqrt(3)/2); -0.5 - 1i*(sqrt(3)/2) ];
    ksiC = ksi.*C;
    r = -1/3/coef(1) .* ( coef(2) + real(ksiC + Delta0 ./ ksiC) );
end
end

