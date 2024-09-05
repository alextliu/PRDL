function [D,Z,patternZ,stat] = FUN_PRDLsca(data,params,D0,Z0)

% This function solves the following phase retrieval with dictinary
% learning problem:
%
% min_{D,Z} 0.5*|| Y - abs(A*D*Z*B) ||^2_F + lambda*||Z||_1
% s.t. {columns of D have unit norm}
%
% Using Successive Convex Approximation
%
% INPUT:
%   data:               measurements, measurement parameters, ground-truth 
%                       values
%       .A:             spatial mixer
%       .B:             temporal mixer
%       .Y:             nosiy magnitude-only measurements
%       .L:             ground-truth sparsity level
%   params:             algorithm parameters
%       .P:             no. of columns in dictionary D
%       .lambda:        sparsity regularization parameter
%       .blockRule:     0: jointly update D and Z in each iteration; 1:
%                       alternatively update D or Z in each iteration.
%       .maxIter:       max no. of iterations
%       .verb:          display algorithm information in each iteration if
%                       verb = 1
%
% OPTIONAL INPUT:
%   D0, Z0 - initialization of variables
%
% OUTPUT:
%   D, Z - solutions for the variables
%   patternZ - support of solution Z
%   stat - algorithm statistics
%
% Author: Tianyi Liu

% disp('Algorithm begins...');

%% Initializiation
timer = tic;
format compact;
tol_rational = 1e-9;    % tolerance for rational approximation

A = data.A; % spatial mixer
B = data.B; % temporal mixer
Y = data.Y; % noisy magnitude-only measurements
L = data.L; % ground-truth sparsity level

if isfield(params,'verb')
    verb = params.verb;
else
    verb = 0;
end
if verb
    var_str = ['DZ'];
end

if isfield(params,'trackConvergence')
    trackConvergence = params.trackConvergence;
else
    trackConvergence = false;
end

% Initializes parameters
N = params.N;   % # rows of D
I = params.I;   % # columns of Z
P = params.P;   % # columns of D
lambda = max(params.lambda,0);  % sparsity regularization parameter
debias = (params.lambda<0);     % debiasing mode if params.lambda<0

if isfield(params,'maxIter')
    maxIter = params.maxIter;
else
    maxIter = 1000;
end
if isfield(params,'tol')
    tol = params.tol;
else
    tol = 1e-4;
end
% number of block variables not to update in each iteration
if isfield(params,'blockRule')
    blockRule = params.blockRule;
else
    blockRule = 0;  % default: joint update
end

% If A and B are matrices, make function handles
if isnumeric(A)
    A = @(XX,forward) op_general_mat(XX,A,forward,1);
end
if isnumeric(B)
    B = @(XX,forward) op_general_mat(XX,B,forward,0);
end

% Remove negative measurements if any
Y(Y<0) = 0;

% initialize dictionary
if( ~exist('D0','var') || isempty(D0) )
    D0 = (randn(N,P) + 1i * randn(N,P));
end
D0 = D0./sqrt(sum(abs(D0).^2,1));   % normalize each column
D = D0;

% initialize code
if( ~exist('Z0','var') || isempty(Z0) )
    Z0 = (randn(P,I) + 1i * randn(P,I))/sqrt(2);
end
Z = Z0;

absZ = abs(Z);
patternZ = (absZ > eps);   % indicators of non-zero entries of Z
if lambda > eps
    lambdaL1normZ = lambda * sum(absZ(:));
else
    lambdaL1normZ = 0;
end
if debias
    if sum(patternZ(:)) == 0
        if verb
            fprintf('Invalid data: debiasing phase with empty support for Z!\n');
        end
        return;
    end
    Z(~patternZ) = 0;
end

[U,S,V] = svd(A(eye(N),1),'econ');    % thin SVD
SUH = diag(S) .* U';    % S*U'
normBrows2 = sum(abs(B(eye(I),1)).^2, 2).';    % squared l2-norm of rows of B, row vector

% Allocate memory
solApproxD = zeros(N,P);
solApproxZ = zeros(P,I);
deltaZ = zeros(P,I);
subgradZ = zeros(P,I);
coef = zeros(5,1);

% Initial intermediate variables, gradient and objective function
AD = A(D,1);
dAD = sum(abs(AD).^2,1).';    % squared l2 norm of each column of A*D, column vector
normDcols2 = sum(abs(D).^2);    % squared l2-norm of columns of D
boundaryDInd = (normDcols2 >= 1 - tol_rational);   % indices of columns of D that are on the boundary of feasible set

ZB = B(Z,1);

Yhat = AD*ZB;
signYhat = sign(Yhat);
signYhat(signYhat == 0) = 1;
Yt = Y .* signYhat;  % Y.*exp(1i*angle(ADZ));
residual = Yt - Yhat; %residual

obj = 0.5*norm( residual ,'fro')^2 + lambdaL1normZ;

% compute gradients of smooth part of upper bound of objective
gradientD = - A( residual * ZB', 0 );   % gradient of smooth part of upper bound of objective w.r.t. D
gradientZ = - B( AD'* residual, 0 );   % gradient of smooth part of upper bound of objective w.r.t. Z

updateInd = true(2,1);  % indicate the variable blocks to be updated
stationaryInd = false(2,1); % indicate the variable blocks (D and Z) that have achieved stationarity with the given tolerance


% Initialize struct to record statistics:
recstats = (nargout > 3) || verb;
if( recstats )
    stat.runtime = zeros(maxIter+1,1);
    stat.objVal = zeros(maxIter+1,1);
%     stat.distance = zeros(maxIter,1);
    stat.improve = zeros(maxIter,1);
    stat.stepsize = zeros(maxIter,1);
    stat.error_subgradD = zeros(maxIter,1);
    stat.error_subgradZ = zeros(maxIter,1);
    stat.sparsityZ = zeros(maxIter+1,1);    % record average sparstiy level of columns of Z
    stat.rmse_sparsityZ = zeros(maxIter+1,1);
    
    stat.objVal(1) = obj;
    stat.sparsityZ(1) = sum(patternZ(:)) / I;
    stat.rmse_sparsityZ(1) = sqrt( mean( abs( sum(patternZ) - L ) .^2 ) );

    stat.runtime(1) = toc(timer);
end

if trackConvergence
    % Evaluate optimality of current point
    % minimum-norm sub-gradient
    if lambda > eps
        subgradZ(:) = 0;
        subgradZ(patternZ) = gradientZ(patternZ) + lambda * sign(Z(patternZ));
        subgradZ(~patternZ) = max( abs(gradientZ(~patternZ))-lambda, 0);    % sign(gradientZ(~nnzInd)) .* max( abs(gradientZ(~nnzInd))-lambda, 0); phase can be dropped since we need only its norm
%         error_subgradZ = norm(subgradZ,'fro') * distanceZ / abs(obj);
%         error_subgradZ = max( abs(subgradZ(:)) );
        error_subgradZ = norm(subgradZ,'fro') / sqrt(numel(Z)) /numel(Y);
    else
        if debias
%             error_subgradZ = max( abs(gradientZ(patternZ)) );
%             error_subgradZ = norm(gradientZ(patternZ)) * distanceZ / abs(obj);
            error_subgradZ = norm(gradientZ(patternZ)) / sqrt(sum(patternZ,'all')) / numel(Y);
        else
%             error_subgradZ = max(abs(gradientZ(:)));
%             error_subgradZ = norm(gradientZ,'fro') * distanceZ / abs(obj);
            error_subgradZ = norm(gradientZ,'fro') / sqrt(numel(Z)) / numel(Y);
        end
    end

    subgradD = gradientD;
    subgradD(:,boundaryDInd) = subgradD(:,boundaryDInd) + sqrt( sum( abs(gradientD(:,boundaryDInd)).^2, 1 ) ) .* D(:,boundaryDInd);
%     error_subgradD = max( abs(subgradD(:)) );
%     error_subgradD = norm(subgradD,'fro') * distanceD / abs(obj);
    error_subgradD = norm(subgradD,'fro') / sqrt(numel(D)) / numel(Y);

    error_subgrad = sqrt( (error_subgradD^2*numel(D) + error_subgradZ^2*numel(Z)) / (numel(D) + numel(Z)) );

    % evaluate and record estimation quality every 10 iterations
    clear sol;
    sol.D = D;
    sol.Z = Z;
    quality = eval_quality(sol,data);
    stat.quality.iterArray = [0;0];
    stat.quality.timeArray = [0; stat.runtime(1)];
    stat.quality.FmeasureZ = [quality.FmeasureZ; quality.FmeasureZ];
    stat.quality.errorD = [quality.errorD; quality.errorD];
    stat.quality.errorZ = [quality.errorZ; quality.errorZ];
    stat.quality.errorDZ = [quality.errorDZ; quality.errorDZ];
    stat.quality.errorX = [quality.errorX; quality.errorX];
    stat.quality.error_subgrad = [error_subgrad; error_subgrad];
    stat.quality.objVal = [obj; obj];
end


% distanceD = 1;
% distanceZ = 1;

% Main algorithm
for t = 1: maxIter
    timer = tic;
    
    if verb
        fprintf('iter %d:\n',t);
    end

    %% Select variable blocks to be updated
    if blockRule
        % Alternatively update each block of variables
        updateInd(1) = logical(mod(t,2));
        updateInd(2) = ~updateInd(1);
    else
        % Jointly update all variables
        updateInd(:) = true;
    end
%     updateInd = updateInd & (~stationaryInd);
    
    if any(updateInd)
        
        %% Compute descent direction for D
        if updateInd(1)
            %     fprintf('\tS_D norm<1 columns: ');
            
            normZBrows2 = sum(abs(ZB).^2,2);  % column vector, squared norm of rows of Z*B
            nnzZBrowInd = (normZBrows2 > eps);     % indicators of non-zero rows of Z*B
            nnzZBrow = sum(nnzZBrowInd);    % # of non-zero rows of Z*B
            
            % 1. For non-zero rows of Z*B, solve the subproblems for corresponding
            % columns of D using Newton method or rational approximation
            if nnzZBrow > 0
                Num = SUH * ( residual*ZB' + AD .* normZBrows2.' ); % SUH * (Res*Z' + bsxfun(@times,AD,normZrows2'));   % each column corresponds to a column of D; zero column if corresponding row of Z is zero
                Den = diag(S).^2 * normZBrows2.';                           % each column corresponds to a column of D; zero column if corresponding row of Z is zero
                Tmp = zeros(size(S,1),P);
                Tmp(:,nnzZBrowInd) = Num(:,nnzZBrowInd) ./ Den(:,nnzZBrowInd);
                normSDcols2 = sum(abs(Tmp).^2); % row vector, squared norm of cols of SD
                nu = zeros(1,P);    % each element corresponds to the Lagragian multiplier of a column of D
                deriv = zeros(1,P);
                for k = 0:maxIter-1
                    activeColInd = nnzZBrowInd' & ( normSDcols2 > 1 + tol_rational); % indicators of columns to be updated
                    if sum(activeColInd) == 0
                        break;
                    end
                    deriv(activeColInd) = -2 * sum( abs(Num(:,activeColInd)).^2 ./ (Den(:,activeColInd)+nu(activeColInd)).^3 );
                    
                    % 1. Newton method
                    %         nu(activeColInd) = nu(activeColInd) + (1-normSDcols2(activeColInd)) ./ deriv(activeColInd);
                    % 2. Rational approx. from left (norm > 1)
                    nu(activeColInd) = nu(activeColInd) + 2 * normSDcols2(activeColInd) ./ deriv(activeColInd) .* (1-sqrt(normSDcols2(activeColInd)));
                    
                    Tmp(:,activeColInd) = Num(:,activeColInd) ./ (Den(:,activeColInd)+nu(activeColInd));
                    normSDcols2(activeColInd) = sum(abs(Tmp(:,activeColInd)).^2);
                end
                solApproxD(:,nnzZBrowInd) = V*Tmp(:,nnzZBrowInd);
            end
            
            % 2. For zero rows of Z, any solution is optimal for corresponding columns of D
            % S_D(:,~nnzRowInd) = ones(N,P-nnzRow);
            % S_D(:,~nnzRowInd) = D(:,~nnzRowInd);    % keep those columns of D unchanged
            solApproxD(:,~nnzZBrowInd) = randn(N,P-nnzZBrow) +1i*randn(N,P-nnzZBrow);
            solApproxD(:,~nnzZBrowInd) = solApproxD(:,~nnzZBrowInd) ./ sqrt(sum(abs(solApproxD(:,~nnzZBrowInd)).^2));    % bsxfun(@rdivide, S_D(:,~nnzRowInd), sqrt(sum(abs(S_D(:,~nnzRowInd)).^2)));

            deltaD = solApproxD - D;
        end
        
        %% Compute descent direction for Z
        if updateInd(2)
            if debias % debiasing mode
                solApproxZ = Z - gradientZ ./ dAD ./ normBrows2;
%                 S_Z(~patternZ) = 0;   % this step is not necessary because S_Z will not be used in the stepsize computation
                deltaZ(:) = 0;
                deltaZ(patternZ) = solApproxZ(patternZ) - Z (patternZ);
            else
                solApproxZ = soft_threshold(lambda, dAD .* Z .* normBrows2 - gradientZ) ./ dAD ./ normBrows2;  % bsxfun(@rdivide, FUN_soft_threshold(lambda, bsxfun(@times,dAD,Z) - gradientZ), dAD);
                deltaZ = solApproxZ - Z;
            end
        end
        
        % ----------------------------------------------------------------
        %% Compute stepsize
        % ----------------------------------------------------------------
        coef(:) = 0;
        if all(updateInd)
            % Jointly update both D and Z in this iteration
            
%             distance = norm([DeltaD(:);DeltaZ(:)]);
            
            AdeltaD = A(deltaD,1);
            deltaZB = B(deltaZ,1);
            %     M0 = -Res; % A*D*Z - Y;
            M1 = AdeltaD * ZB + AD * deltaZB;
            M2 = AdeltaD * deltaZB;
            if lambda > eps
                m4 = lambda * sum(abs(solApproxZ(:))) - lambdaL1normZ;    % lambda * ( sum(abs(S_Z(:))) - sum(abs(Z(:))) );
            else
                m4 = 0;
            end
            % ---------------------------------------------------------------------
            % coeficients of fourth order polynomal
            % 1/4*coef(1) * gamma^4 + 1/3*coef(2) * gamma^3 + 1/2*coef(3) * gamma^2
            % +coef(4) * gamma + coef(5)
            % ---------------------------------------------------------------------
            coef(1) = 2 * norm(M2,'fro')^2;
            coef(2) = 3 * real( M1(:)' * M2(:) );   % 3 * real(dot(M1(:),M2(:)));
            coef(3) = -2 * real( residual(:)' * M2(:) ) + norm(M1,'fro')^2;  % 2 * real(dot(-Res(:),M2(:))) + norm(M1,'fro')^2;
            coef(4) = - real( residual(:)' * M1(:) ) + m4;   % real(dot(-Res(:),M1(:))) + m4;
            %             coef(5) = 0.5*norm(Res,'fro')^2;   % constant term has no effect on min
            
            if coef(4) >= 0
                if verb
                    fprintf('Not descent direction!\n');
                end
                t = t-1;
                break;
            end
            
            if coef(1) > eps
                rts = realcubicroots(coef(1:4)); % all real roots including multiplicity
                rts = rts(rts<=1 & rts>=0);    % real roots in (0,1)
                if isempty(rts) == 1
                    % no root in [0,1]
                    stepsize = 1;
                elseif length(rts) == 1
                    % only one root in [0,1]
                    stepsize = rts;
                else
                    % Compare the function value if there are multiple roots in [0,1]
                    % RTSreal01 = [0;RTSreal;1];
                    rts1 = [rts;1];
                    polyValRTS = polyval(coef./[4:-1:1,1].',rts1);
                    [~,minIdx] = min(polyValRTS);
                    stepsize = rts1(minIdx);
                end
            elseif coef(3) > eps % line search function is quadratic
                stepsize = min( max( -coef(4) / coef(3), 0 ), 1 );
            else % line search function is linearly decreasing
                stepsize = 1;
            end
            
            % Update variables according to line search
            D = D + stepsize * deltaD;
%             distanceD = stepsize * norm(deltaD,'fro');
            if debias
                Z(patternZ) = Z(patternZ) + stepsize * deltaZ(patternZ);
%                 distanceZ = stepsize * norm(deltaZ(patternZ));
            else
                Z = Z + stepsize * deltaZ;
%                 distanceZ = stepsize * norm(deltaZ,'fro');
            end
            
            % Update intermediate variables
            Yhat = Yhat + stepsize * M1 + stepsize^2 * M2;
            
            AD = AD + stepsize * AdeltaD; % A*D;
            dAD = sum(abs(AD).^2,1).';
            normDcols2 = sum(abs(D).^2);    % squared l2-norm of columns of D
            boundaryDInd = (normDcols2 >= 1 - tol_rational);% indices of columns of D that are on the boundary of feasible set

            ZB = ZB + stepsize * deltaZB;   % Z*B
            if ~debias
                absZ = abs(Z);
                patternZ = (absZ > eps);
                if lambda > eps
                    lambdaL1normZ = lambda * sum(absZ(:));
                end
            end
            
        elseif updateInd(1)
            % Update only D in this iteration
            % The line search function is a quadratic polynomial
            % 1/2*coef(3) * gamma^2 - coef(4) * gamma + const.
            
%             distance = norm(DeltaD,'fro');
            
            % stepsize computation
            AdeltaD = A(deltaD,1);
            %     M0 = -Res; % AD*Z - Y;
            M1 = AdeltaD * ZB;
            coef(3) = norm(M1,'fro')^2;
            coef(4) = real( residual(:)' * M1(:) );  % must be > 0 as a descent direction
            if coef(4) <= 0
                if verb
                    fprintf('Not descent direction!\n');
                end
                t = t-1;
                break;
            end
            if coef(3) > eps % line search function is quadratic
                stepsize = min( max( coef(4) / coef(3), 0 ), 1 );
            else % line search function is linearly decreasing
                stepsize = 1;
            end
            
            % Update D according to line search
            D = D + stepsize * deltaD;
%             distanceD = stepsize * norm(deltaD,'fro');
            
            % Update intermediate variables
            Yhat = Yhat + stepsize * M1; % AD*Z;
            
            AD = AD + stepsize * AdeltaD;   % A*D;
            dAD = sum(abs(AD).^2,1).';
            normDcols2 = sum(abs(D).^2);    % squared l2-norm of columns of D
            boundaryDInd = (normDcols2 >= 1 - tol_rational);% indices of columns of D that are on the boundary of feasible set

        else
            % Update only Z in this iteration
            % The line search function is a quadratic polynomial
            % 1/2*coef(3) * gamma^2 - coef(4) * gamma + const.
            
%             distance = norm(DeltaZ,'fro');
            
            % stepsize computation
            deltaZB = B(deltaZ,1);
            %     M0 = -Res; % AD*Z - Y;
            M1 = AD * deltaZB;
            coef(3) = norm(M1,'fro')^2;
            coef(4) = real( residual(:)' * M1(:) );    
            if lambda > eps
                coef(4) = coef(4) - (lambda * sum(abs(solApproxZ(:))) - lambdaL1normZ);
            end
            if coef(4) <= 0 % must be > 0 as a descent direction
                if verb
                    fprintf('Not descent direction!');
                end
                t = t-1;
                break;
            end
            if coef(3) > eps % line search function is quadratic
                stepsize = min( max( coef(4) / coef(3), 0 ), 1 );
            else % line search function is linearly decreasing
                stepsize = 1;
            end
            
            % Update Z and intermediate variables
            if debias
                Z(patternZ) = Z(patternZ) + stepsize * deltaZ(patternZ);
%                 distanceZ = stepsize * norm(deltaZ(patternZ));
            else
                Z = Z + stepsize * deltaZ;
%                 distanceZ = stepsize * norm(deltaZ,'fro');
            end
            
            %     AD = A*D;
            Yhat = Yhat + stepsize * M1; % AD*Z;
            
            ZB = ZB + stepsize * deltaZB;   % Z*B
            if ~debias
                absZ = abs(Z);
                patternZ = (absZ > eps);
                if lambda > eps
                    lambdaL1normZ = lambda * sum(absZ(:));
                end
            end
            
        end
        
        %% Update intermediate variables and objective function value
        signYhat = sign(Yhat);
        signYhat(signYhat == 0) = 1;
        Yt = Y .* signYhat;  % Y.*exp(1i*angle(ADZ));
        residual = Yt - Yhat; %residual
        
        obj = 0.5*norm( residual ,'fro')^2 + lambdaL1normZ;
        
        % compute gradients of smooth part of upper bound of objective
        gradientD = - A( residual * ZB', 0 );   % gradient of smooth part of upper bound of objective w.r.t. D
        gradientZ = - B( AD'* residual, 0 );   % gradient of smooth part of upper bound of objective w.r.t. Z
        
    else
        % No update in this iteration
%         distance = 0;
        stepsize = 0;
    end

    runtime = toc(timer);
    %% Evaluate optimality of current point
    % minimum-norm sub-gradient
    if lambda > eps
        subgradZ(:) = 0;
        subgradZ(patternZ) = gradientZ(patternZ) + lambda * sign(Z(patternZ));
        subgradZ(~patternZ) = max( abs(gradientZ(~patternZ))-lambda, 0);    % sign(gradientZ(~nnzInd)) .* max( abs(gradientZ(~nnzInd))-lambda, 0); phase can be dropped since we need only its norm
%         error_subgradZ = norm(subgradZ,'fro') * distanceZ / abs(obj);
%         error_subgradZ = max( abs(subgradZ(:)) );
        error_subgradZ = norm(subgradZ,'fro') / sqrt(numel(Z)) /numel(Y);
    else
        if debias
%             error_subgradZ = max( abs(gradientZ(patternZ)) );
%             error_subgradZ = norm(gradientZ(patternZ)) * distanceZ / abs(obj);
            error_subgradZ = norm(gradientZ(patternZ)) / sqrt(sum(patternZ,'all')) / numel(Y);
        else
%             error_subgradZ = max(abs(gradientZ(:)));
%             error_subgradZ = norm(gradientZ,'fro') * distanceZ / abs(obj);
            error_subgradZ = norm(gradientZ,'fro') / sqrt(numel(Z)) / numel(Y);
        end
    end

    subgradD = gradientD;
    subgradD(:,boundaryDInd) = subgradD(:,boundaryDInd) + sqrt( sum( abs(gradientD(:,boundaryDInd)).^2, 1 ) ) .* D(:,boundaryDInd);
%     error_subgradD = max( abs(subgradD(:)) );
%     error_subgradD = norm(subgradD,'fro') * distanceD / abs(obj);
    error_subgradD = norm(subgradD,'fro') / sqrt(numel(D)) / numel(Y);

    if recstats
        stat.runtime(t+1) = stat.runtime(t) + runtime;
        stat.objVal(t+1) = obj;
        %         stat.distance(t) = distance;
        stat.stepsize(t) = stepsize;
        stat.improve(t) = (stat.objVal(t)-stat.objVal(t+1)) / stat.objVal(t);
        if debias
            stat.sparsityZ(t+1) = stat.sparsityZ(t);
            stat.rmse_sparsityZ(t+1) = stat.rmse_sparsityZ(t);
        else
            stat.sparsityZ(t+1) = sum(patternZ(:)) / I;
            stat.rmse_sparsityZ(t+1) = sqrt( mean( abs( sum(absZ > eps) - L ) .^2 ) );
        end
        stat.error_subgradZ(t) = error_subgradZ;
        stat.error_subgradD(t) = error_subgradD;
    end

    if trackConvergence
        % evaluate and record estimation quality every 10 iterations
        if mod(t,10) == 0
            clear sol;
            sol.D = D;
            sol.Z = Z;
            quality = eval_quality(sol,data);
            stat.quality.iterArray = [stat.quality.iterArray; t];
            stat.quality.timeArray = [stat.quality.timeArray; stat.runtime(t+1)];
            stat.quality.FmeasureZ = [stat.quality.FmeasureZ; quality.FmeasureZ];
            stat.quality.errorD = [stat.quality.errorD; quality.errorD];
            stat.quality.errorZ = [stat.quality.errorZ; quality.errorZ];
            stat.quality.errorDZ = [stat.quality.errorDZ; quality.errorDZ];
            stat.quality.errorX = [stat.quality.errorX; quality.errorX];
            error_subgrad = sqrt( (error_subgradD^2*numel(D) + error_subgradZ^2*numel(Z)) / (numel(D) + numel(Z)) );
            stat.quality.error_subgrad = [stat.quality.error_subgrad; error_subgrad];
            stat.quality.objVal = [stat.quality.objVal; obj];
        end
    end
    
    if verb
        fprintf('\tUpdate %s\n',var_str(updateInd));
        fprintf('\tvalue = %f, step-size = %f, improve = %e,\n\tsubgradD = %e, subgradZ = %e,\n\tavg. code sparsity = %f, error sparsity = %f, time = %f s.\n',...
            obj,stepsize,stat.improve(t),error_subgradD,error_subgradZ, stat.sparsityZ(t+1), stat.rmse_sparsityZ(t+1),stat.runtime(t));
    end
    
    %% Stopping criterion
    stationaryInd(1) = (error_subgradD <= tol);
    stationaryInd(2) = (error_subgradZ <= tol);
    if all(stationaryInd)
        if verb
            fprintf('Tolerance achieved. Terminating.\n');
        end
        break;
    end
    
    if stepsize == 0
        if verb
            fprintf('Stepsize is zero. Aborting.\n');
        end
        break;
    end
end

if recstats
%     stat.runtime = toc(timer);
    stat.Niter = t;     % # iterations
    stat.optVal = stat.objVal(t+1);  % 0.5*norm( abs( A*D*Z) - Y ,'fro')^2 + lambda*sum(abs(Z(:)));  % optimal value
    if t < maxIter
        if t > 0
            stat.runtime(t+2:end) = stat.runtime(t+1);
        end
        stat.objVal(t+2:end) = stat.optVal;
        stat.error_subgradD(t+1:end) = stat.error_subgradD(t);
        stat.error_subgradZ(t+1:end) = stat.error_subgradZ(t);
        stat.sparsityZ(t+2:end) = stat.sparsityZ(t+1);
        stat.rmse_sparsityZ(t+2:end) = stat.rmse_sparsityZ(t+1);
    end
end

end

