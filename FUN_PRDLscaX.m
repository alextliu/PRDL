function [X,D,Z,patternZ,stat] = FUN_PRDLscaX(data,params,X0,D0,Z0)

% This function solves the following phase retrieval problem with
% dictionary learing problem:
%
% min  0.5*||Y-abs(A*X*B)||^2_F + 0.5*mu*||X-D*Z||^2_F
% s.t. D contains unit columns
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
%       .mu:            weight of dictionary learing term
%       .lambda:        sparsity regularization parameter
%       .blockRule:     0: updateInd all blocks X,D,Z in each iteration;
%                       1: 1 block not updated in each iteration;
%                       2: 2 blocks not updated in each iteration;
%       .maxIter:       max no. of iterations
%       .verb:          display algorithm information in each iteration if
%                       verb = 1
%
% OPTIONAL INPUT:
%   X0, D0, Z0 - initialization of variables
%
% OUTPUT:
%   X, D, Z - solutions for the variables
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

% default values for optional input parameters
if isfield(params,'verb')
    verb = params.verb;
else
    verb = 0;
end
if verb
    var_str = ['XDZ'];
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
mu = params.mu; % sparse approximation regularization parameter
lambda = max(params.lambda,0);  % sparsity regularization parameter
debias = (params.lambda<0);     % debiasing mode if params.lambda<0

if mu == 0 % regularization parameter mu must > 0
    if verb
        disp('Invalid data: mu = 0!');
    end
    return;
end

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

% initialize signals
if( ~exist('X0','var') || isempty(X0) )
    X0 = (randn(N,I) + 1i*randn(N,I))/sqrt(2);
end
X = X0;

% initialize dictionary
if( ~exist('D0','var') || isempty(D0) )
    %     D0 = [eye(s),dict_dct(s1,s)];
    D0 = (randn(N,P) + 1i * randn(N,P));
end
D0 = D0./sqrt(sum(abs(D0).^2,1));   % normalize each column
D = D0;

% initialize code
if( ~exist('Z0','var') || isempty(Z0) )
    Z0 = D0\X0;
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

% Allocate memory
solApproxX = zeros(N,I);
deltaX = zeros(N,I);
solApproxD = zeros(N,P);
deltaD = zeros(N,P);
solApproxZ = zeros(P,I);
deltaZ = zeros(P,I);
subgradZ = zeros(P,I);
coef = zeros(5,1);

% Initial intermediate variables, gradient and objective function
lambda_mu = lambda / mu;
% HessianX = A'*A + mu*eye(N);    % Hessian of smooth part of objective w.r.t. each column of X and same for all columns
% normAcols2 = sum( abs(A).^2 ).';  % squared l2 norm of columns of A, column vector
% normBrows2 = sum( abs(B).^2, 2 ).';   % squared l2 norm of rows of B, row vector
HessianXpart = sum( abs( A(eye(N),1) ).^2 ).' * sum( abs( B(eye(I),1) ).^2, 2 ).' + mu;    % second-order partial derivatives of smooth part of objective w.r.t. each entry of X

normDcols2 = sum(abs(D).^2);    % squared l2-norm of columns of D
boundaryDInd = (normDcols2 >= 1 - tol_rational);% indices of columns of D that are on the boundary of feasible set

Yhat = B( A(X,1), 1 );
DZ = D*Z;
signYhat = sign(Yhat);
signYhat(signYhat == 0) = 1;
Yt = Y .* signYhat;  % Y.*exp(1i*angle(ADZ));
residualData = Yt - Yhat;        % residual of data fitting term
residualApprox = X - DZ;        % residual of sparse approximation term

obj = 0.5*norm( residualData ,'fro')^2 + 0.5*mu*norm(residualApprox,'fro')^2 + lambdaL1normZ;

% compute gradients of smooth part of upper bound of objective
gradientX = - B( A(residualData,0), 0 ) + mu * residualApprox;    % HessianX*X - A'*Yt - mu*DZ;    % gradient of smooth part of objective w.r.t. X
gradientDpart = - residualApprox * Z';
gradientD = mu * gradientDpart;   % mu*(D*Z-X)*Z';   % gradient of smooth part of objective w.r.t. D
gradientZpart = - D' * residualApprox;
gradientZ = mu * gradientZpart;    % mu*D'*(D*Z-X);   % gradient of smooth part of objective w.r.t. Z

updateInd = true(3,1);  % indicate the variable blocks to be updated
stationaryInd = false(3,1); % indicate the variable blocks (X,D,Z) that have achieved stationarity with the given tolerance

% Initialize struct to record statistics:
recstats = (nargout > 4) || verb;
if( recstats )
    stat.runtime = zeros(maxIter+1,1);
    stat.objVal = zeros(maxIter+1,1);
    %     stat.distance = zeros(maxIter,1);
    stat.improve = zeros(maxIter,1);
    stat.stepsize = zeros(maxIter,1);
    stat.error_gradX = zeros(maxIter,1);
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
%     error_gradX = max( abs(gradientX(:)) );
%     error_gradX = norm(gradientX,'fro') * distanceX / abs(obj);
    error_gradX = norm(gradientX,'fro') / sqrt(numel(X)) / numel(Y);
    
    if lambda > eps
        subgradZ(:) = 0;
        subgradZ(patternZ) = gradientZ(patternZ) + lambda * sign(Z(patternZ));
        subgradZ(~patternZ) = max( abs(gradientZ(~patternZ))-lambda, 0);    % sign(gradientZ(~nnzInd)) .* max( abs(gradientZ(~nnzInd))-lambda, 0); phase can be dropped since we need only its norm
%         error_subgradZ = max( abs(subgradZ(:)) );
%         error_subgradZ = norm(subgradZ,'fro') * distanceZ / abs(obj);
        error_subgradZ = norm(subgradZ,'fro') / sqrt(numel(Z)) / numel(Y);
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

    error_subgrad = sqrt( (error_subgradD^2*numel(D) + error_subgradZ^2*numel(Z) + error_gradX*numel(X)) / (numel(D) + numel(Z) + numel(X)) );

    % evaluate and record estimation quality every 10 iterations
    clear sol;
    sol.D = D;
    sol.Z = Z;
    sol.X = X;
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

% distanceX = 1;
% distanceD = 1;
% distanceZ = 1;

% Main algorithm
for t = 1: maxIter
    timer = tic;
    
    if verb
        fprintf('iter %d:\n',t);
    end

    %% Select variable blocks to be updated
    switch blockRule
        case 0
            updateInd(:) = true;
        case 1
            % cyclic
            updateInd(1) = ~(mod(t,3) == 1);
            updateInd(2) = ~(mod(t,3) == 2);
            updateInd(3) = ~(mod(t,3) == 0);
            % random
%             updateInd(:) = true;
%             updateInd(randsample(3,blockRule)) = false;
        case 2
            updateInd(1) = (mod(t,3) == 1);
            updateInd(2) = (mod(t,3) == 2);
            updateInd(3) = (mod(t,3) == 0);
    end
%     updateInd = updateInd & (~stationaryInd);
    
    if any(updateInd)
        
        %% Compute descent direction for X
        if updateInd(1)
            % 1. Jointly, if there is no constraint on X
            %         S_X = X - HessianX \ gradientX;
            %         S_X = HessianX \ (A'*Yt + mu*DZ);
            % 2. Element-wise, if there is element-wise constraint on X
            %             S_X = X - gradientX ./ diag(HessianX);    % =  diag(1./(diag(A'*A)+mu)) * ( (diag(diag(A'*A)) - A'*A)*X + A'*Yt + mu*D*Z );
            solApproxX = X - gradientX ./ HessianXpart;
            deltaX = solApproxX - X;
            deltaYhat = B( A(deltaX,1), 1 );
%         else
%             %         S_X = X;
%             DeltaX(:) = 0;
        end
        
        %% Compute descent direction for D
        if updateInd(2)
            normZrows2 = sum(absZ.^2,2);  % column vector, squared norm of rows of Z
            nnzRowInd = (normZrows2 > eps);     % indicators of non-zero rows of Z
            nnzRow = sum(nnzRowInd);    % # of non-zero rows of Z
            
            % 1. For non-zero rows of Z
            if nnzRow > 0
                solApproxD(:,nnzRowInd) = D(:,nnzRowInd) - gradientDpart(:,nnzRowInd) ./ normZrows2(nnzRowInd).';
                %                 S_D(:,nnzRowInd) = D(:,nnzRowInd) - (D*Z-X)*Z(nnzRowInd,:)' * diag(1./diag(Z(nnzRowInd,:)*Z(nnzRowInd,:)'));
            end
            % 2. For zero rows of Z, any solution is optimal for corresponding columns of D
            % S_D(:,~nnzRowInd) = ones(N,P-nnzRow);
            % S_D(:,~nnzRowInd) = D(:,~nnzRowInd);    % keep those columns of D unchanged
            solApproxD(:,~nnzRowInd) = randn(N,P-nnzRow) +1i*randn(N,P-nnzRow);
            
            normSDcols2 = sum( abs(solApproxD).^2 );   % squared l2-norm of columns of S_D
            ind = (nnzRowInd.' & normSDcols2 <= 1); % indicators of columns of S_D that don't need to be normalized
            solApproxD(:,~ind) = solApproxD(:,~ind) ./ sqrt( normSDcols2(~ind) );    % normalize columns of D
            
            deltaD = solApproxD - D;
%         else
%             %         S_D = D;
%             DeltaD(:) = 0;
        end
        
        
        
        %% Compute descent direction for Z
        if updateInd(3)
            if debias % debiasing mode
                solApproxZ = Z - gradientZpart ./ normDcols2.';
%                 S_Z(~patternZ) = 0;   % this step is not necessary because S_Z will not be used in the stepsize computation
                deltaZ(:) = 0;
                deltaZ(patternZ) = solApproxZ(patternZ) - Z (patternZ);
            else
                solApproxZ = soft_threshold( lambda_mu, normDcols2.' .* Z - gradientZpart ) ./ normDcols2.';
                deltaZ = solApproxZ - Z;
            end
%         else
%             %         S_Z = Z;
%             DeltaZ(:) = 0;
        end
        
%         distance = norm([DeltaX(:); DeltaD(:); DeltaZ(:)]);
        
        %% Compute stepsize
        % ---------------------------------------------------------------------
        % coeficients of fourth-order polynomial
        % 1/4*coef(1) * gamma^4 + 1/3*coef(2) * gamma^3 + 1/2*coef(3) * gamma^2
        % +coef(4) * gamma + coef(5)
        % ---------------------------------------------------------------------
        coef(:) = 0;
        if all(updateInd(2:3))
            % Both D and Z are to be updated -> line search function is a fourth-order polynomial
            M1part = deltaD * Z + D * deltaZ;
            M2 = deltaD * deltaZ;
            if lambda > eps
                m4 = lambda * sum(abs(solApproxZ(:))) - lambdaL1normZ;    % lambda * ( sum(abs(S_Z(:))) - sum(abs(Z(:))) );
            else
                m4 = 0;
            end
            if updateInd(1)
                % X is to be updated
                M1 = M1part - deltaX;
                
                
                coef(1) = 2 * mu * norm(M2,'fro')^2;
                coef(2) = 3 * mu * real( M1(:)' * M2(:) );
                coef(3) = mu * ( -2*real( residualApprox(:)' * M2(:) ) + norm(M1,'fro')^2 ) + norm(deltaYhat,'fro')^2;
                coef(4) = - mu * real( residualApprox(:)' * M1(:) ) + m4 - real( residualData(:)' * deltaYhat(:) );     % must be < 0 as a descent direction
                %     coef(5) = 0.5*mu*norm(Res2,'fro')^2 + 0.5*norm(Res1,'fro')^2;   % constant has no effect on min
            else
                % X is not to be updated
                M1 = M1part;
                
                coef(1) = 2 * mu * norm(M2,'fro')^2;
                coef(2) = 3 * mu * real( M1(:)' * M2(:) );
                coef(3) = mu * ( -2*real( residualApprox(:)' * M2(:) ) + norm(M1,'fro')^2 );
                coef(4) = - mu*real( residualApprox(:)' * M1(:) ) + m4;     % must be < 0 as a descent direction
            end
            
            if coef(4) >= 0
                if verb
                    fprintf('Not descent direction!\n');
                end
                t = t-1;
                break;
            end
            
            if coef(1) > eps
                % The line search function is a fourth-order polynomial
                rts = realcubicroots(coef(1:4)); % real roots including multiplicity
                rts = rts( (rts)<=1 & (rts)>=0);    % real roots in [0,1]
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
                    polyValRTS = polyval(coef./[4:-1:1,1]',rts1);
                    [~,minIdx] = min(polyValRTS);
                    stepsize = rts1(minIdx);
                end
            elseif coef(3) > eps % line search function is quadratic
                stepsize = min( max( -coef(4) / coef(3), 0 ), 1 );
            else % line search function is linearly decreasing
                stepsize = 1;
            end
            
        elseif updateInd(3) && all(~updateInd(1:2))
            % Only Z is to be updated
            % The line search function is a quadratic polynomial
            % 1/2*coef(3) * gamma^2 - coef(4) * gamma + const.
            M1part = D * deltaZ;
            M1 = M1part;
            % 1. Joint line search for Z
%             coef(3) = norm(M1,'fro')^2;
%             coef(4) = real( Res2(:)' * M1(:) ) - (lambda * sum(abs(S_Z(:))) - lambdaL1normZ)/mu;  % must be > 0 as a descent direction
%             if coef(4) < 0
%                 if verb
%                     fprintf('Not descent direction!');
%                 end
%                 t = t-1;
%                 break;
%             end
%             if coef(3) > 0
%                 stepsize = min( max( coef(4) / coef(3), 0 ), 1 );
%             else
%                 stepsize = 0;
%             end
            % 2. Independent line search for each column of Z
            coefZ1 = sum(abs(M1).^2);   % row vector
            coefZ2 = real(dot(residualApprox,M1));   % row vector
            if lambda > eps
                coefZ2 = coefZ2  - lambda_mu * sum(abs(solApproxZ)-absZ);
            end
            if any(coefZ2 < 0)
                if verb
                    fprintf('Not descent direction!\n');
                end
                t = t-1;
                break;
            end
            nnzCoefZ1 = (coefZ1 > eps);  % indicates nonzeros of coefZ1
            stepsize = ones(1,I);
            stepsize(nnzCoefZ1) = min( max( coefZ2(nnzCoefZ1) ./ coefZ1(nnzCoefZ1), 0 ), 1 );
        else
            % The line search function is a quadratic polynomial
            % 1/2*coef(3) * gamma^2 - coef(4) * gamma + const.
            if updateInd(2)
                % D is to be updated
                M1part = deltaD * Z;
                m4 = 0;
            elseif updateInd(3)
                % Z is to be updated
                M1part = D * deltaZ;
                if lambda > eps
                    m4 = lambda * sum(abs(solApproxZ(:))) - lambdaL1normZ;    % lambda * ( sum(abs(S_Z(:))) - sum(abs(Z(:))) );
                else
                    m4 = 0;
                end
            else
                % Neith D nor Z
                M1part = zeros(N,I);
                m4 = 0;
            end
            
            if updateInd(1)
                % X is to be updated
                M1 = M1part - deltaX;
                
                coef(3) = mu * norm(M1,'fro')^2 + norm(deltaYhat,'fro')^2;
                coef(4) = mu*real( residualApprox(:)' * M1(:) ) + real( residualData(:)' * deltaYhat(:) ) - m4;     % must be > 0 as a descent direction
            else
                % X is not to be updated
                M1 = M1part;
                
                coef(3) = mu * norm(M1,'fro')^2;
                coef(4) = mu*real( residualApprox(:)' * M1(:) ) - m4;     % must be > 0 as a descent direction
            end
            
            if coef(4) < 0
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
        end
        
        %% updateInd variables and intermediate variables
        if updateInd(1)
            % X is to be updated
            X = X + stepsize * deltaX;
            Yhat = Yhat + stepsize * deltaYhat;
            signYhat = sign(Yhat);
            signYhat(signYhat == 0) = 1;
            Yt = Y .* signYhat;  % Y.*exp(1i*angle(ADZ));
            residualData = Yt - Yhat;    % residual of the data fitting term
%             distanceX = stepsize * norm(deltaX,'fro');
        end
        
        if all(updateInd(2:3))
            % Both D and Z are to be updated
            D = D + stepsize * deltaD;
            if debias
                Z(patternZ) = Z(patternZ) + stepsize * deltaZ(patternZ);
%                 distanceZ = stepsize * norm(deltaZ(patternZ));
            else
                Z = Z + stepsize * deltaZ;
%                 distanceZ = stepsize * norm(deltaZ);
            end
            
            % Update intermediate variables
            normDcols2 = sum(abs(D).^2);    % squared l2-norm of columns of D
            boundaryDInd = (normDcols2 >= 1 - tol_rational);% indices of columns of D that are on the boundary of feasible set
            
            absZ = abs(Z);
            if ~debias
                patternZ = (absZ > eps);
                if lambda > eps
                    lambdaL1normZ = lambda * sum(absZ(:));
                end
            end
            
            DZ = DZ + stepsize * M1part + stepsize^2 * M2;
        else
            if updateInd(2)
                % D is to be updated and Z is NOT to be updated
                D = D + stepsize * deltaD;
                normDcols2 = sum(abs(D).^2);    % squared l2-norm of columns of D
                boundaryDInd = (normDcols2 >= 1 - tol_rational);% indices of columns of D that are on the boundary of feasible set
                
                DZ = DZ + stepsize * M1part;
            end
            if updateInd(3)
                % Z is to be updated and D is NOT to be updated
                updateZ = stepsize .* deltaZ;
                Z = Z + updateZ;
%                 distanceZ = norm(updateZ,'fro');
                absZ = abs(Z);
                if ~debias
                    patternZ = (absZ > eps);
                    if lambda > eps
                        lambdaL1normZ = lambda * sum(absZ(:));
                    end
                end
                
                DZ = DZ + stepsize .* M1part;
            end
        end
        
        if any(normDcols2 == 0)
            fprintf('STOP!');
        end
        
        %% updateInd intermediate variables and objective function value
        residualApprox = X - DZ;     % residual of the dictionary learing term        
        obj = 0.5*norm( residualData ,'fro')^2 + 0.5*mu*norm(residualApprox,'fro')^2 + lambdaL1normZ;
        
        % compute gradients of smooth part of upper bound of objective
        gradientX = - B( A(residualData,0), 0 ) + mu * residualApprox;    % HessianX*X - A'*Yt - mu*DZ;    % gradient of smooth part of objective w.r.t. X
        gradientDpart = - residualApprox * Z';
        gradientD = mu * gradientDpart;   % mu*(D*Z-X)*Z';   % gradient of smooth part of objective w.r.t. D
        gradientZpart = - D' * residualApprox;
        gradientZ = mu * gradientZpart;    % mu*D'*(D*Z-X);   % gradient of smooth part of objective w.r.t. Z
        
    else
        % No updateInd in this iteration
%         distance = 0;
        stepsize = 0;
    end
    
    runtime = toc(timer);
    
    %% Evaluate optimality of current point
    % minimum-norm sub-gradient
%     error_gradX = max( abs(gradientX(:)) );
%     error_gradX = norm(gradientX,'fro') * distanceX / abs(obj);
    error_gradX = norm(gradientX,'fro') / sqrt(numel(X)) / numel(Y);
    
    if lambda > eps
        subgradZ(:) = 0;
        subgradZ(patternZ) = gradientZ(patternZ) + lambda * sign(Z(patternZ));
        subgradZ(~patternZ) = max( abs(gradientZ(~patternZ))-lambda, 0);    % sign(gradientZ(~nnzInd)) .* max( abs(gradientZ(~nnzInd))-lambda, 0); phase can be dropped since we need only its norm
%         error_subgradZ = max( abs(subgradZ(:)) );
%         error_subgradZ = norm(subgradZ,'fro') * distanceZ / abs(obj);
        error_subgradZ = norm(subgradZ,'fro') / sqrt(numel(Z)) / numel(Y);
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
    
    % record statistics
    if recstats
        stat.runtime(t+1) = stat.runtime(t) + runtime;
        stat.objVal(t+1) = obj;
        %         stat.distance(t) = distance;
        stat.stepsize(t) = max(stepsize);
        stat.improve(t) = (stat.objVal(t)-stat.objVal(t+1)) / stat.objVal(t);
        if debias
            stat.sparsityZ(t+1) = stat.sparsityZ(t);
            stat.rmse_sparsityZ(t+1) = stat.rmse_sparsityZ(t);
        else
            stat.sparsityZ(t+1) = sum(patternZ(:)) / I;
            stat.rmse_sparsityZ(t+1) = sqrt( mean( abs( sum(absZ > eps) - L ) .^2 ) );
        end
        stat.error_gradX(t) = error_gradX;
        stat.error_subgradZ(t) = error_subgradZ;
        stat.error_subgradD(t) = error_subgradD;
    end

    if trackConvergence
        % evaluate and record estimation quality every 10 iterations
        if mod(t,10) == 0
            clear sol;
            sol.D = D;
            sol.Z = Z;
            sol.X = X;
            quality = eval_quality(sol,data);
            stat.quality.iterArray = [stat.quality.iterArray; t];
            stat.quality.timeArray = [stat.quality.timeArray; stat.runtime(t+1)];
            stat.quality.FmeasureZ = [stat.quality.FmeasureZ; quality.FmeasureZ];
            stat.quality.errorD = [stat.quality.errorD; quality.errorD];
            stat.quality.errorZ = [stat.quality.errorZ; quality.errorZ];
            stat.quality.errorDZ = [stat.quality.errorDZ; quality.errorDZ];
            stat.quality.errorX = [stat.quality.errorX; quality.errorX];
            error_subgrad = sqrt( (error_subgradD^2*numel(D) + error_subgradZ^2*numel(Z) + error_gradX*numel(X)) / (numel(D) + numel(Z) + numel(X)) );
            stat.quality.error_subgrad = [stat.quality.error_subgrad; error_subgrad];
            stat.quality.objVal = [stat.quality.objVal; obj];
        end
    end
    
    if verb
        fprintf('\tUpdate %s\n',var_str(updateInd));
        fprintf('\tvalue = %f, step-size = %f, improve = %e,\n\tgradientX = %e, subgradD = %e, subgradZ = %e,\n\tavg. code sparsity = %f, time = %f s.\n',...
            obj,max(stepsize),stat.improve(t),error_gradX, error_subgradD, error_subgradZ, stat.sparsityZ(t+1), stat.runtime(t));
    end

    %% Stopping criterion
    stationaryInd(1) = ( error_gradX <= tol );
    stationaryInd(2) = ( error_subgradD <= tol );
    stationaryInd(3) = ( error_subgradZ <= tol );
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
    stat.Niter = t;     % # iterations
    stat.optVal = stat.objVal(t+1);  % optimal value
    if t < maxIter
        if t > 0
            stat.runtime(t+2:end) = stat.runtime(t+1);
        end
        stat.objVal(t+2:end) = stat.optVal;
        stat.error_gradX(t+1:end) = stat.error_gradX(t);
        stat.error_subgradD(t+1:end) = stat.error_subgradD(t);
        stat.error_subgradZ(t+1:end) = stat.error_subgradZ(t);
        stat.sparsityZ(t+2:end) = stat.sparsityZ(t+1);
        stat.rmse_sparsityZ(t+2:end) = stat.rmse_sparsityZ(t+1);
    end
end

end

