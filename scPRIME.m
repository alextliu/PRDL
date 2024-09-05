function [X,D,Z,patternZ,stat] = scPRIME(data,params,X0,D0,Z0)

% SCPRIME This function solves the following phase retrieval problem with
% dictionary learing problem:
%
% min  0.5*||Y-abs(A*X*B)||^2_F + 0.5*mu*||X-D*Z||^2_F
% s.t. D contains unit columns
%
% Using BSUM
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
            fprintf('Invalid data: debiasing phase with empty support for Z. Aborting.\n');
        end
        return;
    end
    Z(~patternZ) = 0;
end

% Allocate memory
subgradZ = zeros(P,I);

% Initial intermediate variables, gradient and objective function
lambda_mu = lambda / mu;
% hessianXpart = sum( abs( A(eye(N),1) ).^2 ).' * sum( abs( B(eye(I),1) ).^2, 2 ).' + mu;    % second-order partial derivatives of smooth part of objective w.r.t. each entry of X
hessianXmajor = ( svds(A(eye(N),1),1) * svds(B(eye(I),1),1) )^2 + mu;

Yhat = B( A(X,1), 1 );
signYhat = sign(Yhat);
signYhat(signYhat == 0) = 1;
Yt = Y .* signYhat;
residualData = Yt - Yhat;        % residual of data fitting term
residualApprox = X - D*Z;        % residual of sparse approximation term

obj = 0.5*norm( residualData ,'fro')^2 + 0.5*mu*norm(residualApprox,'fro')^2 + lambdaL1normZ;

% compute gradients of smooth part of upper bound of objective
gradientX = - B( A(residualData,0), 0 ) + mu * residualApprox;    % HessianX*X - A'*Yt - mu*DZ;    % gradient of smooth part of objective w.r.t. X
gradientD = - mu * residualApprox * Z';   % mu*(D*Z-X)*Z';   % gradient of smooth part of objective w.r.t. D
gradientZ = - mu * D' * residualApprox;    % mu*D'*(D*Z-X);   % gradient of smooth part of objective w.r.t. Z

stationaryInd = false(3,1); % indicate the variable blocks (X,D,Z) that have achieved stationarity with the given tolerance

% Initialize struct to record statistics:
recstats = (nargout > 4) || verb;
if( recstats )
    stat.runtime = zeros(maxIter+1,1);
    stat.objVal = zeros(maxIter+1,1);
    %     stat.distance = zeros(maxIter,1);
    stat.improve = zeros(maxIter,1);
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
    
    subgradD = gradientD + sqrt( sum( abs(gradientD).^2, 1 ) ) .* D;
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
    
    %% Update X
    updateX = - gradientX ./ hessianXmajor;
    X = X + updateX;
%     distanceX = norm(updateX,'fro');
    
    %% Update D
    Dold = D;
    XZh = X*Z';
    ZZh = Z*Z';
    for p = 1:P
        if ZZh(p,p) > 0
            newdp = (XZh(:,p)-D*ZZh(:,p))/ZZh(p,p) + D(:,p);
        else
            newdp = randn(N,1);
        end
        normNewdp = norm(newdp);
        D(:,p) = newdp / max(1,normNewdp);
    end
%     distanceD = norm(D-Dold,'fro');
    
    %% Update Z
    Zold = Z;
    if debias % debiasing mode
        Z = Z - D'*(D*Z-X) ./ P;
        Z(~patternZ) = 0;
    else
        Z = soft_threshold( lambda_mu, P .* Z - D'*(D*Z-X) ) ./ P;
    end
%     distanceZ = norm(Z-Zold,'fro');
    
    %% update intermediate variables and objective function value
    Yhat = B( A(X,1), 1 );
    signYhat = sign(Yhat);
    signYhat(signYhat == 0) = 1;
    Yt = Y .* signYhat;  % Y.*exp(1i*angle(ADZ));
    residualData = Yt - Yhat;        % residual of data fitting term
    residualApprox = X - D*Z;        % residual of sparse approximation term
    absZ = abs(Z);
    if lambda > eps
        lambdaL1normZ = lambda * sum(absZ(:));
    else
        lambdaL1normZ = 0;
    end
    
    oldobj = obj;
    obj = 0.5*norm( residualData ,'fro')^2 + 0.5*mu*norm(residualApprox,'fro')^2 + lambdaL1normZ;
    improve = (oldobj - obj)/oldobj;

    % compute gradients of smooth part of upper bound of objective
    gradientX = - B( A(residualData,0), 0 ) + mu * residualApprox;    % HessianX*X - A'*Yt - mu*DZ;    % gradient of smooth part of objective w.r.t. X
    
    runtime = toc(timer);
    
    gradientD = - mu * residualApprox * Z';   % mu*(D*Z-X)*Z';   % gradient of smooth part of objective w.r.t. D
    gradientZ = - mu * D' * residualApprox;    % mu*D'*(D*Z-X);   % gradient of smooth part of objective w.r.t. Z

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
    
    subgradD = gradientD + sqrt( sum( abs(gradientD).^2, 1 ) ) .* D;
%     error_subgradD = max( abs(subgradD(:)) );
%     error_subgradD = norm(subgradD,'fro') * distanceD / abs(obj);
    error_subgradD = norm(subgradD,'fro') / sqrt(numel(D)) / numel(Y);
    
    % record statistics
    if recstats
        stat.runtime(t+1) = stat.runtime(t) + runtime;
        stat.objVal(t+1) = obj;
        %         stat.distance(t) = distance;
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
        fprintf('\tvalue = %f, improve = %e,\n\tgradientX = %e, subgradD = %e, subgradZ = %e,\n\tavg. code sparsity = %f, time = %f s.\n',...
            obj,stat.improve(t),error_gradX, error_subgradD, error_subgradZ, stat.sparsityZ(t+1), stat.runtime(t));
%         fprintf('\tdistanceX = %f, distanceD = %f, distanceZ = %f\n',distanceX,distanceD,distanceZ);
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
    
    if improve < 0
        if verb
            fprintf('Not descent update. Aborting.');
        end
        t = t-1;
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

