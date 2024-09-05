function quality = eval_quality(sol,data)
% EVAL_QUALITY evaluates estimation quality

% PRECISION = 1e-9;
spatialMixerType = data.spatialMixerType;
temporalMixerType = data.temporalMixerType;

[N,P] = size(data.D);
I = size(data.Z,2);
if ~isfield(data,'X')
    data.X = data.D * data.Z;
end

D = sol.D;
D_true = data.D;
X_true = data.X;
Z = sol.Z;
Z_true = data.Z;

% To resolve the permutation ambiguity, finds the permutation of columns of
% estimated D that best matches the groundtruth using abs of normalized
% cross correlation between columns of the estimated D and the groundtruth D
xCorrD = abs( D'*D_true ./ sqrt(sum(abs(D_true).^2)) ./ sqrt(sum(abs(D).^2)).' );
permutation = zeros(P,1);
tmp = xCorrD;
for p = 1:P
    [~,ind] = max(tmp(:));
    [ii,jj] = ind2sub([P,P],ind);
    permutation(jj) = ii;
    tmp(ii,:) = -1;
    tmp(:,jj) = -1;
end

D_permut = D(:,permutation);
Z_permut = Z(permutation,:);
X_est = D_permut * Z_permut;    % reconstruct X as D*Z

% Find the optimal complex scale for each columns of D that minimizes the
% squared error between estimate and groundtruth
dictScale = dot(D_permut,D_true) ./ sum(abs(D_permut).^2);
D_est = D_permut .* dictScale;  % apply the optimal complex scale on columns of D
quality.errorD = norm(D_true - D_est,'fro') / norm(D_true,'fro'); % minimium normalized root squared error of D

% Evaluate estimation quality of support of Z
if ~isfield(data,'patternZ')
    data.patternZ = (abs(data.Z)>eps);  % true sparsity pattern of Z
end
patternZ_est = (abs(Z_permut)>eps);    % estimated sparsity pattern of Z
nnzZ_est = sum(patternZ_est(:));
quality.sparsityZ = nnzZ_est / size(Z_permut,2);     % average sparsity of each column of Z
quality.rmse_sparsityZ = sqrt( mean( abs( sum(patternZ_est) - data.L ) .^2 ) ); % root-mean-squared-error of sparsity level of each column compared to the groundtruth value L
if nnzZ_est == 0 % estimated Z is empty
    quality.precisionZ = 0;
    quality.recallZ = 0;
    quality.FmeasureZ = 0;
else
    tp = sum(patternZ_est(:) & data.patternZ(:));    % true positive
    quality.precisionZ = tp / nnzZ_est;
    quality.recallZ = tp / sum(data.patternZ(:));
    quality.FmeasureZ = harmmean([quality.precisionZ, quality.recallZ]);
end

% Compute estimation error of product D*Z
normX_true = norm(X_true,'fro');
if strcmpi(temporalMixerType,'no') || contains(temporalMixerType,'stcdp','IgnoreCase',true)  % no temporal mixer or short-time coded diffraction pattern temporal mixer
    if contains(temporalMixerType,'stcdp','IgnoreCase',true)   % short-time coded diffraction pattern temporal mixer
        X_true = reshape(X_true,N*data.len_stcdp,[]);
        X_est = reshape(X_est,N*data.len_stcdp,[]);
    end
    dotX = dot(X_true, X_est);
    nnzInd = ( abs(dotX) > eps );
    phaseShiftX = ones(size(dotX));
    phaseShiftX(nnzInd) = abs(dotX(nnzInd)) ./ dotX(nnzInd);
    quality.errorDZ = sqrt( norm(X_est,'fro')^2 + normX_true^2 - 2*sum(abs(dotX)) ) / normX_true;
    %     quality.errorDZ = norm( abs(X_est) - abs(X_true), 'fro') / normX_true;
else
    dotX = X_true(:)' * X_est(:);
    if abs(dotX) > eps
        phaseShiftX = abs(dotX) / dotX;
    else
        phaseShiftX = 1;
    end
    quality.errorDZ = sqrt( norm(X_est,'fro')^2 + normX_true^2 - 2*abs(dotX) ) / normX_true;
end

% Compute estimation error of Z
normZ_true = norm(Z_true,'fro');
if contains(temporalMixerType,'stcdp','IgnoreCase',true)   % short-time coded diffraction pattern temporal mixer
    Z_rotate = reshape( reshape(Z_permut,P*data.len_stcdp,[]) .* phaseShiftX, [P,I] ); % rotate Z with the same phase shift as for X
else
    Z_rotate = Z_permut .* phaseShiftX; % rotate Z with the same phase shift as for X
end
Z_est = Z_rotate;
nnzRowInd = sum(patternZ_est,2) > 0; % indices of nonzero rows in Z_est
scaleZ = dot(Z_rotate(nnzRowInd,:),Z_true(nnzRowInd,:),2) ./ sum(abs(Z_rotate(nnzRowInd,:)).^2,2);
Z_est(nnzRowInd,:) = Z_rotate(nnzRowInd,:) .* scaleZ;
quality.errorZ = norm(Z_true - Z_est, 'fro') / normZ_true;  % minimum normalized root squared error of Z
quality.errorAbsZ = norm( abs(Z_true) - abs(Z_est), 'fro' ) / normZ_true;   % normalized root squared error of abs of Z

% For the formulation with variabe X, evaluate also the solution of X
if isfield(sol,'X')
    X_est = sol.X;
    if strcmpi(temporalMixerType,'no') || contains(temporalMixerType,'stcdp','IgnoreCase',true)  % no temporal mixer or short-time coded diffraction pattern temporal mixer
        if contains(temporalMixerType,'stcdp','IgnoreCase',true)   % short-time coded diffraction pattern temporal mixer
            %             X_true = reshape(X_true,N*data.len_stcdp,[]);
            X_est = reshape(X_est,N*data.len_stcdp,[]);
        end
        quality.errorX = sqrt( norm(X_est,'fro')^2 + normX_true^2 - 2*sum( abs( dot(X_true, X_est) ) ) ) / normX_true;
    else
        quality.errorX = sqrt( norm(X_est,'fro')^2 + normX_true^2 - 2*sum( abs( X_true(:)' * X_est(:) ) ) ) / normX_true;
    end
else
    quality.errorX = quality.errorDZ;
end

end

