% This simulation varys the number of random initializations


addpath('tools/');
% clc;
clear all;
% close all;

% diary off;
% diary(['diary_',datestr(date,'yyyy-mm-dd'),'.txt']);

rng('default');

format compact;

% PRECISION = 1e-9;

% specify spatial mixer type
spatialMixerTypeList = {
    'no';
    'gauss';
    'random_phase';     % 3: each entry of spatial mixing matrix has unit norm and uniform random phase
    'cdp_binary';       % 4: coded diffraction pattern with binary mask
    'cdp_ternary';      % 5: coded diffraction pattern with ternary mask
    'cdp_complex';      % 6: coded diffraction pattern with complex mask
    'dft';
    };
spatialMixerTypeNumber = 2;
spatialMixerType = spatialMixerTypeList{spatialMixerTypeNumber};

% specify temporal mixer type
temporalMixerTypeList = {
    'no';
    'gauss';
    'stcdp_binary'      % 3: short-time coded diffraction pattern with binary mask
    'stcdp_ternary'     % 4: short-time coded diffraction pattern with ternary mask
    'stcdp_complex'     % 5: short-time coded diffraction pattern with complex mask
    'stft'              % 6: short-time fourier transform
    };
temporalMixerTypeNumber = 1;
temporalMixerType = temporalMixerTypeList{temporalMixerTypeNumber};

% specify channel type
channelTypeList = {
    'gauss';        % 1: each entry of ground-truth D i.i.d. follows CN(0,1)
    'los';          % 2: D is LOS of ULA
    };
channelTypeNumber = 1;
channelType = channelTypeList{channelTypeNumber};

% specify measurement model
measurementModelList = {
    'gauss';
    'pois';     % disabled
    };
measurementModelNumber = 1;
measurementModel = measurementModelList{measurementModelNumber};

N = 16; %100;       % no. of Rx antennas
I = 16*N;     % no. of snapshots
P = N/2;   % no. of users
L = 2;%0.1*P;     % sparsity level

% parameters for short-time cdp temporal mixing
len_stcdp = I;  % length of short-time coded diffraction pattern

% parameters for STFT temporal mixing
stft_window_len = I/2;    % length of sliding window in STFT; rectangular window
stft_hopsize = I/4;     % hop size in STFT
stft_fft_len = I;       % no. of fft measurements in stft

% sptial and temporal oversampling rates
spatialSampleRate = 4;
temporalSampleRate = 4;
if strcmpi(spatialMixerType,'no')
    spatialSampleRate = 1;
end
if strcmpi(temporalMixerType,'no')
    temporalSampleRate = 1;
elseif strcmpi(temporalMixerType,'stft')
    temporalSampleRate = ( (I + stft_window_len) / stft_hopsize - 1 ) * stft_fft_len / I;
end

% measurement size
M1 = spatialSampleRate * N;
M2 = temporalSampleRate * I;

SNR = 15;

% mu_factor = 10^1.5;     % for sampling rate 2 x 2
mu_factor = 1;

% regularization parameters for the case without temporal mixing
lambda_factor = 0.75^14;
rho_factor = 0.75^14;

% regularization parameters for the case with temporal mixing
% lambda_factor = 0.75^21;
% rho_factor = 0.75^23;

nInstance = 1;     % no. of random instances

initialSet = 1; % [1:5,10,20,30,50];             % test set of no. of random initializations for each instance
initialSize = length(initialSet); % size of set of no. of random initializations to test
maxInitial = max(initialSet);     % max no. of random initializations to test

% set of testing algorithms
algoSet = {
    @FUN_PRDLsca;
    @FUN_PRDLscaX;
    @scPRIME;
    };
nAlgo = length(algoSet);    % no. of algorithms to test
blockRule = zeros(nAlgo,1);  % indicates block update rule, represents no. of blocks NOT to update in each iteration
debiasingOn = 0;

% regularization parameters
approxParamSet = ones(nAlgo,1);
if length(approxParamSet) < nAlgo
    fprintf('Inconsistent simulation parameters. Aborting.\n');
    return;
end
sparsityParamSet = [
    0.75^17;
    0.75^14;
    0.75^14
    ];
if length(sparsityParamSet) < nAlgo
    fprintf('Inconsistent simulation parameters. Aborting.\n');
    return;
end

% set parameters for algorithms
params.maxIter = 2000;
params.tol = 1e-4;
params.verb = 0;   % show information in each iteration if verb = 1
params.N = N;
params.I = I;
params.P = P;
params.L = L;

% varying parameters
paramList = {
    'P';    % 1: no. of users P
    'L';    % 2: sparsity level L
    'M1';   % 3: no. of spatial measurements M1
    'SNR';  % 4: SNR in dB for Gaussian observation model
    'I';    % 5: no. of time-slots I
    };
paramSetList = {
    N; % [N/2 N N*2];        % test set of P
    [2];            % test set of L
    [2 4 8] * N;        % test set of M1
    [10,15,20,40,inf];  % test set of SNR
    [4,8,16,32,64] * N; % test set of I
    };
paramNumber1 = 1;
paramNumber2 = 2;
paramName1 = paramList{paramNumber1};
paramName2 = paramList{paramNumber2};
paramSet1 = paramSetList{paramNumber1};
paramSet2 = paramSetList{paramNumber2};
nParam1 = length(paramSet1);
nParam2 = length(paramSet2);

disp(datetime);
fprintf('Simulation Parameters:\n');
fprintf('Varying parameters: %s = %s, %s = %s, no. of random initializations = %s\n',paramName1,mat2str(paramSet1),paramName2,mat2str(paramSet2),mat2str(initialSet));
fprintf('P = %d, L = %d, N = %d, I = %d, M1 = %d, M2 = %d, SNR = %d dB\n',P,L,N,I,M1,M2,SNR);
fprintf('lambda factor = 0.75^%.0f, ', log(lambda_factor)/log(0.75) );
fprintf('mu factor = %.2e, rho factor = 0.75^%.0f\n', mu_factor, log(rho_factor)/log(0.75));
fprintf('Tolerance = %1.0e, max no. iterations = %d, Monte-Carlo runs = %d\n',params.tol,params.maxIter,nInstance);
fprintf('Channel type: %s\n', channelType);
fprintf('Measurement type: spatial mixer = %s, temporal mixer = %s\n', spatialMixerType, temporalMixerType);
fprintf('spatial oversampling rate = %.2f, temporal oversampling rate = %.2f\n',spatialSampleRate,temporalSampleRate);
if contains(temporalMixerType,'stft','IgnoreCase',true)
    fprintf('Window length = %d, hop size = %d, no. FFT samples = %d\n\n',stft_window_len,stft_hopsize,stft_fft_len);
end
disp('Algorithms:');
disp(algoSet);

%% Data to record
sparsityZ = cell(nParam1,nParam2,nAlgo);
rmse_sparsityZ = cell(nParam1,nParam2,nAlgo);

precisionZ = cell(nParam1,nParam2,nAlgo);
recallZ = cell(nParam1,nParam2,nAlgo);
FmeasureZ = cell(nParam1,nParam2,nAlgo);

errorD = cell(nParam1,nParam2,nAlgo);
errorZ = cell(nParam1,nParam2,nAlgo);
errorAbsZ = cell(nParam1,nParam2,nAlgo);

errorX = cell(nParam1,nParam2,nAlgo);
errorDZ = cell(nParam1,nParam2,nAlgo);

Niter = cell(nParam1,nParam2,nAlgo);
runtime = cell(nParam1,nParam2,nAlgo);

optInitial = cell(nParam1,nParam2,nAlgo);

best.sparsityZ = cell(nParam1,nParam2,nAlgo);
best.rmse_sparsityZ = cell(nParam1,nParam2,nAlgo);

best.precisionZ = cell(nParam1,nParam2,nAlgo);
best.recallZ = cell(nParam1,nParam2,nAlgo);
best.FmeasureZ = cell(nParam1,nParam2,nAlgo);

best.errorD = cell(nParam1,nParam2,nAlgo);
best.errorZ = cell(nParam1,nParam2,nAlgo);
best.errorAbsZ = cell(nParam1,nParam2,nAlgo);
best.errorX = cell(nParam1,nParam2,nAlgo);
best.errorDZ = cell(nParam1,nParam2,nAlgo);

best.debiasedErrorD = cell(nParam1,nParam2,nAlgo);
best.debiasedErrorZ = cell(nParam1,nParam2,nAlgo);
best.debiasedErrorAbsZ = cell(nParam1,nParam2,nAlgo);
best.debiasedErrorX = cell(nParam1,nParam2,nAlgo);
best.debiasedErrorDZ = cell(nParam1,nParam2,nAlgo);

best.Niter = cell(nParam1,nParam2,nAlgo);
best.runtime = cell(nParam1,nParam2,nAlgo);


%% Average results
avgResults.pattern.sparsityZ = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.pattern.rmse_sparsityZ = zeros(nParam1,nParam2,nAlgo,initialSize);

avgResults.pattern.precisionZ = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.pattern.recallZ = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.pattern.FmeasureZ = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.pattern.FmeasureZ_ref = zeros(nParam1,nParam2);

avgResults.error.D = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.error.Z = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.error.AbsZ = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.error.X = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.error.DZ = zeros(nParam1,nParam2,nAlgo,initialSize);

avgResults.debiasedError.D = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.debiasedError.Z = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.debiasedError.AbsZ = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.debiasedError.X = zeros(nParam1,nParam2,nAlgo,initialSize);
avgResults.debiasedError.DZ = zeros(nParam1,nParam2,nAlgo,initialSize);

avgResults.complexity.Niter = zeros(nParam1,nParam2,nAlgo,2,initialSize);
avgResults.complexity.runtime = zeros(nParam1,nParam2,nAlgo,2,initialSize);

%% Simulation starts
for iParam1 = 1:nParam1
    switch paramName1
        case 'P'
            P = paramSet1(iParam1);
            params.P = P;
        case 'L'
            L = paramSet1(iParam1);
            params.L = L;
        case 'M1'
            M1 = paramSet1(iParam1);
            spatialSampleRate = M1/N;
        case 'SNR'
            SNR = paramSet1(iParam1);
        case 'I'
            I = paramSet1(iParam1);
            params.I = I;
            M2 = temporalSampleRate * I;
    end
    
    for iParam2 = 1:nParam2
        switch paramName2
            case 'P'
                P = paramSet2(iParam2);
                params.P = P;
            case 'L'
                L = paramSet2(iParam2);
                params.L = L;
            case 'M1'
                M1 = paramSet2(iParam2);
                spatialSampleRate = M1/N;
            case 'SNR'
                SNR = paramSet2(iParam2);
            case 'I'
                I = paramSet2(iParam2);
                params.I = I;
                M2 = temporalSampleRate * I;
        end
        
        fprintf('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n');
        fprintf('>>>>>> P = %2d, L = %2d, N = %2d, I = %d, M1 = %d, M2 = %d, SNR = %2d dB <<<<<<\n',P,L,N,I,M1,M2,SNR);
        fprintf('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n');
        
        % initialize arrays to record results
        for iAlgo = 1:nAlgo
            
            errorD{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            
            sparsityZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            rmse_sparsityZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            precisionZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            recallZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            FmeasureZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            errorZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            errorAbsZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            
            errorX{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            errorDZ{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            
            Niter{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);
            runtime{iParam1,iParam2,iAlgo} = zeros(nInstance,maxInitial);

            optInitial{iParam1,iParam2,iAlgo} = zeros(nInstance,initialSize);
            
            best.sparsityZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.rmse_sparsityZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            
            best.precisionZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.recallZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.FmeasureZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            
            best.errorD{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.errorZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.errorAbsZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.errorX{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.errorDZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            
            best.debiasedErrorD{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.debiasedErrorZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.debiasedErrorAbsZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.debiasedErrorX{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            best.debiasedErrorDZ{iParam1,iParam2,iAlgo} = zeros(nInstance,1);
            
            best.Niter{iParam1,iParam2,iAlgo} = zeros(nInstance,2);
            best.runtime{iParam1,iParam2,iAlgo} = zeros(nInstance,2);
        end

        for iInstance = 1:nInstance
            
            fprintf('=============================================================================================================================================================================\n');
            fprintf('>>>>>>>>>> Instance %d <<<<<<<<<<\n', iInstance);
            
            %% Measuremt matrix
            % spatial mixer A
            if strcmpi(spatialMixerType,'no')
                A = @(XX,forward) XX;
                sigmaA_min = 1;     % smallest singular value of A
                sigmaA_max = 1;     % largest singular value of A
                sigmaA_minnz = 1;   % smallest nonzero singular value of A
            elseif contains(spatialMixerType,'cdp','IgnoreCase',true)
                mask = squeeze( create_cdp_masks(N,1,spatialSampleRate,spatialMixerType) );
                A = @(XX,forward) op_cdp(XX,mask,forward);
                A_mat = A(eye(N),1);  % spatial mixer A in matrix version
                sigmaA_min = svds(A_mat,1,'smallest');      % smallest singular value of A
                sigmaA_max = svds(A_mat,1);                 % largest singular value of A
                sigmaA_minnz = svds(A_mat,1,'smallestnz');  % smallest nonzero singular value of A
            else
                if strcmpi(spatialMixerType,'gauss') % spatial mixer A is gaussian random
                    A = (randn(M1,N) + 1i*randn(M1,N)) / sqrt(2);
                    %             A = A ./ sqrt(sum(abs(A).^2,1));    % normalize columns
                elseif strcmpi(spatialMixerType,'random_phase')
                    A = (randn(M1,N) + 1i*randn(M1,N));
                    A = sign(A);
                    A(abs(A) <= eps) = 1;
                else % spatial mixer A is dft
                    A = exp(-2i*pi*[0:M1-1].'*[0:N-1]/M1);  % First N columns of M1xM1 DFT matrix. Equivalent to extending length of RX signal to M1 by padding zeros at the end, when M1 > N
                end
                sigmaA_min = svds(A,1,'smallest');      % smallest singular value of A
                sigmaA_max = svds(A,1);                 % largest singular value of A
                sigmaA_minnz = svds(A,1,'smallestnz');  % smallest nonzero singular value of A
                A = @(XX,forward) op_general_mat(XX,A,forward,1);   % make function handle
            end
            
            % temporal mixer B
            if strcmpi(temporalMixerType,'no') % no temporal mixer, B is identity
                B = @(XX,forward) XX;
                sigmaB_min = 1;     % smallest singular value of B
                sigmaB_max = 1;     % largest singular value of B
                sigmaB_minnz = 1;   % smallest nonzero singular value of B
            elseif strcmpi(temporalMixerType,'gauss') % temporal mixer B is gaussian random
                B = (randn(I,M2) + 1i*randn(I,M2)) / sqrt(2);
                %             B = B ./ sqrt(sum(abs(B).^2,2));    % normalize rows
                sigmaB_min = svds(B,1,'smallest');      % smallest singular value of B
                sigmaB_max = svds(B,1);                 % largest singular value of B
                sigmaB_minnz = svds(B,1,'smallestnz');  % smallest nonzero singular value of B
                B = @(XX,forward) op_general_mat(XX,B,forward,0);   % make funciton handle
            elseif contains(temporalMixerType,'stcdp','IgnoreCase',true) % temporal mixer B is short-time cdp
                mask = squeeze( create_cdp_masks(len_stcdp,1,temporalSampleRate,temporalMixerType) );
                %                 mask = [ ones(8,1), [zeros(2,1);ones(6,1)], [zeros(4,1);ones(4,1)], [zeros(6,1);ones(2,1)] ];
                B = @(XX,forward) op_stcdp(XX,mask,forward);
                B_short = op_cdp(eye(len_stcdp),mask,1);
                sigmaB_min = svds(B_short,1,'smallest');
                sigmaB_max = svds(B_short,1);
                sigmaB_minnz = svds(B_short,1,'smallestnz');
            else % temporal mixer B is short-time fourier transform
                B = @(XX,forward) op_stft(XX,stft_window_len,stft_hopsize,stft_fft_len,forward);
                B_mat = B(eye(I),1);
                sigmaB = svd(B_mat);
                sigmaB_max = sigmaB(1);
                sigmaB_min = sigmaB(end);
                sigmaB_minnz = min( sigmaB( sigmaB > eps ) );
            end
            
            %% Generate random channel matrix
            if contains(channelType,'gauss','IgnoreCase',true)
                % 1. Spatial independent Gaussian channel
                D_true = (randn(N,P) + 1i*randn(N,P)) / sqrt(2);
                %             D_true = D_true ./ sqrt( sum(abs(D_true).^2,1) );   % normalize columns
            else
                % 2. Users with equally spaced DoAs, only LOS and equal gain
                direction = linspace(0,180,P+2);    % angle of users (in degree), equally placed
                direction = direction(2:end-1);
                phase = -pi*bsxfun(@times,[0:1:(N-1)]',cos(direction/180*pi));  % ULA with half-wavelength spacing
                D_true = exp(1i*phase);
            end
            % Corr_true = abs(bsxfun(@rdivide,bsxfun(@rdivide,D_true'*D_true,sqrt(sum(abs(D_true).^2,1))),sqrt(sum(abs(D_true).^2,1))')); % normalized correlation among columns of D_true
            
            %% Generate sparse TX signal
            Z_true = zeros(P,I);    % TX signal
            %             for ind = 1 : I/len_stcdp
            %                 supvec = randperm(P,L);% support vector
            %                 Z_true(supvec,((ind-1)*len_stcdp+1):(ind*len_stcdp)) = ( randn(L,len_stcdp) + 1i*randn(L,len_stcdp) ) / sqrt(2);
            %             end
            for ind = 1 : I
                supvec = randperm(P,L);% support vector
                Z_true(supvec,ind) = ( randn(L,1) + 1i*randn(L,1) ) / sqrt(2);
                %                 Z_true(supvec,ind) = ( 10 + randn(L,1) ) .* exp(1i*2*pi*rand(L,1));  % magnitude ~ N(10,1), uniformly distributed phase
                %                 Z_true(supvec,ind) = ( 1  ) .* exp(1i*2*pi*rand(L,1));
            end
            
            %             Z_true = ( rand(P,I) <= L/P );
            %             zeroColInd = true(1,I);
            %             while any(zeroColInd)
            %                 Z_true(:,zeroColInd) = ( rand(P,sum(zeroColInd)) <= L/P );
            % %                 Z_true = ( rand(P,I) <= L/P );
            %                 zeroColInd = ( sum(Z_true) == 0 );
            %             end
            
            patternZ_true = (abs(Z_true)>eps);   % sparsity pattern of TX signal
            Z_true(~patternZ_true) = 0; % remove non-zero values below eps
            %             Z_true_normalized = bsxfun( @rdivide, Z_true, sqrt( sum(abs(Z_true).^2,2)./sum(pattern_true,2) ) );
            
            %% RX signal
            X_true = D_true*Z_true;     % RX signal
            
            %% Magnitude measurements            
            if isinf(SNR)
                Y = abs( B( A(X_true,1), 1) );
            else
                Y = max(awgn(abs( B( A(X_true,1), 1) ),SNR,'measured','dB'),0);     % noisy magnitude-only measurements
            end
            
            % groundtruth values and measurement parameters
            data.D = D_true;
            data.Z = Z_true;
            data.X = X_true;
            data.normD = norm(D_true,'fro');
            data.normZ = norm(Z_true,'fro');
            data.normX = norm(X_true,'fro');
            data.A = A;
            data.B = B;
            data.L = L; % groundtruth sparsity level
            data.len_stcdp = len_stcdp;
            data.spatialMixerType = spatialMixerType;
            data.temporalMixerType = temporalMixerType;
            data.Y = Y;
            
            %% Try several random initializations
            % for recording solutions
            X_sol = cell(nAlgo,maxInitial);
            D_sol = cell(nAlgo,maxInitial);
            Z_sol = cell(nAlgo,maxInitial);
            patternZ_sol = cell(nAlgo,maxInitial);
            optVal = zeros(nAlgo,maxInitial);  % store the returned optimal value
            
            for iInitial = 1:maxInitial
                
                % Generate random initialization
                D0 = (randn(N,P) + 1i * randn(N,P)) / sqrt(2);
                D0 = D0 ./ sqrt(sum(abs(D0).^2,1));   % normalize each column
                X0 = (randn(N,I) + 1i*randn(N,I)) / sqrt(2);
                Z0 = D0 \ X0;     % sparsiest solution; = pinv(D0)*X0 in overdetermined case
                % Z0 = pinv(D0)*X0;
                % Z0 = (randn(P,I) + 1i * randn(P,I))/sqrt(2);
                
                % initialize as ground-truth
                % D0 = D_true ./ sqrt(sum(abs(D_true).^2,1));
                % Z0 = Z_true .* sqrt(sum(abs(D_true).^2,1)).';

                %% Perform estimation using different algorithms
                for iAlgo = 1:nAlgo
                    algoFunc = algoSet{iAlgo};
                    algoName = func2str(algoFunc);

                    params.blockRule = blockRule(iAlgo);

                    mu = mu_factor * sigmaA_minnz^2 * sigmaB_minnz^2;
                    params.mu = mu;

                    % compute upper bound of sparsity parameter
                    switch algoName
                        case 'FUN_PRDLsca'
                            % algorithm w/o variable X, Gaussian observation
                            % model
                            if strcmpi(temporalMixerType,'no')
                                lambda_max = sigmaA_max * sqrt(max(sum(abs(Y).^2)));
                            else
                                lambda_max = sigmaA_max * max( sum( abs(B(eye(I),1)) .* sqrt( sum(abs(Y).^2) ), 2 ) );
                            end
                        case {'FUN_PRDLscaX','scPRIME'}
                            % algorithm with variable X, Gaussian
                            % observation model
                            if strcmpi(temporalMixerType,'no')
                                if M1 >= N
                                    lambda_max = mu/(sigmaA_min^2 + mu) * sigmaA_max * sqrt(max(sum(abs(Y).^2)));
                                else
                                    lambda_max = sigmaA_max * sqrt(max(sum(abs(Y).^2)));
                                end
                            else
                                if M1 >= N && M2 >= I
                                    lambda_max = mu/(sigmaA_min^2 * sigmaB_min^2 + mu) * sigmaA_max * sigmaB_max * norm(Y,'fro');
                                else
                                    lambda_max = sigmaA_max * sigmaB_max * norm(Y,'fro');
                                end
                            end
                    end

                    lambda = lambda_factor * lambda_max;
                    params.lambda = lambda;
                    
%                     fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
%                     fprintf('>>>>>> Algorithm %d: %s <<<<<<\n',iAlgo,algoName);
%                     fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
%                     fprintf(' mu factor |    mu    | lambda factor |  lambda \n');
%                     fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
%                     fprintf(' %9f | %.2e | %13f | %.2e \n',...
%                         mu_factor, mu, lambda_factor, lambda);
%                     fprintf(' No. | mu factor |    mu    | lambda factor |  lambda  | Error X | Error DZ | Error D | Avg. Code Sparsity | RMSE Sparsity | FmeasureZ | Error Z | # iter | Runtime \n');
%                     fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    
                    if nargout(algoFunc) == 4
                        [D,Z,patternZ,stat] = algoFunc(data,params,D0,Z0);
                        X = D*Z;
                    else
                        [X,D,Z,patternZ,stat] = algoFunc(data,params,X0,D0,Z0);
                    end
                    
                    % record solution
                    X_sol{iAlgo,iInitial} = X;
                    D_sol{iAlgo,iInitial} = D;
                    Z_sol{iAlgo,iInitial} = Z;
                    patternZ_sol{iAlgo,iInitial} = patternZ;
                    optVal(iAlgo,iInitial) = stat.optVal;
                    
                    % evaluate and record estimation quality
                    clear sol;
                    sol.D = D;
                    sol.Z = Z;
                    sol.X = X;
                    quality = eval_quality(sol,data);
                    
                    precisionZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.precisionZ;
                    recallZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.recallZ;
                    FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.FmeasureZ;
                    errorX{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.errorX;
                    errorD{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.errorD;
                    errorZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.errorZ;
                    errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.errorAbsZ;
                    errorDZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.errorDZ;
                    sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.sparsityZ;
                    rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iInitial) = quality.rmse_sparsityZ;
                    
                    % record other statistics
                    Niter{iParam1,iParam2,iAlgo}(iInstance,iInitial) = stat.Niter;
                    runtime{iParam1,iParam2,iAlgo}(iInstance,iInitial) = stat.runtime(end);
                    
                    % print simulation info
%                     fprintf(' %3d | %9f | %.2e | %13f | %.2e | %6.5f | %7.6f | %6.5f | %18f | %13f | %9f | %6.5f | %6d | %5.1f s \n',...
%                         iInitial, mu_factor, mu, lambda_factor, lambda, errorX{iParam1,iParam2,iAlgo}(iInstance,iInitial), ...
%                         errorDZ{iParam1,iParam2,iAlgo}(iInstance,iInitial), errorD{iParam1,iParam2,iAlgo}(iInstance,iInitial), ...
%                         sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iInitial), rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iInitial), FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,iInitial),...
%                         errorZ{iParam1,iParam2,iAlgo}(iInstance,iInitial), Niter{iParam1,iParam2,iAlgo}(iInstance,iInitial), runtime{iParam1,iParam2,iAlgo}(iInstance,iInitial));

                end
            end
            
            %% Find the best solution from random initializations and perform debiasing step
            for ii = 1:initialSize
                nInitial = initialSet(ii);

                for iAlgo = 1:nAlgo
                    algoFunc = algoSet{iAlgo};
                    algoName = func2str(algoFunc);

                    % find the best one among the solutions from the frist
                    % nInitial random initializations according to
                    % objective function value
                    [~, optInitialInd] = min( optVal(iAlgo,1:nInitial) );

                    best.FmeasureZ{iParam1,iParam2,iAlgo}(iInstance) = FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);

                    best.precisionZ{iParam1,iParam2,iAlgo}(iInstance) = precisionZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.recallZ{iParam1,iParam2,iAlgo}(iInstance) = recallZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);

                    best.sparsityZ{iParam1,iParam2,iAlgo}(iInstance) = sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance) = rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.Niter{iParam1,iParam2,iAlgo}(iInstance,1) = Niter{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.runtime{iParam1,iParam2,iAlgo}(iInstance,1) = runtime{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);

                    best.errorX{iParam1,iParam2,iAlgo}(iInstance) = errorX{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.errorDZ{iParam1,iParam2,iAlgo}(iInstance) = errorDZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.errorD{iParam1,iParam2,iAlgo}(iInstance) = errorD{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.errorZ{iParam1,iParam2,iAlgo}(iInstance) = errorZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                    best.errorAbsZ{iParam1,iParam2,iAlgo}(iInstance) = errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);

                    optInitial{iParam1,iParam2,iAlgo}(iInstance) = optInitialInd;
                    
                    % perform debiasing if optimal lambda > 0 and support of Z is not empty
                    if debiasingOn && (lambda_factor > 0) && (best.sparsityZ{iParam1,iParam2,iAlgo}(iInstance) > 0) 
                        params.blockRule = blockRule(iAlgo);
                        mu = mu_factor * sigmaA_minnz^2 * sigmaB_minnz^2;
                        params.mu = mu;
                        params.lambda = -1; %params.verb = 1;
                        if nargout(algoFunc) == 4
                            [D,Z,~,stat] = algoFunc(data,params,D_sol{iAlgo,optInitialInd},Z_sol{iAlgo,optInitialInd});
                            X = D*Z;
                        else
                            [X,D,Z,~,stat] = algoFunc(data,params,X_sol{iAlgo,optInitialInd},D_sol{iAlgo,optInitialInd},Z_sol{iAlgo,optInitialInd});
                        end

                        clear sol;
                        sol.X = X;
                        sol.D = D;
                        sol.Z = Z;
                        quality = eval_quality(sol,data);

                        best.debiasedErrorX{iParam1,iParam2,iAlgo}(iInstance) = quality.errorX;
                        best.debiasedErrorDZ{iParam1,iParam2,iAlgo}(iInstance) = quality.errorDZ;
                        best.debiasedErrorD{iParam1,iParam2,iAlgo}(iInstance) = quality.errorD;
                        best.debiasedErrorZ{iParam1,iParam2,iAlgo}(iInstance) = quality.errorZ;
                        best.debiasedErrorAbsZ{iParam1,iParam2,iAlgo}(iInstance) = quality.errorAbsZ;

                        % record other statistics
                        best.Niter{iParam1,iParam2,iAlgo}(iInstance,2) = stat.Niter;
                        best.runtime{iParam1,iParam2,iAlgo}(iInstance,2) = stat.runtime(end);

                    else
                        best.debiasedErrorX{iParam1,iParam2,iAlgo}(iInstance) = errorX{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                        best.debiasedErrorDZ{iParam1,iParam2,iAlgo}(iInstance) = errorDZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                        best.debiasedErrorD{iParam1,iParam2,iAlgo}(iInstance) = errorD{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                        best.debiasedErrorZ{iParam1,iParam2,iAlgo}(iInstance) = errorZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);
                        best.debiasedErrorAbsZ{iParam1,iParam2,iAlgo}(iInstance) = errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd);

                        %                         precisionZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = precisionZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                        %                         recallZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = recallZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                        %                         FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);

                        %                         sparsityZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                        %                         rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);

                        %                 Niter{i_P,i_L}(i_instance,N_lambda+2) = 0;
                        %                 runtime{i_P,i_L}(i_instance,N_lambda+2) = 0;
                    end

                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf('Algorithm %d: %s: best among %d random initializations\n',iAlgo,algoName,initialSet(ii));
%                     fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
%                     fprintf(' No. initial | mu factor |    mu    | lambda factor |  lambda \n');
%                     fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
%                     fprintf(' %11d | %9f | %.2e | %13f | %.2e \n',...
%                         optInitialInd, mu_factor, mu, lambda_factor, lambda);
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(' No. initial | Avg. Code Sparsity | RMSE Sparsity | FmeasureZ | Error X | Error DZ | Error D | Error Z | Error AbsZ |    # iter   | Runtime \n');
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(' %11d | %18f | %13f | %9f | %6.5f | %7.6f | %6.5f | %6.5f | %10.5f | %4.0f + %4.0f | %3.1f + %3.1f s \n',...
                        optInitialInd,sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd), rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd),...
                        FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd), errorX{iParam1,iParam2,iAlgo}(iInstance,optInitialInd), ...
                        errorDZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd), errorD{iParam1,iParam2,iAlgo}(iInstance,optInitialInd), ...
                        errorZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd), errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,optInitialInd),...
                        best.Niter{iParam1,iParam2,iAlgo}(iInstance,1), best.Niter{iParam1,iParam2,iAlgo}(iInstance,2), ...
                        best.runtime{iParam1,iParam2,iAlgo}(iInstance,1), best.runtime{iParam1,iParam2,iAlgo}(iInstance,2));
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(' Debiased Errors:\n Error X | Error DZ | Error D | Error Z | Error AbsZ \n');
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(' %6.5f | %7.6f | %6.5f | %6.5f | %10.5f \n',...
                        best.debiasedErrorX{iParam1,iParam2,iAlgo}(iInstance), best.debiasedErrorDZ{iParam1,iParam2,iAlgo}(iInstance), ...
                        best.debiasedErrorD{iParam1,iParam2,iAlgo}(iInstance),...
                        best.debiasedErrorZ{iParam1,iParam2,iAlgo}(iInstance), best.debiasedErrorAbsZ{iParam1,iParam2,iAlgo}(iInstance));

                end


            end
        end
        
        %% Average results over random instances
        fprintf('=============================================================================================================================================================================\n');
        fprintf('>>>>>>>>>> Average Results over %d Instances <<<<<<<<<<\n',nInstance);
        fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
        
        for iAlgo = 1:nAlgo
            avgResults.pattern.sparsityZ(iParam1,iParam2,iAlgo) = mean(best.sparsityZ{iParam1,iParam2,iAlgo});
            avgResults.pattern.rmse_sparsityZ(iParam1,iParam2,iAlgo) = mean(best.rmse_sparsityZ{iParam1,iParam2,iAlgo});
            
            avgResults.pattern.precisionZ(iParam1,iParam2,iAlgo) = mean(best.precisionZ{iParam1,iParam2,iAlgo});
            avgResults.pattern.recallZ(iParam1,iParam2,iAlgo) = mean(best.recallZ{iParam1,iParam2,iAlgo});
            avgResults.pattern.FmeasureZ(iParam1,iParam2,iAlgo) = mean(best.FmeasureZ{iParam1,iParam2,iAlgo});
            
            avgResults.error.D(iParam1,iParam2,iAlgo) = mean(best.errorD{iParam1,iParam2,iAlgo});
            avgResults.error.Z(iParam1,iParam2,iAlgo) = mean(best.errorZ{iParam1,iParam2,iAlgo});
            avgResults.error.AbsZ(iParam1,iParam2,iAlgo) = mean(best.errorAbsZ{iParam1,iParam2,iAlgo});
            avgResults.error.X(iParam1,iParam2,iAlgo) = mean(best.errorX{iParam1,iParam2,iAlgo});
            avgResults.error.DZ(iParam1,iParam2,iAlgo) = mean(best.errorDZ{iParam1,iParam2,iAlgo});
            
            avgResults.debiasedError.D(iParam1,iParam2,iAlgo) = mean(best.debiasedErrorD{iParam1,iParam2,iAlgo});
            avgResults.debiasedError.Z(iParam1,iParam2,iAlgo) = mean(best.debiasedErrorZ{iParam1,iParam2,iAlgo});
            avgResults.debiasedError.AbsZ(iParam1,iParam2,iAlgo) = mean(best.debiasedErrorAbsZ{iParam1,iParam2,iAlgo});
            avgResults.debiasedError.X(iParam1,iParam2,iAlgo) = mean(best.debiasedErrorX{iParam1,iParam2,iAlgo});
            avgResults.debiasedError.DZ(iParam1,iParam2,iAlgo) = mean(best.debiasedErrorDZ{iParam1,iParam2,iAlgo});
            
            avgResults.complexity.Niter(iParam1,iParam2,iAlgo,1) = mean( sum(Niter{iParam1,iParam2,iAlgo},2), 1 );
            avgResults.complexity.Niter(iParam1,iParam2,iAlgo,2) = mean(best.Niter{iParam1,iParam2,iAlgo}(:,2),1);
            avgResults.complexity.runtime(iParam1,iParam2,iAlgo,1) = mean( sum(runtime{iParam1,iParam2,iAlgo},2), 1);
            avgResults.complexity.runtime(iParam1,iParam2,iAlgo,2) = mean(best.runtime{iParam1,iParam2,iAlgo}(:,2),1);
            
            
            fprintf('>>>>>> Algorithm %d: %s <<<<<<\n',iAlgo,func2str(algoSet{iAlgo}));
            fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(' Avg. Code Sparsity | RMSE Sparsity | FmeasureZ | Error X | Error DZ | Error D | Error Z | Error AbsZ |    # iter   | Runtime \n');
            fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(' %18f | %13f | %9f | %6.5f | %7.6f | %6.5f | %6.5f | %10.5f | %4.0f + %4.0f | %3.1f + %3.1f s \n',...
                avgResults.pattern.sparsityZ(iParam1,iParam2,iAlgo), avgResults.pattern.rmse_sparsityZ(iParam1,iParam2,iAlgo), avgResults.pattern.FmeasureZ(iParam1,iParam2,iAlgo), ...
                avgResults.error.X(iParam1,iParam2,iAlgo), avgResults.error.DZ(iParam1,iParam2,iAlgo), avgResults.error.D(iParam1,iParam2,iAlgo),...
                avgResults.error.Z(iParam1,iParam2,iAlgo), avgResults.error.AbsZ(iParam1,iParam2,iAlgo), ...
                avgResults.complexity.Niter(iParam1,iParam2,iAlgo,1), avgResults.complexity.Niter(iParam1,iParam2,iAlgo,2), avgResults.complexity.runtime(iParam1,iParam2,iAlgo,1), avgResults.complexity.runtime(iParam1,iParam2,iAlgo,2));
            fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(' Debiased Errors:\n Error X | Error DZ | Error D | Error Z | Error AbsZ \n');
            fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(' %6.5f | %7.6f | %6.5f | %6.5f | %10.5f \n',...
                avgResults.debiasedError.X(iParam1,iParam2,iAlgo), avgResults.debiasedError.DZ(iParam1,iParam2,iAlgo),...
                avgResults.debiasedError.D(iParam1,iParam2,iAlgo), avgResults.debiasedError.Z(iParam1,iParam2,iAlgo), avgResults.debiasedError.AbsZ(iParam1,iParam2,iAlgo) );
            fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
        end
        
        avgResults.pattern.FmeasureZ_ref(iParam1,iParam2) = harmmean([1,L/P]);  % Fmeasure of a full matrix as a reference

        
    end
    
    %     save(sprintf('N=%d_I=%d_M1=%d_M2=%d_SNR=%d_P=%s_L=%s_GaussianChannel_measurementType_spatial=%s_temporal=%s_MCruns=%d_%s.mat',...
    %         N,I,M1,M2,SNR,mat2str(P_set(1:i_P)),mat2str(L_set),spatialMixerType,temporalMixerType,N_instance,datestr(datetime,30)));
    
end

diary off;

%% Save data
% save([datestr(datetime,30),'.mat']);
fileName = sprintf('%s=%s_%s=%s_Ninitial=%s_N=%d_channel=%s_measurementType_spatial=%s_temporal=%s_MCruns=%d_%s.mat',...
    paramName1,mat2str(paramSet1),paramName2,mat2str(paramSet2),mat2str(initialSet),N,channelType,spatialMixerType,temporalMixerType,nInstance,datestr(datetime,30));
save(fileName);

%% Plot results
% if strcmpi(paramName2,'M1')
%     paramName2 = 'oversampling rate';
%     paramSet2 = paramSet2 / N;
% end
% 
% % create string for legend
% legend_str = cell(N_param1*N_param2,1);
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         legend_str{(i_param1-1)*N_param2 + i_param2} = sprintf('%s=%d, %s=%d',paramName1,paramSet1(i_param1),paramName2,paramSet2(i_param2));
%     end
% end
% 
% %% Complexity
% figure;
% % plot results of algorithm 1
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, squeeze( sum(avgResults.complexity.Niter(i_param1,i_param2,1,:,:),4) ));
%         hold on;
%     end
% end
% % plot results of algorithm 2
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, squeeze( sum(avgResults.complexity.Niter(i_param1,i_param2,2,:,:),4) ),'--');
%         hold on;
%     end
% end
% hold off;
% grid on;
% xlabel('number of initializations');
% ylabel('number of iterations');
% legend(legend_str);
% matlab2tikz('Niter_varyInitial.tikz');
% 
% figure;
% % plot results of algorithm 1
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, squeeze( sum(avgResults.complexity.runtime(i_param1,i_param2,1,:,:),4) ));
%         hold on;
%     end
% end
% % plot results of algorithm 2
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, squeeze( sum(avgResults.complexity.runtime(i_param1,i_param2,2,:,:),4) ),'--');
%         hold on;
%     end
% end
% hold off;
% grid on;
% xlabel('number of initializations');
% ylabel('CPU time (seconds)');
% legend(legend_str);
% matlab2tikz('cputime_varyInitial.tikz');
% 
% % figure;
% % plot(sum(avgResults.complexity.Niter(:,:,1,:),4).');
% % hold on;
% % plot(sum(avgResults.complexity.Niter(:,:,2,:),4).','--');
% % hold off;
% % grid on;
% % xlabel(paramName2);
% % title('Number of iterations');
% % 
% % figure;
% % plot(sum(avgResults.complexity.runtime(:,:,1,:),4).');
% % hold on;
% % plot(sum(avgResults.complexity.runtime(:,:,2,:),4).','--');
% % hold off;
% % grid on;
% % xlabel(paramName2);
% % title('CPU time');
% 
% 
% %% Fmeasure
% figure;
% % plot results of algorithm 1
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, squeeze( avgResults.pattern.FmeasureZ(i_param1,i_param2,1,:) ));
%         hold on;
%     end
% end
% % plot results of algorithm 2
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, squeeze( avgResults.pattern.FmeasureZ(i_param1,i_param2,2,:) ),'--');
%         hold on;
%     end
% end
% % plot reference Fmeasure of Z
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, squeeze( avgResults.pattern.FmeasureZ_ref(i_param1,i_param2)*ones(initial_size,1)),'-.');
%         hold on;
%     end
% end
% hold off;
% grid on;
% xlabel('number of initializations');
% ylabel('F-measure(Z)');
% legend(legend_str);
% matlab2tikz('FmeasureZ_varyInitial.tikz');
% 
% %% Debiased Error
% figure;
% % plot results of algorithm 1
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, 20*log10(squeeze( avgResults.debiasedError.X(i_param1,i_param2,1,:) )) );
%         hold on;
%     end
% end
% % plot results of DZ of algorithm 2
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, 20*log10(squeeze( avgResults.debiasedError.DZ(i_param1,i_param2,:) )),'--');
%         hold on;
%     end
% end
% % plot results of X of algorithm 2
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, 20*log10(squeeze( avgResults.debiasedError.X(i_param1,i_param2,2,:))),'-.');
%         hold on;
%     end
% end
% hold off;
% grid on;
% xlabel('number of initializations');
% ylabel('MNSE(X) (dB)');
% legend(legend_str);
% matlab2tikz('errorX_varyInitial.tikz');
% 
% figure;
% % plot results of algorithm 1
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, 20*log10(squeeze( avgResults.debiasedError.D(i_param1,i_param2,1,:) )) );
%         hold on;
%     end
% end
% % plot results of algorithm 2
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, 20*log10(squeeze( avgResults.debiasedError.D(i_param1,i_param2,2,:) )),'--');
%         hold on;
%     end
% end
% hold off;
% grid on;
% xlabel('number of initializations');
% ylabel('MNSE(D) (dB)');
% legend(legend_str);
% matlab2tikz('errorD_varyInitial.tikz');
% 
% figure;
% % plot results of algorithm 1
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, 20*log10(squeeze( avgResults.debiasedError.Z(i_param1,i_param2,1,:) )));
%         hold on;
%     end
% end
% % plot results of algorithm 2
% for i_param1 = 1:N_param1
%     for i_param2 = 1:N_param2
%         plot(initial_set, 20*log10(squeeze( avgResults.debiasedError.Z(i_param1,i_param2,2,:) )),'--');
%         hold on;
%     end
% end
% hold off;
% grid on;
% xlabel('number of initializations');
% ylabel('MNSE(Z) (dB)');
% legend(legend_str);
% matlab2tikz('errorZ_varyInitial.tikz');
