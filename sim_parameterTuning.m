% Simulation with different values for the regularization parameters


addpath('tools/');
% clc;
clear all;
% close all;

diary off;
diary(['diary_',datestr(date,'yyyy-mm-dd'),'.txt']);

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
    'gauss';    % Gaussian
    'pois';     % Poisson
    };
measurementModelNumber = 1;
measurementModel = measurementModelList{measurementModelNumber};

N = 64; %100;       % no. of Rx antennas
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
photonScale = 0.1;


% regularization parameters for different algorithms

% % muFactor = 10^1.5;     % for sampling rate 2 x 2
% muFactor = 10;

% lambdaSet = [0.01:0.01:0.2, 0.25:0.05:0.5];
% lambdaSet = [logspace(-3,-2,11)];
% lambdaSet = lambdaSet(1:end-1);
% lambdaSet = [lambdaSet, logspace(-2,-1,3)];
% lambdaSet = logspace(-3,0,16);
% lambdaSet = logspace(-2,-1.8,11);
% lambdaSet = 0.00631;
% lambdaSet = [0.25:0.005:0.3];
% lambdaSet = [0.1];
% lambdaSet = 10^-1.9; % 0.01;
% lambdaSet = [0, 0.75.^[20:-1:11] ];
% lambdaSet = 0.75.^[11];%,10,9];
% lambdaSet = [0, 0.75.^[27:-1:18] ];
% lambdaSet = [0.75.^[27:-1:18] ];
% lambdaSet = [0, 0.75.^[26:-1:20] ];    % for SNR = 15, 20
% lambdaSet = 0.75^22;
lambdaSet = [0.75.^[18:-1:12] ];    % gaussian spatial mixer, no temporal mixer, sampling rate 4 x 1, SNR = 15
% lambdaSet = [ 0.75.^[21:-1:16] ];      % 4 x 4 cdp_complex, SNR = 15
% lambdaSet = 0.75^19;
% lambdaSet = [0, 0.75.^[30:-1:11] ];    % 4 x 4 gauss
% lambdaSet = 0.75^14;
% lambdaSet = [0, 0.75.^[20:-1:19] ];
nLambda = length(lambdaSet);

muSet = 1; % logspace(0,2,11);

% rhoSet = [0, 0.75.^[20:-1:11] ];
% rhoSet = [0, 0.75.^[32:-1:23] ];   % 4 x 4 sampling rate
% rhoSet = [0, 0.75.^[30:-1:25] ];   % for SNR = 15, 20
% rhoSet = 0.75^28;
% rhoSet = logspace(-6,-3,16);
% rhoSet = [0, 0.75.^[24:-1:15] ];
rhoSet = [0.75.^[18:-1:12] ];   % gaussian spatial mixer, no temporal mixer, sampling rate 4 x 1, SNR = 15
% rhoSet = [0, 0.75.^[16:-1:10] ];   % cdp spatial mixer, no temporal mixer, sampling rate 4 x 1, SNR = 15
% rhoSet = [ 0.75.^[28:-1:22] ];     % 4 x 4 cdp_complex, SNR = 15
% rhoSet = 0.75^19;
% rhoSet = [0, 0.75.^[35:-1:16] ];   % 4 x 4 gauss
% rhoSet = [0, logspace(-3,3,26)];
% rhoSet = 0.11;
% rhoSet = [0, 0.75.^[16:-1:10] ];
% rhoSet = 0.75^13;
nRho = length(rhoSet);

% set of testing algorithms
algoSet = {
    @FUN_PRDLsca;
    @FUN_PRDLscaX;
    @scPRIME;
    };
nAlgo = length(algoSet);    % no. of algorithms to test
blockRule = zeros(nAlgo,1);  % indicates block update rule, represents no. of blocks NOT to update in each iteration
debiasingOn = 0;

approxParamSet = {
    1;
    muSet;
    muSet;
    };
if length(approxParamSet) < nAlgo
    fprintf('Inconsistent simulation parameters. Aborting.\n');
    return;
end
sparsityParamSet = {
    lambdaSet;
    rhoSet;
    rhoSet;
    };
if length(sparsityParamSet) < nAlgo
    fprintf('Inconsistent simulation parameters. Aborting.\n');
    return;
end
nApproxParam = zeros(nAlgo,1);
nSparsityParam = zeros(nAlgo,1);
for iAlgo = 1:nAlgo
    nApproxParam(iAlgo) = length(approxParamSet{iAlgo});
    nSparsityParam(iAlgo) = length(sparsityParamSet{iAlgo});
end

nInstance = 1;     % no. of random instances
% N_initial = 1;  % # of random initializations to test

% set parameters for algorithms
params.maxIter = 2000;
params.tol = 1e-5;
params.verb = 0;   % show information in each iteration if verb = 1
params.N = N;
params.I = I;
params.P = P;
params.L = L;
params.trackConvergence = 0;

% varying parameters
paramList = {
    'P';    % 1: no. of users P
    'L';    % 2: sparsity level L
    'M1';   % 3: no. of spatial measurements M1
    'SNR';  % 4: SNR in dB for Gaussian observation model
    'I';    % 5: no. of time-slots I
    };
paramSetList = {
    [N]; % [N/2 N N*2];        % test set of P
    round([0.05, 0.1]*N); % [1:N/2];            % test set of L
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
fprintf('Varying parameters: %s = %s, %s = %s\n',paramName1,mat2str(paramSet1),paramName2,mat2str(paramSet2));
fprintf('P = %d, L = %d, N = %d, I = %d, M1 = %d, M2 = %d\n',P,L,N,I,M1,M2);
fprintf('Tolerance = %1.0e, max no. iterations = %d, Monte-Carlo runs = %d\n',params.tol,params.maxIter,nInstance);
fprintf('Channel type: %s\n', channelType);
fprintf('Measurement type: spatial mixer = %s, temporal mixer = %s\n', spatialMixerType, temporalMixerType);
fprintf('spatial oversampling rate = %.2f, temporal oversampling rate = %.2f\n',spatialSampleRate,temporalSampleRate);
if contains(temporalMixerType,'stft','IgnoreCase',true)
    fprintf('Window length = %d, hop size = %d, no. FFT samples = %d\n\n',stft_window_len,stft_hopsize,stft_fft_len);
end
if contains(measurementModel,'gauss','IgnoreCase',true)
    fprintf('Observation model: Gaussian, SNR = %d dB\n',SNR);
else
    fprintf('Observation model: Poisson, photon scale = %f\n',photonScale);
end
% for iAlgo = 1:nAlgo
%     fprintf('Algorithm %d: %s\nmu factor = %s\nrho factor = [0, 0.75^%s]\n',iAlgo,func2str(algoSet{iAlgo}),mat2str(approxParamSet{iAlgo},3),mat2str(log(sparsityParamSet{iAlgo}(2:end))/log(0.75),3));
% end
for iAlgo = 1:nAlgo
    fprintf('Algorithm %d: %s\nmu factor = %s\nlambda factor = [0.75^%s]\n',iAlgo,func2str(algoSet{iAlgo}),mat2str(approxParamSet{iAlgo},3),mat2str(log(sparsityParamSet{iAlgo})/log(0.75),3));
end

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

best.Params = cell(nParam1,nParam2,nAlgo);

experimentalSNR = cell(nParam1,nParam2);

% objVal = cell(nAlgo,1);
% runtimeArray = cell(nAlgo,1);

%% Average results
avgResults.pattern.sparsityZ = zeros(nParam1,nParam2,nAlgo);
avgResults.pattern.rmse_sparsityZ = zeros(nParam1,nParam2,nAlgo);

avgResults.pattern.precisionZ = zeros(nParam1,nParam2,nAlgo);
avgResults.pattern.recallZ = zeros(nParam1,nParam2,nAlgo);
avgResults.pattern.FmeasureZ = zeros(nParam1,nParam2,nAlgo);
avgResults.pattern.FmeasureZ_ref = zeros(nParam1,nParam2);

avgResults.error.D = zeros(nParam1,nParam2,nAlgo);
avgResults.error.Z = zeros(nParam1,nParam2,nAlgo);
avgResults.error.AbsZ = zeros(nParam1,nParam2,nAlgo);
avgResults.error.X = zeros(nParam1,nParam2,nAlgo);
avgResults.error.DZ = zeros(nParam1,nParam2,nAlgo);

avgResults.debiasedError.D = zeros(nParam1,nParam2,nAlgo);
avgResults.debiasedError.Z = zeros(nParam1,nParam2,nAlgo);
avgResults.debiasedError.AbsZ = zeros(nParam1,nParam2,nAlgo);
avgResults.debiasedError.X = zeros(nParam1,nParam2,nAlgo);
avgResults.debiasedError.DZ = zeros(nParam1,nParam2,nAlgo);

avgResults.complexity.Niter = zeros(nParam1,nParam2,nAlgo,2);
avgResults.complexity.runtime = zeros(nParam1,nParam2,nAlgo,2);

avgResults.experimentalSNR = zeros(nParam1,nParam2);

avgResults.objVal = zeros(nAlgo,params.maxIter+1);
avgResults.runtimeArray = zeros(nAlgo,params.maxIter);

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
            nMu = nApproxParam(iAlgo);
            nLambda = nSparsityParam(iAlgo);
            
            errorD{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            
            sparsityZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            rmse_sparsityZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            precisionZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            recallZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            FmeasureZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            errorZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            errorAbsZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            
            errorX{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            errorDZ{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            
            Niter{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            runtime{iParam1,iParam2,iAlgo} = zeros(nInstance,nMu,nLambda);
            
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
            
            best.Params{iParam1,iParam2,iAlgo} = zeros(nInstance,2);
            
%             objVal{iAlgo} = zeros(nInstance,params.maxIter+1);
%             runtimeArray{iAlgo} = zeros(nInstance,params.maxIter+1);
            
        end
        experimentalSNR{iParam1,iParam2} = zeros(nInstance,1);
        
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
            
            %% Generate intensity measurements
            Yhat_true = B( A(X_true,1), 1);
            if contains(measurementModel,'gauss','IgnoreCase',true)
                % Gaussian observation model
                if isinf(SNR)
                    Y = abs(Yhat_true);
                else
                    Y = max(awgn(abs(Yhat_true),SNR,'measured','dB'),0);     % noisy magnitude-only measurements
                end
                Y = Y.^2;   % convert to intensity measurements
            else
                % Poisson observation model
                Y = poissrnd( photonScale * abs(Yhat_true).^2 ) / photonScale;
                experimentalSNR{iParam1,iParam2}(iInstance) = photonScale * sum(abs(Yhat_true).^4,'all') / norm(Yhat_true,'fro')^2;
                
            end
            nnzY = (Y>0);
            
            % groundtruth values
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
            
            %% Generate random initialization
            X0 = (randn(N,I) + 1i*randn(N,I)) / sqrt(2);
            Yhat0 = B( A(X0,1), 1);
            if contains(measurementModel,'pois','IgnoreCase',true)
                while any(abs(Yhat0(nnzY))<eps)
                    X0 = (randn(N,I) + 1i*randn(N,I)) / sqrt(2);
                    Yhat0 = B( A(X0,1), 1);
                end
            end
            D0 = (randn(N,P) + 1i * randn(N,P)) / sqrt(2);
            D0 = D0 ./ sqrt(sum(abs(D0).^2,1));   % normalize each column
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
                nMu = nApproxParam(iAlgo);
                muSet = approxParamSet{iAlgo};
                nLambda = nSparsityParam(iAlgo);
                lambdaSet = sparsityParamSet{iAlgo};
                
                params.blockRule = blockRule(iAlgo);
                
                D_sol = cell(nMu,nLambda);
                Z_sol = cell(nMu,nLambda);
                X_sol = cell(nMu,nLambda);
                patternZ_sol = cell(1,nLambda);
                
                % >>>>>> Search for lambda
                
                fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                fprintf('>>>>>> Algorithm %d: %s <<<<<<\n',iAlgo,algoName);

                % grid search for regularization parameters
                for iMu = 1:nMu
                    mu = muSet(iMu) * sigmaA_minnz^2 * sigmaB_minnz^2;
                    params.mu = mu;
                    
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(' No. mu | mu factor |    mu    \n');
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(' %6d | %9f | %.2e\n',iMu,muSet(iMu),mu);
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(' No. lambda | lambda factor |  lambda  | Objective | Error X | Error DZ | Error D | Avg. Code Sparsity | RMSE Sparsity | FmeasureZ | Error Z | Error AbsZ | # iter | Runtime \n');
                    fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    
                    % compute upper bound of sparsity parameter
                    switch algoName
                        case 'FUN_PRDLsca'
                            % algorithm w/o variable X, Gaussian observation
                            % model
                            if strcmpi(temporalMixerType,'no')
                                lambda_max = sigmaA_max * sqrt(max(sum(Y,1))); % sigmaA_max * sqrt(max(sum(abs(Y).^2)));
                            else
                                lambda_max = sigmaA_max * max( sum( abs(B(eye(I),1)) .* sqrt( sum(Y,1) ), 2 ) ); % sigmaA_max * max( sum( abs(B(eye(I),1)) .* sqrt( sum(abs(Y).^2) ), 2 ) );
                            end
                        case {'FUN_PRDLscaX','scPRIME'}
                            % algorithm with variable X, Gaussian
                            % observation model
                            if strcmpi(temporalMixerType,'no')
                                if M1 >= N
                                    lambda_max = mu/(sigmaA_min^2 + mu) * sigmaA_max * sqrt(max(sum(Y,1))); % mu/(sigmaA_min^2 + mu) * sigmaA_max * sqrt(max(sum(abs(Y).^2)));
                                else
                                    lambda_max = sigmaA_max * sqrt(max(sum(Y,1))); % sigmaA_max * sqrt(max(sum(abs(Y).^2)));
                                end
                            else
                                if M1 >= N && M2 >= I
                                    lambda_max = mu/(sigmaA_min^2 * sigmaB_min^2 + mu) * sigmaA_max * sigmaB_max * sqrt(sum(Y,'all')); % mu/(sigmaA_min^2 * sigmaB_min^2 + mu) * sigmaA_max * sigmaB_max * norm(Y,'fro');
                                else
                                    lambda_max = sigmaA_max * sigmaB_max * sqrt(sum(Y,'all')); % sigmaA_max * sigmaB_max * norm(Y,'fro');
                                end
                            end
                    end
                    
                    for iLambda = 1:nLambda
                        lambda = lambdaSet(iLambda) * lambda_max;
                        
                        if L == P && lambda > 0
                            continue;
                        end
                        
                        params.lambda = lambda;
                        if nargout(algoFunc) == 4
                            data.Y = sqrt(Y);
                            [D,Z,patternZ,stat] = algoFunc(data,params,D0,Z0);
                            X = D*Z;
                        else
                            data.Y = sqrt(Y);
                            [X,D,Z,patternZ,stat] = algoFunc(data,params,X0,D0,Z0);
                        end
                        
                        % record solution
                        X_sol{iMu,iLambda} = X;
                        D_sol{iMu,iLambda} = D;
                        Z_sol{iMu,iLambda} = Z;
                        patternZ_sol{iMu,iLambda} = patternZ;
                        
                        % evaluate and record estimation quality
                        clear sol;
                        sol.X = X;
                        sol.D = D;
                        sol.Z = Z;
                        quality = eval_quality(sol,data);
                        
                        precisionZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.precisionZ;
                        recallZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.recallZ;
                        FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.FmeasureZ;
                        errorX{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.errorX;
                        errorDZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.errorDZ;
                        errorD{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.errorD;
                        errorZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.errorZ;
                        errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.errorAbsZ;
                        sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.sparsityZ;
                        rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = quality.rmse_sparsityZ;
                        
                        % record other statistics
                        Niter{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = stat.Niter;
                        runtime{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda) = stat.runtime(end);
                        
%                         objVal{iAlgo}(iInstance,:) = stat.objVal;
%                         runtimeArray{iAlgo}(iInstance,:) = stat.runtime;
                        
                        fprintf(' %10d | %13f | %.2e | %.3e | %6.5f | %7.6f | %6.5f | %18f | %13f | %9f | %6.5f | %10.5f | %6d | %5.1f s \n',...
                            iLambda, lambdaSet(iLambda), lambda, stat.optVal, errorX{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda), errorDZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda), errorD{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda),...
                            sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda), rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda), FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda),...
                            errorZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda),errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda), Niter{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda), runtime{iParam1,iParam2,iAlgo}(iInstance,iMu,iLambda));
                    end
                end
                
                % >>>>>> Find the best regularization parameters and
                % perform debiasing for this instance
                % Find the optimal lambda value according to F-measure of
                % estimated sparsity pattern of Z
                maxVal = max( FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,:,:),[],'all');
                best.FmeasureZ{iParam1,iParam2,iAlgo}(iInstance) = maxVal;
                [bestMuInd,bestLambdaInd] = find( permute(FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,:,:),[2 3 1]) == maxVal, 1,'last');
                
                best.precisionZ{iParam1,iParam2,iAlgo}(iInstance) = precisionZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.recallZ{iParam1,iParam2,iAlgo}(iInstance) = recallZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                
                best.sparsityZ{iParam1,iParam2,iAlgo}(iInstance) = sparsityZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance) = rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.Niter{iParam1,iParam2,iAlgo}(iInstance,1) = Niter{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.runtime{iParam1,iParam2,iAlgo}(iInstance,1) = runtime{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                
                best.Params{iParam1,iParam2,iAlgo}(iInstance,:) = [muSet(bestMuInd), lambdaSet(bestLambdaInd)];
                best.ParamsInd{iParam1,iParam2,iAlgo}(iInstance,:) = [bestMuInd,bestLambdaInd];
                
                best.errorX{iParam1,iParam2,iAlgo}(iInstance) = errorX{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.errorDZ{iParam1,iParam2,iAlgo}(iInstance) = errorDZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.errorD{iParam1,iParam2,iAlgo}(iInstance) = errorD{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.errorZ{iParam1,iParam2,iAlgo}(iInstance) = errorZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                best.errorAbsZ{iParam1,iParam2,iAlgo}(iInstance) = errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                
                if debiasingOn && (lambdaSet(bestLambdaInd) > 0) && (best.sparsityZ{iParam1,iParam2,iAlgo}(iInstance) > 0) % perform debiasing if optimal lambda > 0 and support of Z is not empty
                    % debiasing step
                    params.lambda = -1; %params.verb = 1;
                    if nargout(algoFunc) == 4
                        data.Y = sqrt(Y);
                        [D,Z,~,stat] = algoFunc(data,params,D_sol{bestMuInd,bestLambdaInd},Z_sol{bestMuInd,bestLambdaInd});
                        X = D*Z;
                    else
                        data.Y = sqrt(Y);
                        [X,D,Z,~,stat] = algoFunc(data,params,X_sol{bestMuInd,bestLambdaInd},D_sol{bestMuInd,bestLambdaInd},Z_sol{bestMuInd,bestLambdaInd});
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
                    best.debiasedErrorX{iParam1,iParam2,iAlgo}(iInstance) = errorX{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                    best.debiasedErrorDZ{iParam1,iParam2,iAlgo}(iInstance) = errorDZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                    best.debiasedErrorD{iParam1,iParam2,iAlgo}(iInstance) = errorD{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                    best.debiasedErrorZ{iParam1,iParam2,iAlgo}(iInstance) = errorZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                    best.debiasedErrorAbsZ{iParam1,iParam2,iAlgo}(iInstance) = errorAbsZ{iParam1,iParam2,iAlgo}(iInstance,bestMuInd,bestLambdaInd);
                    
                    %                         precisionZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = precisionZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                    %                         recallZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = recallZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                    %                         FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = FmeasureZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                    
                    %                         sparsityZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                    %                         rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,nLambda+1) = rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance,optLambdaInd);
                    
                    %                 Niter{i_P,i_L}(i_instance,N_lambda+2) = 0;
                    %                 runtime{i_P,i_L}(i_instance,N_lambda+2) = 0;
                end
                
                fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                fprintf('Optimal regularization parameters:\n');
                fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                fprintf(' No. mu | mu factor |    mu    | No. lambda | lambda factor | lambda \n')
                fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                fprintf(' %6d | %9f | %.2e | %10d | %13f | %.2e\n',bestMuInd,muSet(bestMuInd),best.Params{iParam1,iParam2,iAlgo}(iInstance,1), ...
                    bestLambdaInd,lambdaSet(bestLambdaInd),best.Params{iParam1,iParam2,iAlgo}(iInstance,2));
                fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                fprintf(' Debiased Error X | Debiased Error DZ | Debiased Error D | Avg. Code Sparsity | RMSE Sparsity | FmeasureZ | Debiased Error Z | Debiased Error AbsZ |    # iter   | Runtime \n');
                fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                fprintf(' %15.5f | %16.6f | %15.5f | %18f | %13f | %9f | %15.5f | %19.5f | %4d + %4d | %3.1f + %3.1f s \n',...
                    best.debiasedErrorX{iParam1,iParam2,iAlgo}(iInstance), best.debiasedErrorDZ{iParam1,iParam2,iAlgo}(iInstance), ...
                    best.debiasedErrorD{iParam1,iParam2,iAlgo}(iInstance), best.sparsityZ{iParam1,iParam2,iAlgo}(iInstance), ...
                    best.rmse_sparsityZ{iParam1,iParam2,iAlgo}(iInstance), best.FmeasureZ{iParam1,iParam2,iAlgo}(iInstance), ...
                    best.debiasedErrorZ{iParam1,iParam2,iAlgo}(iInstance), best.debiasedErrorAbsZ{iParam1,iParam2,iAlgo}(iInstance), ...
                    best.Niter{iParam1,iParam2,iAlgo}(iInstance,1), best.Niter{iParam1,iParam2,iAlgo}(iInstance,2), ...
                    best.runtime{iParam1,iParam2,iAlgo}(iInstance,1), best.runtime{iParam1,iParam2,iAlgo}(iInstance,2));
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
            
            avgResults.complexity.Niter(iParam1,iParam2,iAlgo,:) = mean(best.Niter{iParam1,iParam2,iAlgo},1);
            avgResults.complexity.runtime(iParam1,iParam2,iAlgo,:) = mean(best.runtime{iParam1,iParam2,iAlgo},1);
            
%             avgResults.objVal(iAlgo,:) = mean(objVal{iAlgo},1);
%             avgResults.runtimeArray(iAlgo,:) = mean(runtimeArray{iAlgo},1);
            
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
        avgResults.experimentalSNR(iParam1,iParam2) = mean(experimentalSNR{iParam1,iParam2});
        
        %         save(sprintf('N=%d_I=%d_M1=%d_M2=%d_SNR=%d_P=%s_L=%s_GaussianChannel_GaussianMeasurement_%s.mat',N,I,M1,M2,SNR,mat2str(P_set(1:i_P)),mat2str(L_set(1:i_L)),datestr(datetime,30)));
        
    end
    
    %     save(sprintf('N=%d_I=%d_M1=%d_M2=%d_SNR=%d_P=%s_L=%s_GaussianChannel_measurementType_spatial=%s_temporal=%s_MCruns=%d_%s.mat',...
    %         N,I,M1,M2,SNR,mat2str(P_set(1:i_P)),mat2str(L_set),spatialMixerType,temporalMixerType,N_instance,datestr(datetime,30)));
    
end

diary off;

%% Save data
% save([datestr(datetime,30),'.mat']);
fileName = sprintf('%s=%s_%s=%s_N=%d_channel=%s_measurementType_spatial=%s_temporal=%s_observationModel=%s_MCruns=%d_%s.mat',...
    paramName1,mat2str(paramSet1),paramName2,mat2str(paramSet2),N,channelType,spatialMixerType,temporalMixerType,measurementModel,nInstance,datestr(datetime,30));
save(fileName);
%% Plot
% figure;
% for iAlgo = 1:nAlgo
%     semilogy([0:params.maxIter],avgResults.objVal(iAlgo,1:params.maxIter+1));
%     hold on;
% end
% hold off;
% grid on;
% xlabel('number of iterations');
% ylabel('objective function value');
% legend('SCAphase','SC-PRIME');
% % matlab2tikz('complexityIter.tikz');
% 
% figure;
% for iAlgo = 1:nAlgo
%     semilogy([0, avgResults.runtimeArray(iAlgo,1:params.maxIter)],avgResults.objVal(iAlgo,1:params.maxIter+1));
%     hold on;
% end
% hold off;
% grid on;
% xlabel('CPU time (seconds)');
% ylabel('objective function value');
% legend('SCAphase','SC-PRIME');
% % matlab2tikz('complexityTime.tikz');


% legend_str = cell(1,N_L);
% for i_L = 1:N_L
%     legend_str{i_L} = ['L = ',num2str(L_set(i_L))];
% end
%
% for i_P = 1:N_P
%
%     figure(1 + 10*(i_P-1));
%     for i_L = 1:N_L
%         semilogx(lambda_set,mean(errorX{i_P,i_L},1));
%         hold on;
%     end
%     hold off;
%     xlabel('lambda factor');
%     ylabel('NMRMSE(X)');
%     title(sprintf('N=%d, M=%d, I=%d, SNR=%d dB, P=%d',N,M,I,SNR,P_set(i_P)));
%     grid on;
%     legend(legend_str);
%
%     figure(2 + 10*(i_P-1));
%     for i_L = 1:N_L
%         semilogx(lambda_set,mean(errorD{i_P,i_L},1));
%         hold on;
%     end
%     hold off;
%     xlabel('lambda factor');
%     ylabel('NMRMSE(D)');
%     title(sprintf('N=%d, M=%d, I=%d, SNR=%d dB, P=%d',N,M,I,SNR,P_set(i_P)));
%     grid on;
%     legend(legend_str);
%
%     figure(3 + 10*(i_P-1));
%     for i_L = 1:N_L
%         semilogx(lambda_set,mean(errorZ{i_P,i_L},1));
%         hold on;
%     end
%     hold off;
%     xlabel('lambda factor');
%     ylabel('NMRMSE(Z)');
%     title(sprintf('N=%d, M=%d, I=%d, SNR=%d dB, P=%d',N,M,I,SNR,P_set(i_P)));
%     grid on;
%     legend(legend_str);
%
%     figure(4 + 10*(i_P-1));
%     for i_L = 1:N_L
%         semilogx(lambda_set,mean(errorAbsZ{i_P,i_L},1));
%         hold on;
%     end
%     hold off;
%     xlabel('lambda factor');
%     ylabel('NMRMSE(AbsZ)');
%     title(sprintf('N=%d, M=%d, I=%d, SNR=%d dB, P=%d',N,M,I,SNR,P_set(i_P)));
%     grid on;
%     legend(legend_str);
%
%     figure(5 + 10*(i_P-1));
%     for i_L = 1:N_L
%         semilogx(lambda_set,mean(FmeasureZ{i_P,i_L},1));
%         hold on;
%     end
%     hold off;
%     xlabel('lambda factor');
%     ylabel('F-measure Z');
%     title(sprintf('N=%d, M=%d, I=%d, SNR=%d dB, P=%d',N,M,I,SNR,P_set(i_P)));
%     grid on;
%     legend(legend_str);
%
%     figure(6 + 10*(i_P-1));
%     for i_L = 1:N_L
%         semilogx(lambda_set,mean(sparsityZ{i_P,i_L},1));
%         hold on;
%     end
%     hold off;
%     xlabel('lambda factor');
%     ylabel('sparsity Z');
%     title(sprintf('N=%d, M=%d, I=%d, SNR=%d dB, P=%d',N,M,I,SNR,P_set(i_P)));
%     grid on;
%     legend(legend_str);
%
%     figure(7 + 10*(i_P-1));
%     for i_L = 1:N_L
%         semilogx(lambda_set,mean(rmse_sparsityZ{i_P,i_L},1));
%         hold on;
%     end
%     hold off;
%     xlabel('lambda factor');
%     ylabel('RMSE sparsity Z');
%     title(sprintf('N=%d, M=%d, I=%d, SNR=%d dB, P=%d',N,M,I,SNR,P_set(i_P)));
%     grid on;
%     legend(legend_str);
% end



% figure;
% surf(log10(lambdaSet(1:end)),log10(muSet),permute(mean(FmeasureZ{1,1,3},1),[2 3 1]));
% xlabel('\lambda');
% ylabel('\mu');
% zlabel('F-measure(Z)');
% grid on;
% % set(gca,'xscale','log');
% % set(gca,'yscale','log');
%
% figure;
% surf(lambdaSet,muSet,permute(mean(errorX{1,1,3},1),[2 3 1]));
% xlabel('\lambda');
% ylabel('\mu');
% zlabel('MNSE(X)');
% grid on;
% set(gca,'xscale','log');
% set(gca,'yscale','log');