function M = create_cdp_masks(n1,n2,nM,type)
%function M = create_cdp_masks(n1,n2,nM,type)
%
% Creates nM coded diffraction masks of size n1 x n2 and specified type
% (can be 'ternary' or 'complex', or (experimental) 'binary').
%
% Builds on code accompanying Candes, Li & Soltanolkotabi (2014),
%   "Phase Retrieval via Wirtinger Flow: Theory and Algorithms",
% (arXiv: 1407.1065)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is part of the DOLPHIn package (version 1.10)
% last modified: 02/06/2016, A. M. Tillmann
%
% You may freely use and modify the code for academic purposes, though we
% would appreciate if you could let us know (particularly should you find
% a bug); if you use DOLPHIn for your own work, please cite the paper
%
%    "DOLPHIn -- Dictionary Learning for Phase Retrieval",
%    Andreas M. Tillmann, Yonina C. Eldar and Julien Mairal, 2016.
%    http://arxiv.org/abs/1602.02263
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

masktypes = {'binary';'ternary';'complex'};
% 'binary' is experimental; didn't check admissability...

if nargin < 4
    type = masktypes{3};
end

if contains(type,masktypes{1},'IgnoreCase',true)
    %    M = ones(n1,n2,nM);  % Storage for nM masks, each of dim n1 x n2
    while true
        % Sample magnitudes and make masks
        M = randi(2,n1,n2,nM) - 1; % gives 0 entries with prob. 1/2
        if all(sum(abs(M),3) > 0, 'all')   % each entry should be sampled by at least one mask
            break;
        end
    end
elseif contains(type,masktypes{2},'IgnoreCase',true)
    %    M = ones(n1,n2,nM);  % Storage for nM masks, each of dim n1 x n2
    while true
        % Sample magnitudes and make masks
        temp = rand(n1,n2,nM);
        M = (temp <= 0.25) - (temp > 0.75); % gives 0 entries with prob. 1/2
        if all(sum(abs(M),3) > 0, 'all')   % each entry should be sampled by at least one mask
            break;
        end
    end
elseif contains(type,masktypes{3},'IgnoreCase',true)
    M = zeros(n1,n2,nM);  % Storage for nM masks, each of dim n1 x n2
    
    % Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob.
    for ll = 1:nM, M(:,:,ll) = randsrc(n1,n2,[1i -1i 1 -1]); end
    
    % Sample magnitudes and make masks
    temp = rand(n1,n2,nM);
    M = M .* ( (temp <= 0.2)*sqrt(3) + (temp > 0.2)/sqrt(2) );
end