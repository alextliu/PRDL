function Y = op_stft(X,window_len,hopsize,fft_len,forward)
%OP_STFT Perform short-time fourier transform on rows of X with
%rectangular window
%
%   Forward map performs a STFT on each row of X, where length of DFT is
%   equal to length of row of X but only first fft_len DFT measurements are
%   taken.
%   
%   Note: consider only the case where both window_len and n2 can be
%   divided by hopsize

[n1,n2] = size(X);
overlap_len = window_len - hopsize; % overlap length
if forward % forward map
    nSection = (n2 + overlap_len) / hopsize;  % no. of short-time sections
    dft_matrix = exp(-2i*pi/n2 * [0:fft_len-1].' * [0:window_len-1]);  % the left-upper block of size fft_len x window_len of n2 x n2 non-normalized DFT matrix
    linear_phase = exp( -2i*pi/n2 * [0:fft_len-1].' * [-overlap_len : hopsize : n2-hopsize] ); 
    % for each channel, pad the signal with zeros of length overlap_len at
    % the beginning and end to remove edge-effects, and then extract
    % short-time sections
    Xsections = buffer(reshape([X, zeros(n1,overlap_len)].',[],1),window_len,overlap_len);
    Y = dft_matrix * Xsections;     % calculate DFT for shifted windowed short-time sections
    Y = reshape(Y,fft_len,nSection,n1).* linear_phase;  % compensate linear phase to shift sections back to the orignial time positions
    Y = reshape(Y,[],n1).';
else % adjoint
    nSection = n2 / fft_len;    % no. of short-time sections
    nTime = nSection * hopsize - overlap_len;   % no. of timeslots
    linear_phase = exp( 2i*pi/nTime * [0:fft_len-1].' * [-overlap_len : hopsize : nTime-hopsize] ); 
    idft_matrix = exp(2i*pi/nTime * [0:window_len-1].' * [0:fft_len-1]);  % the left-upper block of size fft_len x window_len of n2 x n2 non-normalized inverse DFT matrix
    Xsections = reshape(X.', fft_len,nSection,n1) .* linear_phase;  % perform linear phase in frequency domain to shift sections to time zero
    Ysections = idft_matrix * reshape(Xsections,fft_len,[],1);      % perform inverse DFT on each section
    Y = sec_reassemble(Ysections,hopsize,nSection,n1);              % reassemble short-time sections for each channel
    % remove the padded part at the beginning and end
    Y(1:overlap_len,:) = [];
    Y(end-overlap_len+1:end,:) = [];
    Y = Y.';
    
    % implement in MATLAB
%     tmp2 = reshape(X.', fft_len,nSection,n1) .* linear_phase;
%     tmp3 = idft_matrix * reshape(tmp2,fft_len,[],1);
%     tmp4 = reshape(tmp3,window_len,nSection,n1);
%     tmp5 = zeros(nTime+2*overlap_len,nSection,n1);
%     sub1 = repmat( [1:window_len].' + [0:hopsize:(nSection-1)*hopsize], 1,1,n1);
%     sub2 = repmat([1:nSection], window_len,1,n1);
%     sub3 = repmat( reshape([1:n1],1,1,n1), window_len,nSection,1);
%     ind = sub2ind(size(tmp5),sub1,sub2,sub3);
%     tmp5(ind) = tmp4;
%     tmp5(1:overlap_len,:,:) = [];
%     tmp5(end-overlap_len+1:end,:,:) = [];
%     Y = squeeze(sum(tmp5,2)).';
end

end

