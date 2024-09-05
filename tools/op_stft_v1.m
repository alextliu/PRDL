function Y = op_stft_v1(X,window_len,hopsize,fft_len,forward)
%OP_STFT performs short-time fourier transform on rows of X with
%rectangular window
%
%   Version 1: Take the windowed section of length 'window_len' and shift
%   it to time zero. Perform a DFT of length 'fft_len'. If fft_len <
%   window_len, pad zeros at the end of the windowed section.

[n1,n2] = size(X);
overlap_len = window_len - hopsize; % overlap length
if forward % forward map
    Y = reshape( stft([zeros(n1,overlap_len),X,zeros(n1,overlap_len)].',...
        'Window',rectwin(window_len),'Overlaplength',overlap_len,...
        'FFTlength',fft_len,'FrequencyRange','twosided') ./ sqrt(fft_len),...
        [],n1).';       % short-time fourier transform with normalized fft
else % adjoint
    oversample_rate = window_len / hopsize; % oversampling rate
%     nTime = n2 * hopsize / fft_len;  % no. of timeslots
%     sum_window = [ reshape( repmat(1:oversample_rate,[hopsize,1]), [],1 );...
%         oversample_rate * ones(nTime-window_len,1)];    % sum of sliding windows at each time
    Y = istft( reshape(X.', fft_len, [], n1),...
        'Window',rectwin(window_len),'Overlaplength',overlap_len,...
        'FFTlength',fft_len,'FrequencyRange','twosided') .* sqrt(fft_len);  % inverse stft with normalized fft
    Y(1:overlap_len,:) = [];
    Y(end-overlap_len+1:end,:) = [];
    Y = (Y.*oversample_rate).';
end

end

