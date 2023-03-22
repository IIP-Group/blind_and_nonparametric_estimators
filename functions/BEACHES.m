% -----------------------------------------------------
%BEACHES   Fast SURE denoiser that assumes noise variance N0 is known.
%   [x_hat,tau_opt,SURE_min,tau,SURE] = BEACHES(y,N0,transform)
%
%   The argument y represents a noisy measurement vector y = x + n, where:
%   - x: noise-free vector that has a sparse representation in some domain
%   - n: white complex Gaussian noise with variance per entry equal to N0
%
%   The argument transform is a string that defines the sparsifying linear 
%   transform that makes x sparse (options: 'none', 'FFT')
%
%   The output x_hat is the estimate of x that results from denoising y 
%   using soft-thresholding, where the threshold is automatically tuned by 
%   minimizing complex SURE given the noise variance N0
%
%   The output tau_opt is the selected denoising threshold, which gives the
%   minimum value of SURE (SURE_min)
%
%   tau and SURE are vectors with (non-optimal) threshold values, and
%   respective values of SURE when applying those thresholds
%
% 2020 (c) ag753@cornell.edu (based on BEACHES algorithm by
% rg548@cornell.edu, studer@cornell.edu, sm2675@cornell.edu, 
% ag753@cornell.edu)
% Edited on March 2023 by ag753@cornell.edu: clean up and comment the code
% for GitHub upload
% -----------------------------------------------------

function [x_hat,tau_opt,SURE_min,tau,SURE] = BEACHES(y,N0,transform)

% length of the vectors:
D = length(y);
% go to the sparse domain:
switch transform
    case 'none'
        y_sparse = y;
    case 'FFT'
        y_sparse = fft(y)/sqrt(D); 
    % more transforms could be added here if needed
    otherwise
        error('transform not implemented')
end
% sort values of |y_sparse|:
y_sparse_sorted = sort(abs(y_sparse),'ascend');
% initialize SURE quantities:
x_hat_sparse = zeros(D,1);
SURE = NaN(D,1);
tau = NaN(D,1);
tau_opt = NaN;
cumsum = 0;
cumsuminv = sum(1./abs(y_sparse(y_sparse~=0))); % from first nonzero index
SURE_min = inf;
tau_low = 0;    
% use SURE to find optimal threshold tau
for dd = find(y_sparse_sorted,1):D % start from first nonzero index
    tau_high = y_sparse_sorted(dd);
    tau(dd) = max(tau_low,min(tau_high,N0/(2*(D-dd+1))*cumsuminv));
    SURE(dd) = cumsum+(D-dd+1)*tau(dd)^2+D*N0-2*N0*(dd-1)-...
        tau(dd)*N0*cumsuminv;
    cumsum = cumsum+y_sparse_sorted(dd).^2;
    cumsuminv = cumsuminv-1/y_sparse_sorted(dd);
    if SURE(dd)<SURE_min
        SURE_min = SURE(dd);
        tau_opt = tau(dd);
    end
    tau_low = y_sparse_sorted(dd);
end    
% compute denoised signal using soft-thresholding with threshold tau_opt:
x_hat_sparse(y_sparse~=0) = (y_sparse(y_sparse~=0)./abs(y_sparse(y_sparse~=0)))...
    .*max(abs(y_sparse(y_sparse~=0))-tau_opt,0); % don't divide by abs(0)
if any(isnan(x_hat_sparse))
    error('denoised vector has NaNs');
end
% go back from the sparse domain to the original domain:
switch transform
    case 'none' 
        x_hat = x_hat_sparse;
    case 'FFT'
        x_hat = ifft(x_hat_sparse)*sqrt(D);
    otherwise
        error('transform not implemented')
end  
SURE_min = SURE_min/D;
SURE = SURE/D;
end


