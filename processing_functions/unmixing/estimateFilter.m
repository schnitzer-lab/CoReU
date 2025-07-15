
%   Estimate time-domain convolutional filter w to approximate signal x1(t) 
%   with the reference x2(t):
%       || x1(t) - (w*x2)(t) ||.^2  -> min
%   The duratin of the filter wn < nt ( x1,x2 \in R^{nt x ns} ) serves as a  
%   regualarizer (reduced-rank regresssion). Assuming x1(t), x2(t), w(t) 
%   are stationary, we can effisiently solve it in Fourier domain. The 
%   length of the filter determines its spectral resolution, so the optimal 
%   wn can be deduced from the spectral properties of the x1 and x2, and 
%   can be further optimized with cross-validation (err_loocv).
%   Signals along the second dimention (ns) are processed independently.
%
%   function w = estimateFilter(x1, x2, wn, dn)
%
%   Inputs:
%       x1          - signal (double/single nt x ns) 
%       x2          - reference (double/single nt x ns)
%       wn          - fft windows size (integer) or
%                     fft window itself (nuttallwin(wn) by default)
%       no          - overlap between adjacent fft windows (integer)
%
%   Outputs
%       w           - filter estimates, double/single (wn x ns) array.
%       err_train   - (optional) mean fit ("train") error computed across all windows
%       err_loocv   - (optional) mean leave-one-out CV ("test") error computed across all windows
%
   
%   by Vasily Kruzhilin
%
%   The initial version of this script was inspired by 
%   "DECONVOLUTION OF TWO DISCRETE TIME SIGNALS IN FREQUENCY DOMAIN" 
%   by Dr. Erol Kalkan, P.E. (availible via MATLAB Central)

function [w,err_train,err_loocv] = estimateFilter(x1, x2, wn, no)
    
    if(length(wn) > 1), win = wn; wn = length(wn);
    else, win = nuttallwin(wn); end
    if(nargin < 4), no = round(0.7*wn); end

    assert(all(size(x1) == size(x2)), "estimateFilter: size(x1) ~= size(x2)")
    assert((round(wn)==wn) && wn > 0, "estimateFilter: wn should be a positive int")
    assert((round(no)==no) && no >= 0 && no < wn, "estimateFilter: no should be an int \in [0, wn)")
    assert(wn <= size(x1,1), "estimateFilter: window size is wn too large")
  
    nt = size(x1,1);       
    indxs = (1:(wn-no):(nt-wn)) + (0:(wn-1))';
    ni = size(indxs,2);

    x1grid = reshape(x1(indxs, :), [wn, ni, size(x1,2)]); clear('x1');
    x2grid = reshape(x2(indxs, :), [wn, ni, size(x2,2)]); clear('x2');

    u1all = fft(win.*x1grid, wn); clear('x1grid');
    u2all = fft(win.*x2grid, wn); clear('x2grid');

    s12 = squeeze(mean(u1all.*conj(u2all), 2));
    % s11 = squeeze(mean(u1all.*conj(u1all), 2)); 
    if(nargout < 2), clear('u1all'); end
    s22 = squeeze(mean(u2all.*conj(u2all), 2)); 
    if(nargout < 2), clear('u2all'); end
    
    s = s12./(s22);
    w = real(fftshift(ifft(s),1));

    if(nargout >= 2)
        % shortcut formula for the linear regression LOOCV estimate
        h = 1/ni*(u2all.*conj(u2all))./(mean(u2all.*conj(u2all), 2));  
        err_loocv = squeeze(mean(abs((u1all - s.*u2all)./(1-h)).^2, 1) ./ mean(abs(u1all).^2, 1));
        err_train = squeeze(mean(abs(u1all - s.*u2all).^2, 1) ./ mean(abs(u1all).^2, 1));
    end
end
