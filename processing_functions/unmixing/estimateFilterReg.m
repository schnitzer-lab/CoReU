 
%   Estimate time-domain convolutional filter w (regularized) to approximate 
%   signal x1(t) with the reference x2(t):
%       || x1(t) - (w*x2)(t) ||.^2 + ||\lambda*w||.^2 -> min
%   The duratin of the filter nw < nt ( x1,x2 \in R^{nt x ns} ) serves as a  
%   regualarizer (reduced-rank regresssion), and the L2 regularizer \lambda = 
%   \labmda(x1,x2,nw) is chosen adaptiveley. Assuming x1(t), x2(t), w(t) 
%   are stationary, we can effisiently solve it in Fourier domain. The 
%   length of the filter determines its spectral resolution, so the optimal 
%   nw can be deduced from the spectral properties of the x1 and x2, and the 
%   optimal  value of ridge regularizer eps_ridge is expected to be ~ 1; 
%   both can be further optimized with cross-validation (err_loocv).
%   Signals along the second dimention (ns) are processed independently.
% 
%  function [w,err,loocv] = estimateFilterReg(x1, x2, nw)
%   estimateFilterReg(x1, x2, nw, no, ...
%       eps_ridge, s0abs, fref, eps_tikh, noise2_sd, fnoise)
%
%   Inputs:
%       x1          - signal (double/single nt x ns) 
%       x2          - reference (double/single nt x ns)
%       nw          - fft windows size (integer) or
%                     fft window itself (nuttallwin(nw) by default)
%       no          - (optional) overlap between adjacent fft windows (integer, default 0.7*nw)
%       eps_ridge   - (optional) ridge regularization weight (1 x ns, default 1)
%       s0abs       - (optional) expected value for |s| (nw x ns, default [] - estimated)
%       fref        - (optional) reference frequency to use for s0abs
%                     estimate (default []). Overwritten by s0abs
%       eps_tikh    - (optional) tikhnov regularization weight (1 x ns, default 1e-3)
%       noise2_sd   - (optional) x2 gaussian noise std (1 x ns, default [] - estimated)
%       fnoise      - (optional) low-frequency cut-off used for noise level 
%                     estimation (rel., default 0.75). Overwritten by noise_sd
%
%   Outputs
%       w           - filter estimates, double/single (nw x ns) array.
%       err_train   - (optional) mean fit ("train") error computed across all windows
%       err_loocv   - (optional) mean leave-one-out CV ("test") error computed across all windows
   
%   by Vasily Kruzhilin
%
%   The initial version of this script was inspired by 
%   "DECONVOLUTION OF TWO DISCRETE TIME SIGNALS IN FREQUENCY DOMAIN" 
%   by Dr. Erol Kalkan, P.E. (availible via MATLAB Central)

function [w,err,loocv] = estimateFilterReg(x1, x2, nw, no, ...
    eps_ridge, s0abs, fref, eps_tikh, noise2_sd, fnoise)

    if(length(nw) > 1), win = nw; nw = length(nw);
    else, win = nuttallwin(nw); end

    if(nargin < 4), no = round(0.7*nw); end
    
    if(nargin < 5 || isempty(eps_ridge)), eps_ridge = 1; end
    if(nargin < 6), s0abs = []; end
    if(nargin < 7), fref = []; end

    if(nargin < 8 || isempty(eps_tikh)), eps_tikh = 1e-3; end
    if(nargin < 9), noise2_sd = []; end    
    if(nargin < 10 || isempty(fnoise)), fnoise = 0.75; end    

    assert(all(size(x1) == size(x2)), "estimateFilterReg: size(x1) ~= size(x2)")
    assert((round(nw)==nw) && nw > 0, "estimateFilterReg: nw should be a positive int")  
    assert((round(no)==no) && no >= 0 && no < nw, "estimateFilter: no should be an int \in [0, nw)")
    assert(nw <= size(x1,1), "estimateFilterReg: window size is nw too large")
    assert(all(eps_ridge >= 0), "estimateFilterReg: eps_ridge should be >= 0")
    assert(all(eps_tikh >= 0), "estimateFilterReg: eps_tikh should be >= 0")
    assert(fnoise >= 0 && fnoise < 1, "estimateFilterReg: fnoise shoud be within [0,1)")
  
    nt = size(x1,1);  
    ns = size(x1,2);       
    indxs = (1:(nw-no):(nt-nw)) + (0:(nw-1))';
    ni = size(indxs,2);

    x1grid = reshape(x1(indxs, :), [nw, ni, ns]); clear('x1');
    x2grid = reshape(x2(indxs, :), [nw, ni, ns]); clear('x2');

    u1all = fft(win.*x1grid, nw); clear('x1grid');
    u2all = fft(win.*x2grid, nw); clear('x2grid');

    s12 = mean(u1all.*conj(u2all), 2);
    s22 = mean(u2all.*conj(u2all), 2); 

    % if no noise amplitude provided, estimate it, assuming the signal is
    % white noise dominated at frequencies above fnoise
    if(~isempty(noise2_sd)), s22noiselevel = noise2_sd.^2*nw*mean(win.^2,1);
    else, s22noiselevel = median(s22(round(fnoise/2*nw):round((1-fnoise/2)*nw),:,:,:), 1); end
    assert(all(s22noiselevel(~isnan(s22noiselevel)) >= 0), "estimateFilterReg: s22noiselevel should be >= 0")
    
    % to prevent division by "0" (for example, if x2 was pre-filtered)
    lambda_tikh = eps_tikh*s22noiselevel;
    
    % to use for adaptive lambda_ridge, compute minimally-regularized 
    % solution and use it to estimate x1 "noise" - that is, all the variance
    % that cannot be reconstructed from x2
    s0 = s12./(s22 + lambda_tikh);
    s11noise = (ni)./(ni-1)*mean(abs(u1all - s0.*u2all).^2,2);
    if(nargout < 2), clear('u1all', 'u2all'); end
    
    % to use for adaptive lambda_ridge, compute an estimate s0abs for |s|
    if(isempty(s0abs) &&  isempty(fref)), s0abs = median(abs(s0),1); end
    if(isempty(s0abs) && ~isempty(fref))
        assert(fref >= 0 && fref <= 1, "estimateFilterReg: fref should be \in [0,1]")
        fs = linspace(0,1,length(s0))';
        [~,ind_fref] = min(abs(fs-fref));
        s0abs = max(abs(s0(ind_fref)), median(abs(s0),1));
    end
    assert(all(s0abs(~isnan(s0abs)) >= 0), "estimateFilterReg: s0abs should be >= 0")

    % optimal ridge regression regularizer under bayesian hypothesis
    % u1noise ~ N(0,s11noise) and s ~ N(0,s0abs)
    lambda_ridge = eps_ridge*s11noise./s0abs.^2*1/ni ;
    
    % total L2 regularizer
    lambda = lambda_ridge + lambda_tikh; 
    
    s = s12./(s22 + lambda); 
    w = real(fftshift(ifft(squeeze(s)),1));
    
    if(nargout >= 2)
        % shortcut formula for the linear regression LOOCV estimate
        h = 1/ni*(u2all.*conj(u2all))./(mean(u2all.*conj(u2all), 2) + lambda);  
        loocv = squeeze(mean(abs((u1all - s.*u2all)./(1-h)).^2, 1) ./ mean(abs(u1all).^2, 1));
        err = squeeze(mean(abs(u1all - s.*u2all).^2, 1) ./ mean(abs(u1all).^2, 1));
        if(size(loocv,1) == 1), loocv = loocv'; err = err'; end
    end    
end