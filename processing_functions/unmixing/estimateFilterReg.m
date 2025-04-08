
%%

function w = estimateFilterReg(x1, x2, wn, dn, eps, noise_levels)
      
    if(nargin < 5), eps = 1; end
    if(nargin < 6), noise_levels = []; end    
    
    if(length(wn) > 1)
        W = wn;
        wn = length(wn);
    else
        % if(mod(wn,2) == 0), wn = wn+1; end
        W = hann(wn);
    end


    nt = length(x1);
        
    indxs = (1:dn:(nt-wn)) + (0:(wn-1))';
    U1fall = fft(W.*x1(indxs), wn);
    U2fall = fft(W.*x2(indxs), wn);

    u12f = mean(U1fall.*conj(U2fall),2);
    u22f = mean(U2fall.*conj(U2fall),2);
    u11f = mean(U1fall.*conj(U1fall),2);
    
    if(isempty(noise_levels))
        u11noise = median(u11f(round(3/8*wn):round(5/8*wn)));
        u22noise = median(u22f(round(3/8*wn):round(5/8*wn)));
        noise_levels = [sqrt(u11noise/wn/mean(W.^2, 1)),...
                        sqrt(u22noise/wn/mean(W.^2, 1))];
    end

    xn1 = randn(size(x1))*noise_levels(1);
    xn2 = randn(size(x2))*noise_levels(2);
    u12noise = mean( fft(W.*xn1(indxs),wn).*conj(fft(W.*xn2(indxs),wn)) ,2);

    s1 = u12f./(u22f+eps*mean(abs(u12noise)));
    s = s1;
   
    % s2 = u11f./(u12f+eps*mean(abs(u11noise)));
    % s = (s1+s2)/2;

    w = real(fftshift(ifft(s)));
end
