function [x1,x2,x11,x12] = generateTestUnmixingSignal(nt, sd1, sd2, sd12, nslow)

    % nt = 1e5;
    
    u1 = 2*exp(-((1:nt)' - nt/20).^2./(nt/30).^2) + ...
        0.5*exp(-((1:nt)' - nt/8).^2./(nt/100).^2) + ...
        0.25*exp(-((1:nt)' - nt/3).^2./(nt/30).^2) + ...
        0.5*exp(-((1:nt)' - 6*nt/16).^2./(nt/1000).^2);

    u1s = 0.5*exp(-((1:nt)' - nt/15).^2./(nt/60).^2) + ...
        1*exp(-((1:nt)' - nt/6).^2./(nt/150).^2) + ...
        0.05*exp(-((1:nt)' - 0.45*nt).^2./(nt/150).^2) ;
    
    u2a = 2*exp(-((1:nt)' - nt/20).^2./(nt/30).^2).*exp(1i*pi/8) + ...
        0.6*exp(-((1:nt)' - nt/8).^2./(nt/100).^2).*exp(1i*pi/16) + ...
        0.2*exp(-((1:nt)' - nt/3).^2./(nt/30).^2).*exp(-1i*pi/16) + ...
        0.6*exp(-((1:nt)' - 6*nt/16).^2./(nt/1000).^2).*exp(-1i*pi/8) ;


    u2b = 2*exp(-((1:nt)' - nt/20).^2./(nt/30).^2).*exp(-1i*pi/16) + ...
        0.4*exp(-((1:nt)' - nt/8).^2./(nt/100).^2).*exp(-1i*pi/32) + ...
        0.3*exp(-((1:nt)' - nt/3).^2./(nt/30).^2).*exp(+1i*pi/32) + ...
        0.4*exp(-((1:nt)' - 6*nt/16).^2./(nt/1000).^2).*exp(+1i*pi/16) ;
    
    z = fft(randn(size(u1)));
    z1s = fft(randn(size(u1)));
    
    
    x1 = (real(fft(u1.*z + u1s.*z1s)));
    x11 = (real(fft(u1s.*z1s)));
    x12 = (real(fft(u1.*z)));

    alpha = sin((1:nt)/nt*2*pi*nslow)';
    x2 = (real(fft((alpha.*u2a + (1-alpha).*u2b).*z)));
    
    x11 = x11./(std(x1)+eps);
    x12 = x12./(std(x1)+eps);
    x1 = x1./(std(x1)+eps);
    x2 = x2./(std(x2)+eps);

    shared_noise = sd12*randn(size(x1)); 
    x1 = x1 + shared_noise + sd1*randn(size(x1));
    x2 = x2 + shared_noise + sd2*randn(size(x2));
end