% function w = limitFilter(w, varargin)
%
%   Inputs:
%       w        - time-domail filter to limit, double/single (wn x 1) column
%       varargin - see parseInputs(varargin) below for the list of 
%                  parameters used to limit the filter
%
%   Outputs
%       w        - time-domail filter to limit with limits applied, 
%                  double/single (wn x 1) column

function w = limitFilter(w, varargin)
    
    assert(iscolumn(w), "limitFilter: w should be a column-vector")
    options = parseInputs(varargin{:});
    
    wn = length(w);

    if(isfinite(options.max_delay))
        if(mod(wn, 2))
            phase_delay = [linspace(0, options.max_delay*pi, wn/2), pi, ...
                     linspace(options.max_delay*pi, 0, wn/2)];
        else
            phase_delay = [linspace(0, options.max_delay*pi, wn/2), ...
                     linspace(options.max_delay*pi, 0, wn/2)];
        end

        options.max_phase = min(phase_delay, options.max_phase)';
        options.max_phase(options.max_phase>pi) = pi;
    end

    s = fft(ifftshift(w));
    fs = linspace(0,1,length(s))';
    
    if(~isempty(options.fref))
        [~,ind_fref] = min(abs(fs-options.fref));
        options.max_amp = options.max_amp_rel*max(abs(s(ind_fref)), median(abs(s)));
    end
    
    if(~isempty(options.max_amp))
        % limit max spectral amplitude for fs < flim_max, but preserve the phase
        % note: does not commute with max_phase
        a_large = (abs(s) > abs(options.max_amp));
        in_flim = (fs<=options.flim_max & fs<=0.5) | (flip(fs)<=options.flim_max & flip(fs)<=0.5);
        s(a_large&in_flim) = options.max_amp*(s(a_large&in_flim)./abs(s(a_large&in_flim)));
    end

    % project s(f) on exp(i*max_phase(f)) direction if phase is too large
    % note: does not commute with max_amp
    phase_too_large_p = angle(s) >  options.max_phase;
    phase_too_large_n = angle(s) < -options.max_phase;
    s(phase_too_large_p) = abs(s(phase_too_large_p)).*...
        max(cos(angle(s(phase_too_large_p))-options.max_phase(phase_too_large_p)), 0).*...
        exp(1.i*options.max_phase(phase_too_large_p));
    s(phase_too_large_n) = abs(s(phase_too_large_n)).*...
        max(cos(-angle(s(phase_too_large_n))-options.max_phase(phase_too_large_n)), 0).*...
        exp(-1.i*options.max_phase(phase_too_large_n));




    w = real(fftshift(ifft(s)));
end
%%

function [options,p] = parseInputs(varargin)
    
    p = inputParser();
    isinrange = @(x,a,b) isnumeric(x)&(x>=a)&(x<=b) ;

    % maximimal spectral amplitude of a filter
    p.addParameter('max_amp', [], @(x)isinrange(x,0,Inf)) 
    % maximal frequency to apply amplitude correction, \in [0,0.5]
    p.addParameter('flim_max', 0.5, @(x)isinrange(x,0,0.5))
    
    % reference frequency to determine max_amp, overwrites max_amp, \in [0,0.5]
    p.addParameter('fref',[], @(x)isinrange(x,0,0.5)) 
    % scaling of the options.fref amplitude determine max_amp
    p.addParameter('max_amp_rel',1.2, @(x)isinrange(x,0,Inf)) 

    % max phase delay (uniform across frequencies), [-\pi,pi]
    p.addParameter('max_phase',pi,@(x)isinrange(x,-pi,pi))
    % max time delay (biderectional; relative: max_delay(s)*fps(Hz)), overwrites options.max_phase
    p.addParameter('max_delay',Inf,@(x)isinrange(x,0,Inf))
    
    p.parse(varargin{:});
    options = p.Results;
end