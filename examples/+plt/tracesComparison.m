function tracesComparison(traces, varargin)
    
    options = parseInputs(varargin{:});
    
    if(options.fw ~= 0)
        options.nw = options.fw*size(traces,1)/options.fps/2;
    end
    options.nw = max([options.nw, 1.25]); % min value matlab's pmtm accepts
    
    spacings = [0, std(traces(:, 2:end)) + std(traces(:, 1:(end-1)))];
    
    if(options.spectra)
        if(options.horizontal) subplot(1,2,1)
        else subplot(2,1,1) 
        end
    end
    ts = (0:(size(traces,1)-1))/options.fps + options.t0;
    xs = traces- mean(traces)*options.nomean + cumsum(spacings.*options.spacebysd); 
    plot(ts, xs*options.x_plot_scale, 'LineWidth', options.linewidth); xlim([min(ts), max(ts)]);
%     if(~isempty(options.labels)) legend(options.labels, 'Interpreter', 'none','FontSize', 6); end
    title("Time trace"); xlabel("t, s"); ylabel('x'); grid(); 
    
    if(options.spectra)
        if(options.horizontal) subplot(1,2,2)
        else subplot(2,1,2) 
        end

        z0 = pmtm(traces - mean(traces)*options.nomean_psd, options.nw);
        fs = linspace(0, options.fps/2, size(z0,1)); 
        norm0 = 1/pi * options.fps/2; % sum(traces.^2)/length(traces) / sum(z0/norm0*mean(diff(fs))) == 1 % ylabel("[x^2]/Hz")
        % norm0 = options.fps^2/(2*pi*size(traces,1)); % sum(traces.^2*1/options.fps) / sum(z0/norm0*mean(diff(fs))) == 1 % ylabel("[x^2]·s/Hz")

        z0(fs < options.f0,:) = NaN;
        z0(end,:) = NaN;

        semilogy(fs, z0/norm0, 'LineWidth', options.linewidth);  grid(); 
        if(~isempty(options.labels)) legend(options.labels, 'Interpreter', 'none','FontSize', 6); end
        title("PSD"); xlabel("f, Hz"); ylabel("[x^2]/Hz")
        xlim([min(fs), max(fs)]); 
        ylim([0.9,2].*[...
            min(z0(fs >= options.f0,:)/norm0, [], 'all'), ...
            max(z0(fs >= options.f0,:)/norm0, [], 'all')])
    end
end

function [options,p] = parseInputs(varargin)
    

    p = inputParser();
    p.addParameter('spectra', true);

    p.addParameter('fw', 0);
    p.addParameter('nw', 1.25);
    
    p.addParameter('fps', 1);
    p.addParameter('x_plot_scale', 1);
    p.addParameter('labels', []);
    p.addParameter('nomean', true);
    p.addParameter('nomean_psd', true);
    p.addParameter('spacebysd', false);
    p.addParameter('horizontal', false);
%     options.colors = [];
    
    p.addParameter('f0', 0)
    p.addParameter('t0', 0);
    
    p.addParameter('linewidth', 1);

    p.parse(varargin{:});
    options = p.Results;
end