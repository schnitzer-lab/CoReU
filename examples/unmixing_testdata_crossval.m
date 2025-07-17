
% Generate test signals:
% (signal channel) ch1_meas = ch1_sig + ch1_extra + ch1_noise 
% (reference channel) ch2_meas = ch2_extra + ch2_noise
% where ch1_extra and ch2_extra are coherent, but have minor phase delays
% and slightly varying mixing proportions across the frequencies
[ch1_meas, ch2_meas, ch1_sig, ch1_cont] = ...
    unmixing_generateTestSignals(1e5, 0.05, 0.1, 0, 0);

% Sampling frequency.
% just so the absolute values are more familiar
Fs = 120; 
%%

% run the unmixing for different window lenths and use the built-in error
% estimation to accumulate the errors

% use non-overlaping windows
no = 0;

errors = [];  

window_lenths = (10.^(-1:0.1:2.5));

for wn = round(Fs*window_lenths)
    %%  
    % estimating ch1 contamination via regularized convolutional filter:
    [~, err_train, err_loocv] = estimateFilter(ch1_meas, ch2_meas, wn, no);
   
    errors = [errors; mean(err_train), mean(err_loocv)];
    %%
end
%%
% The optimal window length is about 10s in this case

loglog(window_lenths, errors, '. -'); 
legend(["train", "test"]);

xlabel('window length, s')
ylabel('error, rel')
%%

% We could also use a regularized version of the filter estimation script.
% Note that the explicit regularization only makes a difference when the 
% measured traces ch1_meas and ch2_meas are not long enough (compared to the 
% lenght of the convolution window wn) for regularization to happen
% naturally through averaging
[ch1_meas, ch2_meas, ch1_sig, ch1_cont] = ...
    unmixing_generateTestSignals(1e4, 0.05, 0.1, 0, 0);
%%

% run the unmixing for different values of the regularizer and use the 
% built-in error estimation to accumulate the errors

wn = 10*Fs;
no = 0;
errors = [];  

epss = (10.^(-4:0.1:1));

for eps_ridge = epss
    %%  
    % estimating ch1 contamination via regularized convolutional filter:
    [~, err_train, err_loocv] = estimateFilterReg(ch1_meas, ch2_meas, wn, no, eps_ridge);
    
    errors = [errors; mean(err_train), mean(err_loocv)];
    %%
end
%%

% we can see that the regularizer has a marginal effect and the changes
% in the error amplitude are small, but the optimal 
% value is around eps_ridge ~= 1, as designed internaly
loglog(epss, errors, '.-'); 
legend(["train", "test"]);

xlabel('ridge regularizer, rel')
ylabel('error, rel')