
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

% looking into the signals, 'measured' ch1_meas contains signal component at 
% ~8, ~20, and ~55Hz, together with contamination at ~5Hz (that effectiveley 
% masks the ~8Hz signal), 15Hz, 40Hz and 45Hz. 'Measured' ch2_meas has only
% the aforementioned contamination componens. The flat areas of the PSDs
% indicate the white noise floor limits 
plt.tracesComparison([ch1_meas, ch2_meas, ch1_sig, ch1_cont], ...
    'spacebysd', 3, 'labels', [...
    "signal+contamination ('measured')", "reference ('measured')", ...
    "clean signal (ground truth)", "contamination (ground truth)"], ...
    'fps', Fs, 'fw', 0.25)

set(gca, 'ylim',[1*1e-6,100])
subplot(2,1,1)
set(gca, 'xlim', [0,8])
%%
% As mentioned, ch1_meas and ch2_meas coherently share some signals, however
% the amplitude and phase are different across frequencies

[Cxy, ~] = mscohere(ch1_meas, ch2_meas,round(2*Fs),round(1.5*Fs),[], Fs);
[Pxy, fs] = cpsd(ch1_meas, ch2_meas(1:end), round(2*Fs), round(1.5*Fs),[], Fs);

plot(fs, Cxy, 'linewidth', 1); 
ylabel('coherence (voltage and reference channel)')

yyaxis right
a = angle(Pxy);
a((Cxy < 0.1)) = NaN;
% plot(F,a/pi, 'linewidth', 2);
plot(fs,a./fs, 'linewidth', 2);
ylabel('relative phase (x \pi)')
yline(0)
grid on
ylim([-0.1,0.1])
%%

% sliding filter estimation time window length
wn = round(2*Fs);

% estimating ch1 contamination via convolutional filter:
[w,~,~] = estimateFilter(ch1_meas, ch2_meas, wn, no);
ch1_cont_filt = conv(ch2_meas, w, 'same');
ch1_sig_filt = ch1_meas - ch1_cont_filt;

% estimating ch1 contamination via simple linear regression:
ch1_cont_reg = ch2_meas*(ch2_meas\ch1_meas);
ch1_sig_reg = ch1_meas - ch1_cont_reg;

%%

% Comparing the unmixed signlas, we see that both filtering and simple
% regression are impoving the measurement, however, the regression is
% unable to fully unmix the ~5, 15 or 45Hz components of the contamination
% due to small phase shifts and differences in relative contribution. Also
% note that the regression unmixing increases the white-noise level (flat
% parts of the PSD), to a point where the ~55Hz signal peak is completely
% masked. The filtering succsesfuly removed all the contamination
% components and preserved the initial noise level of the ch1. The
% filtering procedure is however still limited by the noise level in the
% ch2, that is it cannot completeley unmix the contaminations that were
% below the ch2 noise level (and thus effectiveley not measured in ch2),
% but above ch1 noise level.
plt.tracesComparison([ch1_sig, ch1_sig_filt, ch1_sig_reg], ...
    'spacebysd', 3, 'labels', [...
    "clean signal (ground truth)", "clean signal (filtering)", ...
    "clean signal (regression)"],...
    'fps', Fs, 'fw', 0.25)

set(gca, 'ylim',[1*1e-6,100])
subplot(2,1,1)
set(gca, 'xlim', [0,8])
%%

% We could also use a regularized version of the filter estimation script.
% Note that the explicit regularization only makes a difference when the 
% measured traces ch1_meas and ch2_meas are not long enough (compared to the 
% lenght of the convolution window wn) for regularization to happen
% naturally through averaging
[ch1_meas, ch2_meas, ch1_sig, ch1_cont] = ...
    unmixing_generateTestSignals(1e3, 0.05, 0.1, 0, 0);
%%

% let's force to overlap of the filter estimation time windows to be zero
% to reduce the effects of the averaging
no = 0;

% estimating ch1 contamination via convolutional filter:
[w0,~,~] = estimateFilter(ch1_meas, ch2_meas, wn, no);
ch1_cont_filt = conv(ch2_meas, w0, 'same');
ch1_sig_filt = ch1_meas - ch1_cont_filt;

% estimating ch1 contamination via regularized convolutional filter:
[wr,~,~] = estimateFilterReg(ch1_meas, ch2_meas, wn, no, 1, [], fref);
ch1_cont_filtr = conv(ch2_meas, wr, 'same');
ch1_sig_filtr = ch1_meas - ch1_cont_filtr;

% estimating ch1 contamination via simple linear regression:
ch1_cont_reg = ch2_meas*(ch2_meas\ch1_meas);
ch1_sig_reg = ch1_meas - ch1_cont_reg;
%%

% The unregularized estimation of the contamination suffers from the
% overfitting of the ch2 noise to the strong true signal in ch1 at ~20Hz,
% but the regularized one does not
plt.tracesComparison([ch1_cont, ch1_cont_reg, ch1_cont_filt, ch1_cont_filtr], ...
    'spacebysd', 3, 'labels', [...
    "contamination (ground truth)", "contamination (regression)", ...
    "contamination (filtering)", "contamination (filtering reg)" ], ...
    'fps', Fs, 'fw', 0.5)

set(gca, 'ylim',[1*1e-8,100])
subplot(2,1,1)
set(gca, 'xlim', [0,8])
%%

% Howerver, the ch1_signal estimation remains virtualy the same as the true
% signal at 20Hz overpowers the overfitted component
plt.tracesComparison([ch1_sig, ch1_sig_reg, ch1_sig_filt, ch1_sig_filtr], ...
    'spacebysd', 3, 'labels', [...
    "clean signal (ground truth)", "clean signal (regression)", ...,
    "clean signal (filtering)", "clean signal (filtering reg)"],...
    'fps', Fs, 'fw', 0.5)

set(gca, 'ylim',[1*1e-6,1])
subplot(2,1,1)
set(gca, 'xlim', [0,8])
