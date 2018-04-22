
%% Rx processing params

FFT_OFFSET                    = 4;           % Number of CP samples to use in FFT (on average)
LTS_CORR_THRESH               = 0.8;         % Normalized threshold for LTS correlation

% Repeat the following code for each packet

%% Packet Detection

% Cross correlation of received signal with LTS or STS
% With LTS: Can use the thesis for reference
% With STS: You can read the paper 'Robust Frequency and Timing
% Synchronization for OFDM' by Timothy M. Schmidl and Donald C. Cox

%% CFO estimation and correction


%% Channel estimation


%% Rx payload processing


%% SFO estimation and correction using pilots


%% Phase Error Correction


%% Demodulation

% FEC decoder
Demap_out = demapper(rx_syms,mod_type,1);

% viterbi decoder
rx_data=vitdec(Demap_out,trel,7,'trunc','hard');

% rx_data is the final output corresponding to tx_data, which can be used
% to calculate BER

[number,ber] = biterr(tx_data,rx_data);
