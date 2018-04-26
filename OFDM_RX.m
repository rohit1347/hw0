
%% Rx processing params

rx_data = raw_rx_dec;          % run OFDM tx code to get raw_rx_dec
LTS_CORR_THRESH = 0.8;         % Normalized threshold for LTS correlation

% Repeat the following code for each packet

%% Packet Detection
% Cross correlation of received signal with LTS


%% CFO estimation and correction
% Use two copies of LTS for cross-correlation


%% CP Removal


%% FFT


%% Channel estimation and correction


%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.
% Extract the pilot tones and "equalize" them by their nominal Tx values

%% Phase Error Correction
% Extract the pilots and calculate per-symbol phase error


%% Demodulation

% FEC decoder
Demap_out = demapper(rx_syms,mod_type,1);

% viterbi decoder
rx_data=vitdec(Demap_out,trel,7,'trunc','hard');

% rx_data is the final output corresponding to tx_data, which can be used
% to calculate BER

[number,ber] = biterr(tx_data,rx_data);
