
 function [tx_data,decoded_data]= MyOfdmReceiver(data);

 
 %% run transmitter code to load sts and lts and other parameters 
 OFDM_TX; 
 
%% Rx processing params

rx_data = raw_rx_data;          % run OFDM tx code to get raw_rx_dec
STS_CORR_THRESH =1/MOD_ORDER;         % Normalized threshold for LTS correlation
LTS_CORR_THRESH =2/MOD_ORDER;         % Normalized threshold for LTS correlation
% Usage: Find all peaks whose magnitude is greater than 0.8 times
% the maximum magnitude after cross correlation (Packet Detection)

% Repeat the following code for each packet

%% Packet Detection

% ideas: Cross correlation of received signal with LTS or use STS to detect the packet?

length_samples= length(rx_data);
sample=16;

sts_corr=xcorr(rx_data,sts_t);
lts_corr=xcorr(rx_data,lts_t);
rx_data_acorr=xcorr(rx_data,rx_data);
sts_corr=sts_corr(length_samples:end);
sts_corr=sts_corr/max(sts_corr);
lts_corr=lts_corr(length_samples:end);
lts_corr=lts_corr/max(lts_corr);
sts_max_corr=find(abs(sts_corr)>STS_CORR_THRESH);
sts_start=sts_max_corr(1);
lts_max_corr=find(abs(lts_corr)>LTS_CORR_THRESH);
lts_start=lts_max_corr(1);
f_data=rx_data(sts_start:end);
figure(4)
plot(abs(f_data))
title('Data after Packet Detection')
xlabel('Symbol Index')

%%
figure(7)
subplot(2,1,1)
plot(abs(sts_corr))
title('STS Correlation in Presence of CFO')
xlabel('Symbol Index')
xlim([1 1000])
yline(STS_CORR_THRESH,'--k','label','Correlation threshold','LabelHorizontalAlignment','left')
subplot(2,1,2)
plot(abs(lts_corr))
title('LTS Correlation in Presence of CFO')
xlabel('Symbol Index')
xlim([1 1000])
yline(STS_CORR_THRESH,'--k','label','Correlation threshold','LabelHorizontalAlignment','left')
%%
% 
% while(sample< length_samples)
% 	% sts correlation
% 
% 	output(sample) = (rx_data(sample-length(sts_t)+(1:length(sts_t)))* rx_data(sample+(1:length(sts_t)))')./rms(rx_data(sample+(1:length(sts_t))));
% 
% 	sample= sample+1;
% end
% 
% output= output./max(abs(output));
% 
% figure;
% plot((abs(output))); 
% hold on 
% plot(abs(rx_data)); 
% hold off

% Output: Single packet extracted from rx_data
% with knowledge of preamble (LTS) indices and payload vector indices


%% CFO estimation and correction
% Use two copies of LTS for cross-correlation (Reference: Thesis)
index1=length(sts_t)*30+2*CP_LEN+1;
rx_lts1=f_data(index1:index1+N_SC-1);
rx_lts2=f_data(index1+N_SC:index1+2*N_SC-1);
cfo=0;
for n=1:N_SC
   cfo=cfo+imag(rx_lts2(n)/rx_lts1(n))/(2*pi*N_SC);
end
cfo=cfo/N_SC;               %Taking mean
%%
rx_data_cfo_offset=find(rx_data>0);
cfo_offset=rx_data_cfo_offset(1);
rx_data_wocfo=zeros(size(rx_data));
rx_data_wocfo(cfo_offset:end)=rx_data(cfo_offset:end).*exp(-1j*2*pi*(0:length(rx_data(cfo_offset:end))-1)*cfo);

%%Redoing Packet Detection after CFO removal
sts_corr2=xcorr(rx_data_wocfo,sts_t);
sts_corr2=sts_corr2(length_samples:end);
sts_corr2=sts_corr2/max(sts_corr2);
sts_max_corr2=find(abs(sts_corr2)>STS_CORR_THRESH);
sts_start2=sts_max_corr2(1);
f_data2=rx_data_wocfo(sts_start2:end);      %Data with STS and LTS after CFO correction
figure(5)
plot(abs(f_data2))
title('Data after Packet Detection and CFO Removal')
xlabel('Symbol Index')

%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help
index2=length(sts_t)*30+2*CP_LEN+2*length(lts_t)+1;
rx_mat_wcp=reshape(f_data2(index2:index2+(CP_LEN+N_SC)*N_OFDM_SYMS-1),CP_LEN+N_SC,N_OFDM_SYMS);
rx_mat_nocp=rx_mat_wcp(CP_LEN+[1:N_SC],:);
% Output: CP free payload matrix of size (N_SC * N_OFDM_SYMS)


%% FFT
% Refer to IFFT perfomed at TX
payload=fft(rx_mat_nocp, N_SC, 1);
% Output: Symbol matrix in frequency domain of same size


%% Channel estimation and correction
% Use the two copies of LTS and find channel estimate (Reference: Thesis)
% Convert channel estimate to matrix form and equlaize the above matrix
rx_lts1f=fft(f_data2(index1:index1+N_SC-1));
rx_lts2f=fft(f_data2(index1+N_SC:index1+2*N_SC-1));
temp1=(rx_lts2f./lts_f);
temp2=(rx_lts1f./lts_f);
H=(temp1+temp2)/2;
H=H(:,SC_IND_DATA)';
% H_dash=inv(H'*H)*H';
% H_dash=1./H;
% H_dash=H_dash(:,SC_IND_DATA)';
% Output : Symbol equalized matrix in frequency domain of same size


%% Advanced topics: 
%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.

% Output: Symbol equalized matrix with pilot phase correction applied


%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error
pfo=zeros(size(pilots));
mean_cfo = mean(angle(payload(SC_IND_PILOTS,:)./pilots),1);
cfo_fine = repmat(unwrap(mean_cfo),N_SC,1);
payload = payload.*(exp(-1i*cfo_fine));

% Output: Symbol equalized matrix with pilot phase correction applied
% Remove pilots and flatten the matrix to a vector rx_syms

payload_syms_mat = payload(SC_IND_DATA,:);
payload_syms_mat=payload_syms_mat./H;
rx_syms=reshape(payload_syms_mat,1,numel(payload_syms_mat));
%% Demodulation

figure(6);
scatter(real(rx_syms), imag(rx_syms),'filled');
title(' Signal Space of received bits');
xlabel('I'); ylabel('Q');

% FEC decoder
Demap_out = demapper(rx_syms,MOD_ORDER,1);

% viterbi decoder
decoded_data = vitdec(Demap_out,trel,7,'trunc','hard');

% decoded_data is the final output corresponding to tx_data, which can be used
% to calculate BER
%% Bit error
[num,ratio]=biterr(tx_data,decoded_data)