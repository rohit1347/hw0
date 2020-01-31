
% function [tx_data,decoded_data]= MyOfdmReceiver(data);


%% run transmitter code to load sts and lts and other parameters
OFDM_TX;
% rx_data=raw_rx_data;
raw_rx_data=load('packet_2_QPSK.mat');

%% Rx processing params
sts_pkt_dtt = 1;
sts_cfo_flag = 1;
pilots_cfo_cflag=1;
sfo_flags=1;

rx_data = raw_rx_data.data;          % run OFDM tx code to get raw_rx_dec
STS_CORR_THRESH =1/MOD_ORDER;         % Normalized threshold for LTS correlation
LTS_CORR_THRESH =1/MOD_ORDER;         % Normalized threshold for LTS correlation
% Usage: Find all peaks whose magnitude is greater than 0.8 times
% the maximum magnitude after cross correlation (Packet Detection)

% Repeat the following code for each packet
% %% Packet detection using auto corr
length_samples= length(rx_data);
sample=16;
%% Packet detection using auto-correlation

rx_data_acorr=zeros(1,length_samples);

for m=sample+1:length_samples-sample
    rx_data_acorr(m-sample)=rx_data(m:m+sample-1)*rx_data(m-sample:m-1)';
end

rx_data_acorr=rx_data_acorr/max(rx_data_acorr);
sts_max_acorr=find(abs(rx_data_acorr)>STS_CORR_THRESH);
sts_acorr_start=sts_max_acorr(1);


%% Packet Detection using xcorr

% ideas: Cross correlation of received signal with LTS or use STS to detect the packet?
sts_corr=xcorr(rx_data,sts_t);
lts_corr=xcorr(rx_data,lts_t);

sts_corr=sts_corr(length_samples:end);
sts_corr=sts_corr/max(sts_corr);
lts_corr=lts_corr(length_samples:end);
lts_corr=lts_corr/max(lts_corr);
sts_max_corr=find(abs(sts_corr)>STS_CORR_THRESH);
sts_start=sts_max_corr(1);
lts_max_corr=find(abs(lts_corr)>LTS_CORR_THRESH);
lts_start=lts_max_corr(1);
if sts_pkt_dtt>0
    f_data=rx_data(sts_start:end);
else
    f_data=rx_data(lts_start:end);
end
%%
%%
figure(8)
subplot(2,1,1)
plot(abs(rx_data_acorr))
title('Auto-correlation')
axis([1 1000 0 1])
yline(LTS_CORR_THRESH,'--k','label','Correlation threshold','LabelHorizontalAlignment','left')
grid on
subplot(2,1,2)
plot(abs(sts_corr))
title('Cross-correlation')
sgtitle('Packet detection using various correlation methods')
axis([1 1000 0 1])
yline(LTS_CORR_THRESH,'--k','label','Correlation threshold','LabelHorizontalAlignment','left')
grid on
%% STS xcorr packet detection
figure(4)
plot(abs(f_data))
title('Data after Packet Detection')
xlabel('Symbol Index')

%%
figure(7)
subplot(2,1,1)
plot(abs(sts_corr))
title(['STS Correlation with CFOFLAG = ',num2str(CFO_FLAG)])
xlabel('Symbol Index')
xlim([1 1000])
yline(STS_CORR_THRESH,'--k','label','Correlation threshold','LabelHorizontalAlignment','left')
subplot(2,1,2)
plot(abs(lts_corr))
title('LTS Correlation in Presence of CFO')
xlabel('Symbol Index')
xlim([1 1000])
yline(LTS_CORR_THRESH,'--k','label','Correlation threshold','LabelHorizontalAlignment','left')
%%


% Output: Single packet extracted from rx_data
% with knowledge of preamble (LTS) indices and payload vector indices


%% CFO estimation and correction
% Use two copies of LTS for cross-correlation (Reference: Thesis)
rx_data_cfo_offset=find(abs(rx_data)>0);
cfo_offset=rx_data_cfo_offset(1);

if sts_pkt_dtt>0
    index1=length(sts_t)*30+2*CP_LEN+1;
else
    index1=1;
end
cfo=0;
if sts_cfo_flag>0
    
    rx_sts30 = rx_data(lts_start-(2*CP_LEN+sample):lts_start-2*CP_LEN);
    rx_sts29 = rx_data(lts_start-(3*CP_LEN+sample):lts_start-3*CP_LEN);
    for n=1:sample
        cfo=cfo+imag(rx_sts30(n)/rx_sts29(n))/(2*pi*sample);
    end
    cfo=cfo/sample;               %Taking mean
end
rx_data_wocfo=zeros(size(rx_data));
rx_data_wocfo(cfo_offset:end)=rx_data(cfo_offset:end).*exp(-1j*2*pi*(0:length(rx_data(cfo_offset:end))-1)*cfo);

if sts_pkt_dtt>0
    f_data_wocfo=rx_data_wocfo(sts_start:end);
else
    f_data_wocfo=rx_data_wocfo(lts_start:end);
end
cfo=0;
rx_lts1=f_data_wocfo(index1:index1+N_SC-1);
% rx_lts1=f_data_wocfo(index1-CP_LEN:index1-CP_LEN+N_SC-1);
rx_lts2=f_data_wocfo(index1+N_SC:index1+2*N_SC-1);
% rx_lts2=f_data_wocfo(index1+N_SC-CP_LEN:index1+2*N_SC-1-CP_LEN);

for n=1:N_SC
    cfo=cfo+imag(rx_lts2(n)/rx_lts1(n))/(2*pi*N_SC);
    %             cfo=cfo+imag(rx_lts2(n)/lts_t(n))/(2*pi*N_SC);
    %             %Uncomment for CFOo estimation using xcorr
end
cfo=cfo/N_SC;       % Taking mean
%     end
rx_data_wocfo(cfo_offset:end)=rx_data_wocfo(cfo_offset:end).*exp(-1j*2*pi*(0:length(rx_data_wocfo(cfo_offset:end))-1)*cfo);

figure(9)
plot(unwrap(angle(rx_lts1)))
hold on
plot(unwrap(angle(rx_lts2)))
plot(unwrap(angle(rx_data_wocfo(sts_start+index1-1:sts_start+index1+N_SC-2))))
plot(unwrap(angle(rx_data_wocfo(sts_start+index1+N_SC-1:sts_start+index1+2*N_SC-2))))
legend('wcfo\_lts1','wcfo\_lts2','wocfo\_lts1','wocfo\_lts2')
hold off
%%
%%Redoing Packet Detection after CFO removal
if sts_pkt_dtt>0
    sts_corr2=xcorr(rx_data_wocfo,sts_t);
    sts_corr2=sts_corr2(length_samples:end);
    sts_corr2=sts_corr2/max(sts_corr2);
    sts_max_corr2=find(abs(sts_corr2)>STS_CORR_THRESH);               % Commented since not using sts correlation
    sts_start2=sts_max_corr2(1);
    f_data2=rx_data_wocfo(sts_start2:end);      %Data with STS and LTS after CFO correction
else
    f_data2= rx_data_wocfo(lts_start:end);                                          % Comment if not using lts corr for pkt detection
end
figure(5)
plot(abs(f_data2))
title('Data after Packet Detection and CFO Removal')
xlabel('Symbol Index')

%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help
if sts_pkt_dtt>0
    index2=length(sts_t)*30+2*CP_LEN+2*length(lts_t)+1;   % in case of using STS corr
else
    index2=2*length(lts_t)+1;       %In case of using LTS corr
end
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
temp1=(rx_lts2f./fft(lts_t));
temp2=(rx_lts1f./fft(lts_t));
H=(temp1+temp2)/2;
% H=temp1;        % If only one channel estimate is to be used
% H=temp2;        % If only one channel estimate is to be used
H_full=H.';
H=H(:,SC_IND_DATA).';

% H_dash=inv(H'*H)*H';
% H_dash=1./H;
% H_dash=H_dash(:,SC_IND_DATA)';
% Output : Symbol equalized matrix in frequency domain of same size

payload=payload./H_full;
%% Advanced topics:
%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.
if sfo_flags>0
    sfo=angle(payload(SC_IND_PILOTS,:)./pilots);
    
    sfo_slope=zeros(length(SC_IND_PILOTS),1);
    for pnum=1:length(SC_IND_PILOTS)
        sfo(pnum,:)=unwrap(sfo(pnum,:)) ;
        %    sfo_slope(pnum)=(sfo(pnum,1)+sfo(pnum,end))/(N_OFDM_SYMS-1);
    end
    mean_sfo_intercept=mean(sfo(:,1));
    % mean_sfo_slope=mean(sfo_slope);
    mean_sfo_slope=(mean(sfo(:,end))-mean(sfo(:,1)))/(N_OFDM_SYMS-1);
    mean_sfo_vec=repmat(mean_sfo_slope,1,N_OFDM_SYMS).*(1:N_OFDM_SYMS)+mean_sfo_intercept;
    mean_sfo_mat=repmat(mean_sfo_vec,N_SC,1);
    sfo_rec=exp(-1i*mean_sfo_mat);
    payload=payload.*sfo_rec;
    
    figure(10)
    plot(sfo(1,:))
    hold on
    plot(mean_sfo_vec)
    plot(sfo(2,:))
    plot(sfo(3,:))
    plot(sfo(4,:))
    plot(mean(sfo))
    legend('sfo1','rec','sfo2','sfo3','sfo4','mean')
    hold off
end
% Output: Symbol equalized matrix with pilot phase correction applied


%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error
if pilots_cfo_cflag>0
    pfo=zeros(size(pilots));
    mean_cfo = mean(angle(payload(SC_IND_PILOTS,:)./pilots),1);
    cfo_fine = repmat(unwrap(mean_cfo),N_SC,1);
    payload = payload.*(exp(-1i*cfo_fine));
end

% Output: Symbol equalized matrix with pilot phase correction applied
% Remove pilots and flatten the matrix to a vector rx_syms

payload_syms_mat = payload(SC_IND_DATA,:);
% payload_syms_mat=payload_syms_mat./H;
rx_syms=reshape(payload_syms_mat,1,numel(payload_syms_mat));
%% Demodulation

% FEC decoder
Demap_out = demapper(rx_syms,MOD_ORDER,1);

% viterbi decoder
decoded_data = vitdec(Demap_out,trel,7,'trunc','hard');

% decoded_data is the final output corresponding to tx_data, which can be used
% to calculate BER
%% Bit error
tx_data=load('tx_data_QPSK.mat');
tx_data=tx_data.data;
[num,ratio]=biterr(tx_data,decoded_data)

h=figure(6);
subplot(1,2,1)
scatter(real(rx_syms), imag(rx_syms),'filled');
xline(0,'k--')
yline(0,'k--')
title(' Signal Space of received symbols');
xlabel('I'); ylabel('Q');
subplot(1,2,2)
scatter(real(tx_syms), imag(tx_syms),'filled');
title(' Signal Space of transmitted symbols');
xlabel('I'); ylabel('Q');
sgtitle(['Modulation Order: ',num2str(MOD_ORDER),' | BER: ',num2str(ratio)])
saveas(h, 'sspace_16QAM_txrx.pdf')
% ! pdfcrop sspace_16QAM_txrx.pdf sspace_16QAM_txrx.pdf