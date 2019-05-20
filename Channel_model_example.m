%% Transmission over Channel Model with Delay Profile CDL-D
% Transmit waveform through a Clustered Delay Line (CDL) channel model with 
% delay profile CDL-D from TR 38.901 Section 7.7.1.
% 
% Define the channel configuration structure using an |nrCDLChannel| System 
% object. Use delay profile CDL-D, a delay spread of 10 ns, and UT velocity of 
% 15 km/h:

clc; close all;
v = 200;                    % UT velocity in km/h
fc = 3.5e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd = (v*1000/3600)/c*fc;     % UT max Doppler frequency in Hz
cdl = nrCDLChannel;
%cdl.SampleDensity = 64;
cdl.DelaySpread = 10e-6;
cdl.CarrierFrequency = fc;
cdl.MaximumDopplerShift = fd;

%% 
% Configure the transmit array as [M N P Mg Ng] = [2 2 2 1 1], representing 
% 1 panel (Mg=1, Ng=1) with a 2-by-2 antenna array (M=2, N=2) and P=2 polarization 
% angles. Configure the receive antenna array as [M N P Mg Ng] = [1 1 2 1 1], 
% representing a single pair of cross-polarized co-located antennas.
NTx = 64; NRx = 64;
cdl.TransmitAntennaArray.Size = [sqrt(NTx) sqrt(NTx) 1 1 1];
cdl.ReceiveAntennaArray.Size = [2 2 1 4 4];
cdl.ReceiveAntennaArray.ElementSpacing = [0.5 0.5 10 10];

%% 
% Create a random waveform of 1 subframe duration with 8 antennas.

SR = 15.36e6;
T = SR * 1e-3;
cdl.SampleRate = SR;
cdlinfo = info(cdl);
Nt = cdlinfo.NumTransmitAntennas;
Nsub = 2048;                                     %No of subcarriers
txWaveform = complex(randn(T,Nt),randn(T,Nt));

 
%% Generating the wireless channel

[rxWaveform, pathgains] = cdl(txWaveform);
N_train = size(pathgains,1);
temp = zeros(N_train,Nsub,NTx,NRx); Hf = zeros(size(temp));
for j=1:N_train
temp(j,:,:,:) = cat(1,squeeze(pathgains(j,:,:,:)),zeros(Nsub-size(pathgains,2),size(pathgains,3),size(pathgains,4)));
Hf(j,:,:,:) = fftshift(fft(squeeze(temp(j,:,:,:))));
end
H1 = squeeze(Hf(:,1,:,:));                      %Channel for one subcarrier
%% Generating uplink and downlink channels

Hup = zeros(size(squeeze(Hf(:,1,:,:)))); Hdl = zeros(size(Hup));
% Generating RF chain gains
tb1 = diag(randn(NTx,1) + 1j*randn(NTx,1));
rb1 = diag(randn(NTx,1) + 1j*randn(NTx,1));
tu1 = diag(randn(NRx,1) + 1j*randn(NRx,1));
ru1 = diag(randn(NRx,1) + 1j*randn(NRx,1));

for j=1:N_train
tb = tb1 + eye(size(tb1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
rb = rb1 + eye(size(rb1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
tu = tu1 + eye(size(tu1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
ru = ru1 + eye(size(ru1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
Hup(j,:,:) = rb*squeeze(H1(j,:,:))*tu;
Hdl(j,:,:) = ru*conj(squeeze(H1(j,:,:)))'*tb;
end

%% Building the training data for 1 subcarrier

H = reshape(Hup,N_train,NTx*NRx);
%S = svd(H);
H_real = real(H); H_imag = imag(H);
H_final = cat(2,H_real,H_imag);
%% PCA
 % Run PCA
[U, S] = pca(H_final);
H_compressed = H_final*U(:,1:10);


%% Save matrices
 save('H_uplink.mat',Hup);
 save('H_downlink.mat',Hdl);
 save('Training_data.mat',H_final);
% % %% SVD compression
% % H1 = squeeze(Hf(1,1,:,:));
% % [U1,S1,Vtemp] = svd(H1); V1 = (Vtemp)';
% % check = diag(S1) > max(diag(S1))/1000;
% % r1 = nnz(check); 
% % %r1=16;
% % H2 = U1(:,1:r1)*S1(1:r1,1:r1)*V1(1:r1,:);
% % svd_loss = norm(abs(H2-H1));