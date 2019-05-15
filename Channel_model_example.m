%% Transmission over Channel Model with Delay Profile CDL-D
% Transmit waveform through a Clustered Delay Line (CDL) channel model with 
% delay profile CDL-D from TR 38.901 Section 7.7.1.
% 
% Define the channel configuration structure using an |nrCDLChannel| System 
% object. Use delay profile CDL-D, a delay spread of 10 ns, and UT velocity of 
% 15 km/h:

clc; close all;
v = 180.0;                    % UT velocity in km/h
fc = 3.5e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd = (v*1000/3600)/c*fc;     % UT max Doppler frequency in Hz
cdl = nrCDLChannel;
cdl.SampleDensity = 64;
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
cdl.ReceiveAntennaArray.Size = [sqrt(NRx) sqrt(NRx) 1 1 1];
cdl.ReceiveAntennaArray.ElementSpacing = [0.5 0.5 1 1];

%% 
% Create a random waveform of 1 subframe duration with 8 antennas.

SR = 15.36e6;
T = SR * 1e-3;
cdl.SampleRate = SR;
cdlinfo = info(cdl);
Nt = cdlinfo.NumTransmitAntennas;
Nsub = 2048;                                     %No of subcarriers
txWaveform = complex(randn(T,Nt),randn(T,Nt));

%% 
% Transmit the input waveform through the channel.

[rxWaveform, pathgains] = cdl(txWaveform);
N_train = size(pathgains,1);
temp = zeros(N_train,Nsub,NTx,NRx); Hf = zeros(size(temp));
for j=1:N_train
temp(j,:,:,:) = cat(1,squeeze(pathgains(j,:,:,:)),zeros(Nsub-size(pathgains,2),size(pathgains,3),size(pathgains,4)));
Hf(j,:,:,:) = fft(squeeze(temp(j,:,:,:)));
end

%% Building the training data for 1 subcarrier

H = reshape(squeeze(Hf(:,1,:,:)),N_train,NTx*NRx);
S = svd(H);
 
