clc; close all;
NTx = 32; NRx = 4; N_train = 100;  %NTx - No of BTS antenna, NRx - No of UE antenna

%% Generating correlated wireless channel matrix
H = zeros(NTx,NRx); Hup = zeros(N_train,NTx,NRx); Hdl = zeros(N_train,NRx,NTx); a = randi([0 1],1,NTx-2);
for j = 1:NTx
    for k = 1:NRx
 
    H(j,k) = (1/sqrt(2)).*(randn + 1j*randn);   %Generate normally distributed random numbers
    end
    %H(i,j,:) = H(i,j,:)./(norm(max(H(i,j,:))));   % Normalise them such that maximum norm is 1
end
for i=1:N_train
%Generating RF chain gains
tb = diag(randn(NTx,1) + 1j*randn(NTx,1));
rb = diag(randn(NTx,1) + 1j*randn(NTx,1));
tu = diag(randn(NRx,1) + 1j*randn(NRx,1));
ru = diag(randn(NRx,1) + 1j*randn(NRx,1));
% tb = tb1 + eye(size(tb1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
% rb = rb1 + eye(size(rb1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
% tu = tu1 + eye(size(tu1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
% ru = ru1 + eye(size(ru1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
Hup(i,:,:) = rb*H*tu;
Hdl(i,:,:) = ru*transpose(H)*tb;
end

Hup1 = reshape(Hup,N_train,NTx*NRx);
Hdl1 = reshape(Hdl,N_train,NTx*NRx);

%% Zero mean


[Hup_norm, mu] = featureNormalize(Hup1);
[Hdl_norm, mu1] = featureNormalize(Hdl1);


%% Uplink and downlink correlation matrix

Upcorr = (1/N_train)*(Hup_norm*(Hup_norm'));
Dlcorr = (1/N_train)*(Hdl_norm*(Hdl_norm'));
diff = Upcorr - Dlcorr;


%% Uplink Channel estimation
%For the uplink channel estimate, no of transmit antennas is NRx and no of
%receive antennas in NTx
Nsamples = NRx;
% Generating pilot waveform
X = zeros(N_train,NRx,Nsamples);
for i=1:N_train
X(i,:,:) = eye(NRx,Nsamples).*(rand(1,1)+1j*rand(1,1)); 
end

Y = zeros(N_train,NTx,Nsamples); Hest = zeros(N_train,NTx,NRx);
error = zeros(size(Hup)); e = zeros(N_train,1); 
SNR = 30;              

for i=1:N_train
    Y(i,:,:) = squeeze(Hup(i,:,:))*squeeze(X(i,:,:)); % Y = HX
    for j=1:NRx                                       %Adding noise
    a = squeeze(Y(i,j,:));
    p = mean(power(abs(squeeze(Y(i,j,:))),2));        %Power of input signal
    Y(i,j,:) = awgn(squeeze(Y(i,j,:)),SNR,(1/length(a))*(a'*a));
    end
    Hest(i,:,:) = squeeze(Y(i,:,:))*pinv(squeeze(X(i,:,:)));
    error(i,:,:) = power(abs(squeeze(Hest(i,:,:))-squeeze(Hup(i,:,:))),2);
    e(i) = norm(squeeze(error(i,:,:)),'fro');         %Frobenius norm of error
end