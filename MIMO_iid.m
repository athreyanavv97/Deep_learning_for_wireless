clc; close all;
NTx = 32; NRx = 4; N_train = 100;  %NTx - No of BTS antenna, NRx - No of UE antenna
%% Generating RF chain gains
tb1 = diag(randn(NTx,1) + 1j*randn(NTx,1));
rb1 = diag(randn(NTx,1) + 1j*randn(NTx,1));
tu1 = diag(randn(NRx,1) + 1j*randn(NRx,1));
ru1 = diag(randn(NRx,1) + 1j*randn(NRx,1));

%% Generating correlated wireless channel matrix
H = zeros(N_train,NTx,NRx); Hup = zeros(N_train,NTx,NRx); Hdl = zeros(N_train,NRx,NTx); a = randi([0 1],1,NTx-2);
for i=1:N_train
tb = tb1 + eye(size(tb1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
rb = rb1 + eye(size(rb1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
tu = tu1 + eye(size(tu1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
ru = ru1 + eye(size(ru1)).*(0.01*(rand(1,1)+1j*rand(1,1)));
for j = 1:NTx
    for k = 1:NRx
 
    H(i,j,k) = (1/sqrt(2)).*(randn + 1j*randn);   %Generate normally distributed random numbers
    end
    %H(i,j,:) = H(i,j,:)./(norm(max(H(i,j,:))));   % Normalise them such that maximum norm is 1
end
Hup(i,:,:) = rb*squeeze(H(i,:,:))*tu;
Hdl(i,:,:) = ru*conj(squeeze(H(i,:,:)))'*tb;
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

X = eye(NRx,NRx).*(rand(1,1)+1j*rand(1,1)); Y = zeros(N_train,NTx,NRx); Hest = zeros(N_train,NTx,NRx);
error = zeros(size(Hup));
SNR = 100;              
for i=1:N_train
    Y(i,:,:) = squeeze(Hup(i,:,:))*X;
    for j=1:NRx                                     %Adding noise
    a = squeeze(Y(i,j,:));
    p = mean(power(abs(squeeze(Y(i,j,:))),2));      %Power of input signal
    Y(i,j,:) = awgn(Y(i,j,:),SNR,(1/length(a))*(a'*a));
    end
    Hest(i,:,:) = squeeze(Y(i,:,:))*pinv(X);
    error(i,:,:) = power(abs(squeeze(Hest(i,:,:))-squeeze(Hup(i,:,:))),2);
end