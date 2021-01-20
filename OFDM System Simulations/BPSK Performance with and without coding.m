clc, clear, close all
% Generate Random Bit stream bk
L = 1e6;
Eb_BPSK = 1;
bk = randi([0 1], 1, L);

% Generate BPSK symbols based on the bit stream xk
xk = bk;
xk(xk == 1) = sqrt(Eb_BPSK);
xk(xk == 0) = -sqrt(Eb_BPSK);

% Generate the complex channel vector and complex noise vector
N = 1:10;
N0 = Eb_BPSK ./ (10 .^ (N ./ 10)); % assumption
BER_v = zeros(1, 10);
SNR = zeros(1, 10);
for i = 1:10
    hr = normrnd(0, sqrt(0.5), 1, L);
    hi = normrnd(0, sqrt(0.5), 1, L);
    channel_BPSK = hr + 1j * hi;
    nc = normrnd(0, sqrt(N0(i) / 2), 1, L);
    ns = normrnd(0, sqrt(N0(i) / 2), 1, L);
    Noise_BPSK = nc + 1j * ns;
    % Compute the received symbol vector yk
    yk = channel_BPSK .* xk + Noise_BPSK;
    % Compensate for the channel gain at the receiver
    bk_telda = yk ./ channel_BPSK;                
    bk_telda(bk_telda > 0) = 1;
    bk_telda(bk_telda < 0) = 0;
    
    % Compute the bit-error rate (BER) for SNR
    BER = 0;
    for n = 1:L
        if(bk_telda(n) ~= bk(n))
            BER = BER + 1;
        else
            continue;
        end
    end
    BER_v(i) = BER / L;
    SNR(i) = 10 * log10(Eb_BPSK / N0(i));
end
semilogy(SNR, BER_v);
title('BER vs. SNR for BPSK before repetition');
xlabel('SNR (dB)');
ylabel('BER');
grid on;
hold on;

% % Theoritacl Calculation
% BER_TH = zeros(length(N), 20);
% BER_TH(:, 1) = berfading(N, 'psk', 4, 1);
% semilogy(N, BER_TH)
% legend('Practical BER','Theoritical BER')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Repeat the code
clc, clear
% Generate Random Bit stream bk
L = 1e6;
Eb_BPSK = 1;
bk = randi([0 1], 1, L);

% Generate BPSK symbols based on the bit stream xk
xk = bk;
xk(xk == 1) = sqrt(Eb_BPSK);
xk(xk == 0) = -sqrt(Eb_BPSK);
xk_rep = repelem(xk, 3);
xk_ver = zeros(1, L);
BER_ver = zeros(1, 10);
% Generate the complex channel vector and complex noise vector
N = 1:10;
N0 = Eb_BPSK ./ (10 .^ (N ./ 10)); % assumption
SNR = zeros(1, 10);
for i = 1:10
    hr = normrnd(0, sqrt(0.5), 1, 3 * L);
    hi = normrnd(0, sqrt(0.5), 1, 3 * L);
    channel_BPSK = hr + 1j * hi;
    nc = normrnd(0, sqrt(N0(i) / 2), 1, 3 * L);
    ns = normrnd(0, sqrt(N0(i) / 2), 1, 3 * L);
    Noise_BPSK = nc + 1j * ns;
    % Compute the received symbol vector yk
    yk = channel_BPSK .* xk_rep + Noise_BPSK;
    % Compensate for the channel gain at the receiver
    bk_telda = yk ./ channel_BPSK;                
    bk_telda(bk_telda > 0) = 1;
    bk_telda(bk_telda < 0) = 0;
    
    % make the decsion
    counter = 1;
    for j = 1:3:3*L
        if((bk_telda(j) + bk_telda(j+1) + bk_telda(j+2)) > 1)
            xk_ver(counter) = 1;
        else
            xk_ver(counter) = 0;
        end
        counter = counter + 1;
    end
    
    % Compute the bit-error rate (BER) for SNR
    BER = 0;
    for n = 1:L
        if(xk_ver(n) ~= bk(n))
            BER = BER + 1;
        else
            continue;
        end
    end
    BER_ver(i) = BER / L;
    SNR(i) = 10 * log10(Eb_BPSK / N0(i));
    %     fprintf("BER = %f\n", BER);
    %     fprintf("SNR = %f\n", SNR);
end

semilogy(SNR, BER_ver);
title('BER vs. SNR for BPSK after repetition');
xlabel('SNR (dB)');
ylabel('BER');
grid on;

%% This section For Graphs
title('BER vs. SNR for BPSK before and after repetaion');
legend('BPSK before repition', 'BPSK after repetition')
