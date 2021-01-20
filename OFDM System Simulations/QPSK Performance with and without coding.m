clc, clear;
% Generate Random Bit stream bk
L = 1e6;
Eb_QPSK = 1;
bk = randi([0, 1], 1, L);

% Generate QPSK symbols based on the bit stream bk_QPSK
bk_QPSK = zeros(1, L / 2);
for i= 1:2:L
    if (bk(i) == 0 && bk(i + 1) == 0)
        bk_QPSK(floor(i / 2) + 1) = -1 - 1i;
    elseif (bk(i) == 0 && bk(i + 1) == 1)
        bk_QPSK(floor(i / 2) + 1) = -1 + 1i;
    elseif (bk(i) == 1 && bk(i + 1) == 1)
        bk_QPSK(floor(i / 2) + 1) = 1 + 1i;
    elseif (bk(i) == 1 && bk(i + 1) == 0)
        bk_QPSK(floor(i / 2) + 1) = 1 - 1i;
    end
end

% scatterplot(bk_QPSK)

% Generate the complex channel vector and complex noise vector
N = 1:10;
N0 = Eb_QPSK ./ (10 .^ (N ./ 10));      % assumption
BER_v = zeros(1, 10);
SNR = zeros(1, 10);
for i = 1:10
    hr = normrnd(0, sqrt(0.5), 1, L / 2);
    hi = normrnd(0, sqrt(0.5), 1, L / 2);
    channel_QPSK = hr + 1j * hi;
    nc = normrnd(0, sqrt(N0(i) / 2), 1, L / 2);
    ns = normrnd(0, sqrt(N0(i) / 2), 1, L / 2);
    Noise_QPSK = nc + 1j * ns;
    % Compute the received symbol vector yk
    yk = channel_QPSK .* bk_QPSK + Noise_QPSK;
    % Compensate for the channel gain at the receiver
    bk_telda = yk ./ channel_QPSK;                
    
    % make rough estimation
    bk_telda_QPSK = zeros(1, L);
    j = 1;
    for k = 1:length(bk_telda)
        if (real(bk_telda(k)) > 0 && imag(bk_telda(k)) > 0)
            bk_telda_QPSK(j)        = 1;
            bk_telda_QPSK(j + 1)    = 1;
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) > 0)
            bk_telda_QPSK(j)        = 0;
            bk_telda_QPSK(j + 1)    = 1;
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) < 0)
            bk_telda_QPSK(j)        = 0;
            bk_telda_QPSK(j + 1)    = 0;
        elseif (real(bk_telda(k)) > 0 && imag(bk_telda(k)) < 0)
            bk_telda_QPSK(j)        = 1;
            bk_telda_QPSK(j + 1)    = 0;
        end
        j = j + 2;
    end
    % Compute the bit-error rate (BER) for SNR
    BER = 0;
    for n = 1:L
        if(bk_telda_QPSK(n) ~= bk(n))
            BER = BER + 1;
        else
            continue;
        end
    end
    BER_v(i) = BER / L;
    SNR(i) = 10 * log10(Eb_QPSK / N0(i));
    %     fprintf("BER = %f\n", BER);
    %     fprintf("SNR = %f\n", SNR);
end
figure
semilogy(SNR, BER_v);
title('BER vs. SNR for QPSK');
xlabel('SNR (dB)');
ylabel('BER');
grid on;
hold on;

% % Theoritacl Calculation
% BER_TH = zeros(length(N), 20);
% BER_TH(:, 1) = berfading(N, 'psk', 4, 1);
% semilogy(N, BER_TH)
% legend('Practical BER','Theoritical BER')

%% Repeat the code
clc, clear
% Generate Random Bit stream bk
L = 1e6;
Eb_QPSK = 1;
bk = randi([0 1], 1, L);

% Generate QPSK symbols based on the bit stream xk
xk = bk;
xk_rep = repelem(xk, 3);
% Generate QPSK symbols based on the bit stream bk_QPSK
bk_QPSK = zeros(1, 3 * (L / 2));
for i= 1:2:3*L
    if (xk_rep(i) == 0 && xk_rep(i + 1) == 0)
        bk_QPSK(floor(i / 2) + 1) = -1 - 1i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 1)
        bk_QPSK(floor(i / 2) + 1) = -1 + 1i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 1)
        bk_QPSK(floor(i / 2) + 1) = 1 + 1i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 0)
        bk_QPSK(floor(i / 2) + 1) = 1 - 1i;
    end
end

xk_ver = zeros(1, L);
BER_ver = zeros(1, 10);
% Generate the complex channel vector and complex noise vector
N = 1:10;
N0 = Eb_QPSK ./ (10 .^ (N ./ 10)); % assumption
SNR = zeros(1, 10);
for i = 1:10
    hr = normrnd(0, sqrt(0.5), 1, 3 * L / 2);
    hi = normrnd(0, sqrt(0.5), 1, 3 * L / 2);
    channel_QPSK = hr + 1j * hi;
    nc = normrnd(0, sqrt(N0(i) / 2), 1, 3 * L / 2);
    ns = normrnd(0, sqrt(N0(i) / 2), 1, 3 * L / 2);
    Noise_QPSK = nc + 1j * ns;
    % Compute the received symbol vector yk
    yk = channel_QPSK .* bk_QPSK + Noise_QPSK;
    % Compensate for the channel gain at the receiver
    bk_telda = yk ./ channel_QPSK;                
    j = 1;
    for k = 1:length(bk_telda)
        if (real(bk_telda(k)) > 0 && imag(bk_telda(k)) > 0)
            bk_telda_QPSK(j)        = 1;
            bk_telda_QPSK(j + 1)    = 1;
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) > 0)
            bk_telda_QPSK(j)        = 0;
            bk_telda_QPSK(j + 1)    = 1;
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) < 0)
            bk_telda_QPSK(j)        = 0;
            bk_telda_QPSK(j + 1)    = 0;
        elseif (real(bk_telda(k)) > 0 && imag(bk_telda(k)) < 0)
            bk_telda_QPSK(j)        = 1;
            bk_telda_QPSK(j + 1)    = 0;
        end
        j = j + 2;
    end
   
    % make the decsion
    counter = 1;
    for j = 1:3:3*L
        if((bk_telda_QPSK(j) + bk_telda_QPSK(j + 1) + bk_telda_QPSK(j + 2)) > 1)
            xk_ver(counter) = 1;
        else
            xk_ver(counter) = 0;
        end
        counter = counter + 1;
    end
    
    % Compute the bit-error rate (BER) for SNR
    BER = 0;
    for n = 1:L
        if(xk_ver(n) ~= xk(n))
            BER = BER + 1;
        else
            continue;
        end
    end
    BER_ver(i) = BER / L;
    SNR(i) = 10 * log10(Eb_QPSK / N0(i));
end
semilogy(SNR, BER_ver);
title('BER vs. SNR for QPSK after repetition');
xlabel('SNR (dB)');
ylabel('BER');
grid on;

%% This section For Graphs
title('BER vs. SNR for QPSK before and after repetaion');
legend('BPSK before repetition', 'QPSK after repetition');
