clc, clear;

% Generate Random Bit stream bk
L = 1e6;
Eb_16QAM = 2.5;
bk = randi([0, 1], 1, L);
bk_16QAM = zeros(1, L / 4);


% Generate 16QAM symbols based on the bit stream bk_16QAM
k = 1;
for i = 1 : 4 : L
    if (bk(i) == 0 && bk(i + 1) == 0 && bk(i + 2) == 0 && bk(i + 3) == 0)       %0000
        bk_16QAM(k) = -3 - 3i;
    elseif (bk(i) == 0 && bk(i + 1) == 0 && bk(i + 2) == 0 && bk(i + 3) == 1)   %0001
        bk_16QAM(k) = -3 - 1i;
    elseif (bk(i) == 0 && bk(i + 1) == 0 && bk(i + 2) == 1 && bk(i + 3) == 0)   %0010
        bk_16QAM(k) = -3 + 3i;
    elseif (bk(i) == 0 && bk(i + 1) == 0 && bk(i + 2) == 1 && bk(i + 3) == 1)   %0011
        bk_16QAM(k) = -3 + 1i;
    elseif (bk(i) == 0 && bk(i + 1) == 1 && bk(i + 2) == 0 && bk(i + 3) == 0)   %0100
        bk_16QAM(k) = -1 - 3i;
    elseif (bk(i) == 0 && bk(i + 1) == 1 && bk(i + 2) == 0 && bk(i + 3) == 1)   %0101
        bk_16QAM(k) = -1 - 1i;
    elseif (bk(i) == 0 && bk(i + 1) == 1 && bk(i + 2) == 1 && bk(i + 3) == 0)   %0110
        bk_16QAM(k) = -1 + 3i;
    elseif (bk(i) == 0 && bk(i + 1) == 1 && bk(i + 2) == 1 && bk(i + 3) == 1)   %0111
        bk_16QAM(k) = -1 + 1i;
    elseif (bk(i) == 1 && bk(i + 1) == 0 && bk(i + 2) == 0 && bk(i + 3) == 0)   %1000
        bk_16QAM(k) = 3 - 3i;
    elseif (bk(i) == 1 && bk(i + 1) == 0 && bk(i + 2) == 0 && bk(i + 3) == 1)   %1001
        bk_16QAM(k) = 3 - 1i;
    elseif (bk(i) == 1 && bk(i + 1) == 0 && bk(i + 2) == 1 && bk(i + 3) == 0)   %1010
        bk_16QAM(k) = 3 + 3i;
    elseif (bk(i) == 1 && bk(i + 1) == 0 && bk(i + 2) == 1 && bk(i + 3) == 1)   %1011
        bk_16QAM(k) = 3 + 1i;
    elseif (bk(i) == 1 && bk(i + 1) == 1 && bk(i + 2) == 0 && bk(i + 3) == 0)   %1100
        bk_16QAM(k) = 1 - 3i;
    elseif (bk(i) == 1 && bk(i + 1) == 1 && bk(i + 2) == 0 && bk(i + 3) == 1)   %1101
        bk_16QAM(k) = 1 - 1i;
    elseif (bk(i) == 1 && bk(i + 1) == 1 && bk(i + 2) == 1 && bk(i + 3) == 0)   %1110
        bk_16QAM(k) = 1 + 3i;
    elseif (bk(i) == 1 && bk(i + 1) == 1 && bk(i + 2) == 1 && bk(i + 3) == 1) 	%1111
        bk_16QAM(k) = 1 + 1i;
    end
    k = k + 1;
end
% scatterplot(bk_16QAM)

% Generate the complex channel vector and complex noise vector
N = 1:10;
N0 = Eb_16QAM ./ (10 .^ (N ./ 10));      % assumption
BER_v = zeros(1, 10);
SNR = zeros(1, 10);
for i = 1:10
    hr = normrnd(0, sqrt(0.5), 1, L / 4);
    hi = normrnd(0, sqrt(0.5), 1, L / 4);
    channel_16QAM = hr + 1j * hi;
    nc = normrnd(0, sqrt(N0(i) / 2), 1, L / 4);
    ns = normrnd(0, sqrt(N0(i) / 2), 1, L / 4);
    Noise_16QAM = nc + 1j * ns;
    % Compute the received symbol vector yk
    yk = channel_16QAM .* bk_16QAM + Noise_16QAM;
    % Compensate for the channel gain at the receiver
    bk_telda = yk ./ channel_16QAM;                % D7kware
    
    % make rough estimation
    bk_telda_16QAM = zeros(1, L);
    j = 1;
    for k = 1:length(bk_telda)
        if (real(bk_telda(k)) > 0 && imag(bk_telda(k)) > 0)     % first quarter
            if (real(bk_telda(k)) > 2 && imag(bk_telda(k)) > 2) % 1010
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > 2 && imag(bk_telda(k)) < 2) % 1011
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) > 2) % 1110
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) < 2) % 1111
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            end
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) > 0)     % second quarter
            if (real(bk_telda(k)) < -2 && imag(bk_telda(k)) > 2)    % 0010
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < -2 && imag(bk_telda(k)) < 2) % 0011
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) > 2) % 0110
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) < 2) % 0111
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            end
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) < 0)     % third quarter
            if (real(bk_telda(k)) < -2 && imag(bk_telda(k)) < -2)    % 0000
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < -2 && imag(bk_telda(k)) > -2) % 0001
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) < -2) % 0100
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) > -2) % 0101
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            end
        elseif (real(bk_telda(k)) > 0 && imag(bk_telda(k)) < 0)     % fourth quarter
            if (real(bk_telda(k)) > 2 && imag(bk_telda(k)) < -2)    % 1000
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > 2 && imag(bk_telda(k)) > -2) % 1001
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) < -2) % 1100
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) > -2) % 1101
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            end
        end
        j = j + 4;
    end
    
    % Compute the bit-error rate (BER) for SNR
    BER = 0;
%     for n = 1:L
%         if(bk_telda_16QAM(n) ~= bk(n))
%             BER = BER + 1;
%         else
%             continue;
%         end
%     end
    BER = biterr(bk_telda_16QAM, bk);
    BER_v(i) = BER / L;
    SNR(i) = 10 * log10(Eb_16QAM / N0(i));
    %     fprintf("BER = %f\n", BER);
    %     fprintf("SNR = %f\n", SNR);
end
semilogy(SNR, BER_v);
title('BER vs. SNR for 16QAM');
xlabel('SNR (dB)');
ylabel('BER');
grid on;
hold on;

% % Theoritical calculation
% N = 1:10;
% BER_TH = zeros(length(N), 20);
% BER_TH(:, 1) = berfading(N, 'qam', 16, 1);
% semilogy(N, BER_TH)
% legend('Practical BER', 'Theoritical BER')

%% Repeat the code
clc, clear
% Generate Random Bit stream bk
L = 1e6;
Eb_16QAM = 1;
bk = randi([0 1], 1, L);

% Generate 16QAM symbols based on the bit stream xk
xk = bk;
xk_rep = repelem(xk, 3);

bk_16QAM = zeros(1, 3 * L / 4);


% Generate 16QAM symbols based on the bit stream bk_16QAM
k = 1;
for i = 1 : 4 : 3 * L
    if (xk_rep(i) == 0 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 0)       %0000
        bk_16QAM(k) = -3 - 3i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 1)   %0001
        bk_16QAM(k) = -3 - 1i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 0)   %0010
        bk_16QAM(k) = -3 + 3i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 1)   %0011
        bk_16QAM(k) = -3 + 1i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 0)   %0100
        bk_16QAM(k) = -1 - 3i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 1)   %0101
        bk_16QAM(k) = -1 - 1i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 0)   %0110
        bk_16QAM(k) = -1 + 3i;
    elseif (xk_rep(i) == 0 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 1)   %0111
        bk_16QAM(k) = -1 + 1i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 0)   %1000
        bk_16QAM(k) = 3 - 3i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 1)   %1001
        bk_16QAM(k) = 3 - 1i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 0)   %1010
        bk_16QAM(k) = 3 + 3i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 0 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 1)   %1011
        bk_16QAM(k) = 3 + 1i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 0)   %1100
        bk_16QAM(k) = 1 - 3i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 0 && xk_rep(i + 3) == 1)   %1101
        bk_16QAM(k) = 1 - 1i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 0)   %1110
        bk_16QAM(k) = 1 + 3i;
    elseif (xk_rep(i) == 1 && xk_rep(i + 1) == 1 && xk_rep(i + 2) == 1 && xk_rep(i + 3) == 1) 	%1111
        bk_16QAM(k) = 1 + 1i;
    end
    k = k + 1;
end

xk_ver = zeros(1, L);
BER_ver = zeros(1, 10);
% Generate the complex channel vector and complex noise vector
N = 1:10;
N0 = Eb_16QAM ./ (10 .^ (N ./ 10)); % assumption
SNR = zeros(1, 10);
for i = 1:10
    hr = normrnd(0, sqrt(0.5), 1, 3 * L / 4);
    hi = normrnd(0, sqrt(0.5), 1, 3 * L / 4);
    channel_16QAM = hr + 1j * hi;
    nc = normrnd(0, sqrt(N0(i) / 2), 1, 3 * L / 4);
    ns = normrnd(0, sqrt(N0(i) / 2), 1, 3 * L / 4);
    Noise_16QAM = nc + 1j * ns;
    % Compute the received symbol vector yk
    yk = channel_16QAM .* bk_16QAM + Noise_16QAM;
    % Compensate for the channel gain at the receiver
    bk_telda = yk ./ channel_16QAM;                
    % make rough estimation
    bk_telda_16QAM = zeros(1, L);
    j = 1;
    for k = 1:length(bk_telda)
        if (real(bk_telda(k)) > 0 && imag(bk_telda(k)) > 0)     % first quarter
            if (real(bk_telda(k)) > 2 && imag(bk_telda(k)) > 2) % 1010
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > 2 && imag(bk_telda(k)) < 2) % 1011
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) > 2) % 1110
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) < 2) % 1111
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            end
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) > 0)     % second quarter
            if (real(bk_telda(k)) < -2 && imag(bk_telda(k)) > 2)    % 0010
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < -2 && imag(bk_telda(k)) < 2) % 0011
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) > 2) % 0110
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) < 2) % 0111
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 1;
                bk_telda_16QAM(j + 3)    = 1;
            end
        elseif (real(bk_telda(k)) < 0 && imag(bk_telda(k)) < 0)     % third quarter
            if (real(bk_telda(k)) < -2 && imag(bk_telda(k)) < -2)    % 0000
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < -2 && imag(bk_telda(k)) > -2) % 0001
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) < -2) % 0100
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > -2 && imag(bk_telda(k)) > -2) % 0101
                bk_telda_16QAM(j)        = 0;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            end
        elseif (real(bk_telda(k)) > 0 && imag(bk_telda(k)) < 0)     % fourth quarter
            if (real(bk_telda(k)) > 2 && imag(bk_telda(k)) < -2)    % 1000
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) > 2 && imag(bk_telda(k)) > -2) % 1001
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 0;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) < -2) % 1100
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 0;
            elseif (real(bk_telda(k)) < 2 && imag(bk_telda(k)) > -2) % 1101
                bk_telda_16QAM(j)        = 1;
                bk_telda_16QAM(j + 1)    = 1;
                bk_telda_16QAM(j + 2)    = 0;
                bk_telda_16QAM(j + 3)    = 1;
            end
        end
        j = j + 4;
    end
    
    % make the decsion
    counter = 1;
    for j = 1:3:3*L
        if((bk_telda_16QAM(j) + bk_telda_16QAM(j + 1) + bk_telda_16QAM(j + 2)) > 1)
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
    SNR(i) = 10 * log10(Eb_16QAM / N0(i));
end
semilogy(SNR, BER_ver);
title('BER vs. SNR for 16QAM after repetition');
xlabel('SNR (dB)');
ylabel('BER');
grid on;

%% This section For Graphs
title('BER vs. SNR for 16QAM before and after repetaion');
legend('16QAM before repetition', '16QAM after repetition');
