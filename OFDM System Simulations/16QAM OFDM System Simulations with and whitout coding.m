clc, clear, close all;

% Parameters
Eb_16QAM = 1;
CYC_PERFIX = 8;
L = 1e3;                    % Number of symbols (1000)
N = L * (32 + CYC_PERFIX);  % Number of bits (40000)
SNR = 1:10;                 % Assumption SNR in (dB)
OFDM_SYM = [];              % the OFDM symblos
bk_16QAM_ifft_CYC = [] ;     % the output from the cyclic extension
bk_receivied = [];          % the recevied data at the receiver
BER_16QAM = [];

% Ask for repetition (NO == 1) or (YES == 2) and check for validity
repetition = menu('Do you need repetition', 'NO', 'YES');
if repetition == 1          % NO repetiton
    fprintf('No Repetion will happen\n');
elseif repetition == 2      % YES repetiotn
    repetition = 3;         % set it to 3 to use it later for the number of repetiton
    fprintf('1/3 code repetition will happen\n');
else                        % wrong choise
    fprintf('Please choose YES or NO\n');
end

for x1 = 1:L
    % Generate the random bits
    if repetition == 3      % YES repetition will happen
        NO_OF_BITS = 42;
    else                    % NO repetition will happen
        NO_OF_BITS = 128;
    end
    bk = randi([0, 1], 1, NO_OF_BITS);
    
    % Coding the data if needed (This part is for the repetition)
    bk_rep = repelem(bk, repetition);
    if (repetition == 3)
        bk_rep = [bk_rep, zeros(1,2)];          % Concatinate a zero in the end
    end
    
    % construct the OFDM symbols
    OFDM_SYM = [OFDM_SYM, bk];
    
    % passing the data througth the interleaver block
    bk_intrlv = matintrlv(bk_rep, 8, 16);
    
    % mapping the data into 16QAM
    bk_16QAM = [];
    x2 = 1;
    for x3 = 1 : 4 : 64
        if (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 0)       %0000
            bk_16QAM(x2) = -3 - 3i;
        elseif (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 1)   %0001
            bk_16QAM(x2) = -3 - 1i;
        elseif (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 0)   %0010
            bk_16QAM(x2) = -3 + 3i;
        elseif (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 1)   %0011
            bk_16QAM(x2) = -3 + 1i;
        elseif (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 0)   %0100
            bk_16QAM(x2) = -1 - 3i;
        elseif (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 1)   %0101
            bk_16QAM(x2) = -1 - 1i;
        elseif (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 0)   %0110
            bk_16QAM(x2) = -1 + 3i;
        elseif (bk_intrlv(x3) == 0 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 1)   %0111
            bk_16QAM(x2) = -1 + 1i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 0)   %1000
            bk_16QAM(x2) = 3 - 3i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 1)   %1001
            bk_16QAM(x2) = 3 - 1i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 0)   %1010
            bk_16QAM(x2) = 3 + 3i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 0 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 1)   %1011
            bk_16QAM(x2) = 3 + 1i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 0)   %1100
            bk_16QAM(x2) = 1 - 3i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 0 && bk_intrlv(x3 + 3) == 1)   %1101
            bk_16QAM(x2) = 1 - 1i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 0)   %1110
            bk_16QAM(x2) = 1 + 3i;
        elseif (bk_intrlv(x3) == 1 && bk_intrlv(x3 + 1) == 1 && bk_intrlv(x3 + 2) == 1 && bk_intrlv(x3 + 3) == 1) 	%1111
            bk_16QAM(x2) = 1 + 1i;
        end
        x2 = x2 + 1;
    end
    
    % Passing the data tjrougth the ifft block
    bk_16QAM_ifft = ifft(bk_16QAM, 32);
    
    % Add the cyclic extension
    bk_16QAM_ifft_CYC = [bk_16QAM_ifft_CYC, bk_16QAM_ifft(32 - CYC_PERFIX + 1:32) bk_16QAM_ifft];
end


% Ask for the channel type (Flat Fadding == 1) or (Frequncy Selective == 2)
CH_TP = menu('Choose the channle type', 'Flat Fadding', 'Frequency Selective');

% Passing the data througth the channel, then design the receiver
% Parameters

for x4 = SNR
    SNR_LIN = 10 ^ (x4 / 10);
    N0 = Eb_16QAM / SNR_LIN;
    
    % Model the noise
    nc = normrnd(0, sqrt(N0 / 2), 1, N);
    ns = normrnd(0, sqrt(N0 / 2), 1, N);
    Noise_16QAM = nc + 1j * ns;
    
    
    if CH_TP == 1               % Flat Fadding channel
        %---------------------------------------------------------------%
        % CHANNEL_RAYLEIGH_FLAT_FADING_CHANNEL
        hr = normrnd(0, sqrt(0.5), 1, L);
        hi = normrnd(0, sqrt(0.5), 1, L);
        channel_16QAM = hr + 1j * hi;
        
        % Passing the data througth the channel
        yk = [];                % the output of the data after the channel and noise
        yk_telda = [];          % the ouput after the channel compensator
        x5 = 1;
        for x6 = 1 : (CYC_PERFIX + 32) : L * (CYC_PERFIX + 32)
            yk = [yk, bk_16QAM_ifft_CYC(x6:x6 + 31 + CYC_PERFIX) * channel_16QAM(x5)];
            x5 = x5 + 1;
        end
        
        % Adding the noise
        yk_noise = yk + Noise_16QAM;
        
        % Design the receiver
        % First remove the channel effect
        x7 = 1;
        for x8 = 1 : (CYC_PERFIX + 32) : L * (CYC_PERFIX + 32)
            yk_telda = [yk_telda yk_noise(x8:x8 + CYC_PERFIX +31) / channel_16QAM(x7)];
            x7 = x7 + 1;
        end
        %----------------------------------------------------------------------%
    elseif CH_TP == 2                   % Frequncy Selective cghannel
        %----------------------------------------------------------------------%
        %FREQUENCY_SELECTIVE_FADING_CHANNEL
        hr = transpose(normrnd(0, sqrt(0.5), 8, L));
        hi = transpose(normrnd(0, sqrt(0.5), 8, L));
        channel_16QAM = hr + 1i * hi;
        
        % Passing the data througth the channel
        yk = [];
        yk_telda_intermediate = [];
        x14 = 1;
        bk_channel_sel = bk_16QAM_ifft_CYC + Noise_16QAM;
        for x6 = 1 : (CYC_PERFIX + 32) : L * (CYC_PERFIX + 32)
            temp1 = conv(bk_channel_sel(x6:x6 + CYC_PERFIX +31), channel_16QAM(x14, :));
            yk = [yk, temp1];
            x14 = x14 + 1;
        end
        yk_copy = yk;
        
        x15 = 1;
        for x8 = 1 : (CYC_PERFIX + 32 + 7) : L * (CYC_PERFIX + 32 +7)
            temp2 = deconv(yk_copy(x8:x8 + CYC_PERFIX + 31 +7), channel_16QAM(x15, :));
            yk_telda_intermediate = [yk_telda_intermediate, temp2];
            x15 = x15 + 1;
        end
        yk_telda = yk_telda_intermediate;
        %----------------------------------------------------------------------%
    else                % Wrong Choise
        fprintf('Please choose the channel type\n');
        return;
    end
    
    % Remove the cyclic effect, then de-map the signal, and pass throuth
    % the de-interleaver
    % Parameters
    yk_telda_remove_CYC = [];
    bk_received = [] ;
    received_bits = zeros(1,32);
    
    % Remove the cyclic extension
    for x9 = 1 : (CYC_PERFIX + 32) : (L * (CYC_PERFIX + 32))
        yk_telda_remove_CYC = [yk_telda_remove_CYC yk_telda(x9 + CYC_PERFIX:x9 + (CYC_PERFIX + 31))];
    end
    
    % Passing througth the fft block to invert the effect the ifft
    for x10 = 1:32:32 * L
        received_bits = yk_telda_remove_CYC(x10:x10 + 31);
        received_bits_fft = fft(received_bits, 32);
        bits_fft = received_bits_fft;
        
        %-----------------------------------------------------------------%
        % passing the data througth the de mapper
        bit_demapped = [];
        x11 = 1;
        for x12 = 1:32
            if (real(bits_fft(x12)) > 0 && imag(bits_fft(x12)) > 0)     % first quarter
                if (real(bits_fft(x12)) > 2 && imag(bits_fft(x12)) > 2) % 1010
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) > 2 && imag(bits_fft(x12)) < 2) % 1011
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 1;
                elseif (real(bits_fft(x12)) < 2 && imag(bits_fft(x12)) > 2) % 1110
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) < 2 && imag(bits_fft(x12)) < 2) % 1111
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 1;
                end
            elseif (real(bits_fft(x12)) < 0 && imag(bits_fft(x12)) > 0)     % second quarter
                if (real(bits_fft(x12)) < -2 && imag(bits_fft(x12)) > 2)    % 0010
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) < -2 && imag(bits_fft(x12)) < 2) % 0011
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 1;
                elseif (real(bits_fft(x12)) > -2 && imag(bits_fft(x12)) > 2) % 0110
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) > -2 && imag(bits_fft(x12)) < 2) % 0111
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 1;
                    bit_demapped(x11 + 3)    = 1;
                end
            elseif (real(bits_fft(x12)) < 0 && imag(bits_fft(x12)) < 0)     % third quarter
                if (real(bits_fft(x12)) < -2 && imag(bits_fft(x12)) < -2)    % 0000
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) < -2 && imag(bits_fft(x12)) > -2) % 0001
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 1;
                elseif (real(bits_fft(x12)) > -2 && imag(bits_fft(x12)) < -2) % 0100
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) > -2 && imag(bits_fft(x12)) > -2) % 0101
                    bit_demapped(x11)        = 0;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 1;
                end
            elseif (real(bits_fft(x12)) > 0 && imag(bits_fft(x12)) < 0)     % fourth quarter
                if (real(bits_fft(x12)) > 2 && imag(bits_fft(x12)) < -2)    % 1000
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) > 2 && imag(bits_fft(x12)) > -2) % 1001
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 0;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 1;
                elseif (real(bits_fft(x12)) < 2 && imag(bits_fft(x12)) < -2) % 1100
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 0;
                elseif (real(bits_fft(x12)) < 2 && imag(bits_fft(x12)) > -2) % 1101
                    bit_demapped(x11)        = 1;
                    bit_demapped(x11 + 1)    = 1;
                    bit_demapped(x11 + 2)    = 0;
                    bit_demapped(x11 + 3)    = 1;
                end
            end
            x11 = x11 + 4;
        end
        %-----------------------------------------------------------------%
        
        % Passing the data througth the de-interleaver
        bk_de_intrlv = matintrlv(bit_demapped, 16, 8);
        bk_de_intrlv_copy = bk_de_intrlv;
        
        % This section if the repetition happend
        if repetition == 3
            temp = [] ;
            x11 = 1;
            for x13 = 1:3:126
                NO_OF_ZEROS = nnz(~bk_de_intrlv(x13:x13 + 2)); % Get the number of zeros
                if NO_OF_ZEROS > 1
                    temp(x11) = 0;
                else
                    temp(x11) = 1;
                end
                x11 = x11 + 1;
            end
            bk_de_intrlv_copy = temp;
        end
        bk_received = [bk_received, bk_de_intrlv_copy];
    end
    [~, Rat1] = biterr(OFDM_SYM , bk_received);
    BER_16QAM = [BER_16QAM, Rat1];
    
    
    
end
% plotiing the data
semilogy(SNR, BER_16QAM)
grid on;
if CH_TP == 1
    fprintf('With Flad Fadding Channel\n');
    if repetition == 1
        title('BER vs. SNR for 16QAM Flat Fadding channel, without coding');
        xlabel('SNR (dB)');
        ylabel('BER');
    elseif repetition == 3
        title('BER vs. SNR for 16QAM Flat Fadding channel, with coding');
        xlabel('SNR (dB)');
        ylabel('BER');
    end
elseif CH_TP == 2
    fprintf('With Frequncy Selective Channel\n');
    if repetition == 1
        title('BER vs. SNR for 16QAM Frequency Selective channel, without coding');
        xlabel('SNR (dB)');
        ylabel('BER');
    elseif repetition == 3
        title('BER vs. SNR for 16QAM Frequency Selective channel, with coding');
        xlabel('SNR (dB)');
        ylabel('BER');
    end
end
