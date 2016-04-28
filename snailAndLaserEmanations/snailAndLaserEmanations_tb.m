%
% The snail and laser emanations
%
% Mo Alloulah
%
% 
%
% 270416 - v1.0  - Initial issue version. 
%
%% OFDM Paramters (DVB-T2 inspired)

Nsym = 3;                                   % # of OFDM symbols
Nu = 27841;                                 % # of atctive carriers (utilization)
N  = 32768;                                 % FFT size
Ng = 32768*19/256;                          % guard interval (samples)
Ns = N + Ng;                                % OFDM symbol size
eqLimit = 3;                                % equalizer magnitude limit
constellation = 'BPSK';                     % constellation
constMapTable = [1 -1];                     % constellation map

%% Channel Parameters

Ts = 7/64e6;                                % elementary period
echoDelays_us = Ng*Ts/1e-6 + [0:150-1];     % echo dealys (us)
echoDelays = round(echoDelays_us*1e-6./Ts); % echo delays (samples)
echoAtt = 0.0;                              % echo attenulation (dB)
chSNRdB = inf;                              % channel SNR (dB)

%% Tx

tx_t = zeros(Ns*Nsym, 1);
refData_f = zeros(Nu, Nsym);
for n = 1:Nsym
    % 1 impulse subcarrier, 0 rest
    ndx = floor([1; 0.5*ones(Nu-1,1);]*length(constMapTable));
    txData_f = constMapTable(ndx).';
    txData_f(2:end) = zeros(Nu-1,1);
    refData_f(:,n) = txData_f;
    
    % zero padding inactive subcarriers
    txDataZeroPadded_f = [txData_f; zeros(N-Nu,1)];
    txDataZeroPadded_f = circshift(txDataZeroPadded_f, -(Nu-1)/2);

    % IFFT
    txData_t = (N/sqrt(Nu)) * ifft(txDataZeroPadded_f);

    % cyclic prefix
    txDataGi_t = [txData_t(end-Ng+1:end); txData_t];

    tx_t(1+(n-1)*Ns:n*Ns) = txDataGi_t;
end

%% Rx - iterate over channel delays

iciErr_mtrx = zeros(Nu, length(echoDelays));

for dd = 1:length(echoDelays)
    % channel delay
    echoDly = echoDelays(dd);
    
    % add multipath
    ch_t = tx_t + circshift(tx_t, echoDly) * 10^(-echoAtt/20);
    ch_t = ch_t / sqrt(1+10^(-echoAtt/10));

    % add AWGN
    rx_t = ch_t + 10^(-chSNRdB/20)*(randn(size(ch_t)) + 1i*randn(size(ch_t)))/sqrt(2) * sqrt(N/Nu);

    % CIR estimation
    CIR = zeros(N, 1);
    CIR(1) = 1;
    CIR(1+echoDly) = 10^(-echoAtt/20);
    CIR = CIR/sqrt(1+10^(-echoAtt/10));
    CFR = fft(CIR);
    CFRu = circshift(CFR, +(Nu-1)/2);
    CFRu = CFRu(1:Nu);
    
    % EOG-induced ICI
    iciErr = zeros(Nu, Nsym-2);
    
    % loop over symbols
    for ss = 2:Nsym-1

        % GI remove and align to right guard edge    
        rxSym_t = rx_t((ss-1)*Ns+Ng+1:ss*Ns);
        
        % hard-decision of previous symbol
        hdData = refData_f(:, ss-1);
        
        % inter-block (IBI) cancellation
        rxSym0_t = rxSym_t;
        if (echoDly-Ng > 0)
            % synthesize IBI
            preData_f = circshift([hdData; zeros(N-Nu,1)], -(Nu-1)/2);
            preData_t = ifft(preData_f);
            preData_t(1:N/2-1) = 0;
            preData_f = fft(preData_t);
            preRxData_f = preData_f .* CFR;
            preRxData_t = (N/sqrt(Nu)) * ifft(preRxData_f);
            
            % remove spill-over from the previous symbol
            rxSym0_t(1:echoDly-Ng) = rxSym_t(1:echoDly-Ng) - preRxData_t(1+Ng:echoDly);
        end
        
        % demod to examine ISI-induced ICI
        % FFT
        rxData_f = (N/sqrt(Nu))^-1 * fft(rxSym0_t);
        
        % equalize [gently] with rotation
        eqData_f = (rxData_f .* exp(-1i*angle(CFR))) ./ (abs(CFR) + eps);
        
        % limit maximum magnitude
        contourNdx = find( abs(eqData_f) > eqLimit );
        eqData_f(contourNdx) = eqLimit * exp(-1i*angle(eqData_f(contourNdx)));
        
        % extract active carriers
        eqDataAc_f = circshift(eqData_f, +(Nu-1)/2);
        eqDataAc_f = eqDataAc_f(1:Nu);
        
        % measure ICI 
        ici_err0 = eqDataAc_f - refData_f(:,ss);
        ici_err1 = ici_err0 .* abs(CFRu);
        iciErr(:, ss-1) = ici_err1;
    
    end % for symbols
    
    % pack ICI error
    iciErr_mtrx(:, dd) = iciErr;

end % for delays

%% visualization - here we go

echoNdx = 2;                                % snail index
y = iciErr_mtrx(:, echoNdx);

figure();
plot( y, '.' )
ttl = title(gca, 'snail');
set(ttl, 'FontSize', 10);

[N, L] = size(iciErr_mtrx);
Z_mtrx = zeros(N, L);
for k = 1:L
    Z_mtrx(:, k) = 20*log10(abs(fftshift(fft(iciErr_mtrx(:,k)))));
end
Z0 = Z_mtrx(:, 2:end);
surfView = [0 -90];
prnt0 = figure();
h0 = subplot(1, 1, [1]);
hz0 = surf( h0, Z0 );
axis(h0, [1 L-1 1 N]);
view(h0, surfView);
set(hz0, 'LineStyle', 'none');

xlbl = xlabel(h0, 'intensity');
set(xlbl, 'FontSize', 10);
ylbl = ylabel(h0, 'broadside');
set(ylbl, 'FontSize', 10);
ttl = title(h0, 'laser emanations');
set(ttl, 'FontSize', 10);
colormap(h0, 'hsv')                         % or 'hot'
