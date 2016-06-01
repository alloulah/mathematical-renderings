%
% Rise and Decay
%
% Mo Alloulah
%
% 
%
% 010616 - v1.0  - Initial issue version. 
%
%% Offset Multisine parameters

M = 7;
Fc = 0;
m = [-(M-1)/2:(M-1)/2];
Fg = 2e3;
G = 200;
Fd = G*Fg;
p = [0 -6 -11 -14 -11 -1 43];

fm_hat = Fc + m*Fd + p*Fg;

BW = 2.4e6;
N = BW/Fg;

%% synthesis

toneNdx_vect = fm_hat/Fg;

negIndices = find( toneNdx_vect < 0 );
toneNdx_vect(negIndices) = toneNdx_vect(negIndices) + N-1;     % translate negative indices

R = 19;
Nu = N;

pilotHd = 10^(5/20);

offsetMultisine_f = zeros(N, 1);
offsetMultisine_f(toneNdx_vect+1) = pilotHd;
offsetMultisine0_f = [offsetMultisine_f(1:N/2); repmat(zeros(N,1), R, 1); offsetMultisine_f(N/2+1:N);];
offsetMultisine_t  = ((R+1)*N/sqrt(Nu)) * ifft(offsetMultisine0_f);

%% animate rise and decay

fontSize = 9;

prnt0 = figure(256+0);
h0 = subplot(1, 1, [1]);

% two degrees of freedom to control animation: sweepStride & sweepRefreshDelay
ccSweep = [0.0125:0.0125:16];
sweepStride = 4;
sweepPace = [1:sweepStride:length(ccSweep)];

sweepRefreshDelay = 10e-6;

for cc = ccSweep(sweepPace)
cla(h0);
x = offsetMultisine_t;
y = 4 * ( 1./(1+exp(-cc*x)) - 0.5 );

y0_f = (sqrt(Nu)/((R+1)*N)) * fft( y );
y0_f = y0_f/pilotHd;
sigPwr_dB = 10*log10(y0_f'*y0_f/length(y0_f));
y1_f = 10^((0-sigPwr_dB)/20)*y0_f;
y_f = abs(fftshift(y1_f));

X = [-length(x)/2:length(x)/2-1];
Y0 = 20*log10(y_f);
Y1 = 4 * ( 1./(1+exp(-cc*X/length(X))) - 0.5 );
[ax0, h00, h01] = plotyy(h0, X, Y0, X/length(X), Y1, 'plot');
set(ax0, {'ycolor'}, {'b'; 'r'});
set(ax0(1), 'XLim', [min(X) max(X)]);
set(ax0(1), 'YLim', [-350 max(Y0)]);
set(ax0(2), 'XLim', [min(X) max(X)]/length(X));
set(ax0(2), 'YLim', [-2 2]);
yTick = [-350 -250 -150 -50 50];
set(ax0(1), 'YTick', yTick);
yTick = [-2 -1 0 1 2];
set(ax0(2), 'YTick', yTick);
set(get(ax0(1), 'Ylabel'), 'String', 'noise (dB)', 'FontSize', fontSize)
set(get(ax0(2), 'Ylabel'), 'String', 'nonlinearity', 'FontSize', fontSize)
set(h00, 'LineStyle', 'none', 'Marker', '.', 'MarkerEdgeColor', 'b', 'MarkerSize', fontSize)
set(h01, 'LineStyle', '-', 'Marker', 'none', 'MarkerEdgeColor', 'r', 'MarkerSize', fontSize)
ttl = title(h0, 'rise and decay');
set(ttl, 'FontSize', fontSize+1);

grid(h0);
pause(sweepRefreshDelay)
end
