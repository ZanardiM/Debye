clc;clear; format long

%-----Frequency array--------%
freq = logspace(-70,20,1000);   % size and gap of array 
minfreq = min(freq);
maxfreq = max(freq);
fw = length(freq);
w = 2*pi*freq;  % angular frequency array

%-----Diameter array--------%
diam = transpose(linspace(10,22,13));  % diameter array size and gap in nm
tamanho = length(diam);
raio = (diam/2)*1e-9;

% Volume
volume = (4/3)*pi*(raio.^3);

% Constants
t0 = 1e-9;  % s
Ms = 2.39e+5;   % A/m
m0 = 1.26e-6;   % H/m
Hc = 3.66e+4;   % A/m
Hk = Hc * 2.09; % A/m
Kb = 1.38e-23;  % J/K
T = 30;    % K (this is where you chance T value)
KbT = Kb*T;
K = 1.15e+4;    % J/m3

% KV (J)
KV = K*volume;

% B
B = KV/KbT;

% t
t = t0*exp(B);

% Xhf
Xhf = (m0*volume*Ms)/(3*KbT);

% Xlf
Xlf = (2*Ms)/(3*Hk);

% ΔX
deltaX = Xlf - Xhf;

% wt e wt²
wt = w.*t;
wt2 = wt.^2;

% X' using X' = Xhf + ΔX/(1+(wt)²)
a = 1./(1.+wt2);
b = deltaX .* a;
Xr = zeros(tamanho,fw);
for i=1:tamanho
    for j=1:fw
        Xr(i,j) = Xhf(i) + b(i,j);
    end
end

% X' curves for  d=10, d=13, d=16, d=19 e d=22  
figure,loglog(freq,Xr(1,:),freq,Xr(4,:),freq,Xr(7,:),freq,Xr(9,:),freq,Xr(13,:))
xlabel("Frequência (Hz)")
xlim([minfreq maxfreq])
fontsize(gca,12,"points")
ylabel("X' (SI)")
legend('d = 10nm', 'd = 13nm', 'd = 16nm', 'd = 19nm', 'd = 22nm',Location='northeast')
grid on
title("T = 30K")

% X'' using X'' = ΔXwt/(1+(wt)²)
ai = wt./(1.+wt2);
bi = deltaX .* ai;
Xi = zeros(tamanho,fw);
for i=1:tamanho
    for j=1:fw
        Xi(i,j) = bi(i,j);
    end
end

% X'' curves for d=10, d=13, d=16, d=19 e d=22
figure,loglog(freq,Xi(1,:),freq,Xi(4,:),freq,Xi(7,:),freq,Xi(9,:),freq,Xi(13,:))
xlabel("Frequência (Hz)")
xlim([minfreq maxfreq])
ylim([1e-7 1e+1])
fontsize(gca,12,"points")
ylabel("X'' (SI )")
legend('d = 10nm', 'd = 13nm', 'd = 16nm', 'd = 19nm', 'd = 22nm',Location='southwest')
grid on