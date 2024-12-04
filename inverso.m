clc;clear;format long;

%-------Synthetic data calculation----------%

% Frequency array
freq = logspace(1,20,2000);   % frequency array
minfreq = min(freq);
maxfreq = max(freq);
fw = length(freq);
w = 2*pi*freq;  % angular frequency array

% Particle sizes
diam = [10e-9, 13e-9, 16e-9, 19e-9, 22e-9];
raio = diam/2;
raio1 = diam(1)/2;
raio2 = diam(2)/2;
raio3 = diam(3)/2;
raio4 = diam(4)/2;
raio5 = diam(5)/2;

% Volume
volume = (4/3)*pi*(raio.^3);

% Constants
t0 = 1e-9;      % segundos
Ms = 2.39e+5;   % A/m
m0 = 1.26e-6;   % H/m
Hc = 3.66e+4;   % A/m
Hk = Hc * 2.09; % A/m
Kb = 1.38e-23;  % J/K
T = 300;        % K
KbT = Kb*T; % J
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

% wt and (wt)²
transT = transpose(t);
wt = w.*transT;
wt2 = wt.^2;

% X' real using X' = Xhf + DX/(1+wt2)
a = 1./(1.+wt2);
ta = transpose(a);
b = deltaX .* ta;
Xr = Xhf + b;

% X'' using X'' = ΔXwt/(1+(wt)²)
ai = wt./(1.+wt2);
tai = transpose(ai);
bi = deltaX .* tai;
Xi = bi;

%------------------Inversion method-------------------%

% Time array
tk = logspace(-10,1,1000);   
tamanhotk = length(tk);     
tkmin = min(tk);
tkmax = max(tk);

% Mr e Mi matrices
Mr = zeros(fw,tamanhotk);
Mi = zeros(fw,tamanhotk);
for i=1:fw
    for j=1:tamanhotk
        denominador = (w(i)*tk(j))^2;
        denominador = denominador + 1;
        Mr(i,j) = 1/denominador;
        Mi(i,j) = (w(i)*tk(j))/denominador;
    end
end

% ΔX parameter 
di = Xi;
di1 = Xi(:,1);
di2 = Xi(:,2);
di3 = Xi(:,3);
di4 = Xi(:,4);
di5 = Xi(:,5);
q = linsolve(Mi,di);
q1 = lsqnonneg(Mi,di1);
q2 = lsqnonneg(Mi,di2);
q3 = lsqnonneg(Mi,di3);
q4 = lsqnonneg(Mi,di4);
q5 = lsqnonneg(Mi,di5);
dXinv = norm(q);
dXinv1 = norm(q1);
dXinv2 = norm(q2);
dXinv3 = norm(q3);
dXinv4 = norm(q4);
dXinv5 = norm(q5);

% Xhf parameter
dr = Xr;
r = dr - (Mr*q);
Xhfinv = mean(r);

% p distribution
p = q./dXinv;
p1 = q1./dXinv1;
p2 = q2./dXinv2;
p3 = q3./dXinv3;
p4 = q4./dXinv4;
p5 = q5./dXinv5;

% p distributin figure using inversion data
figure, semilogx(tk,p1)
hold on
semilogx(tk,p2)
semilogx(tk,p3)
semilogx(tk,p4)
semilogx(tk,p5)
hold off
xlim([1e-9 1e-1])
xlabel('Tempo de relaxação (s)')
ylabel('P')
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 13nm', 'd = 16nm', 'd = 19nm', 'd = 22nm',Location='southwest')

% cross-plot of calculated and simulated ΔX
dxinv = [dXinv1,dXinv2,dXinv3,dXinv4,dXinv5];
figure, loglog(deltaX,dxinv, "ro")
xlabel('delta X calculado')
ylabel('delta X simulado')


% X' using inversion data
a = 1./(1.+wt2);
ta = transpose(a);
b = dXinv .* ta;
Xrinv = Xhfinv + b;

% X' curve comparing calculated and simulated values
figure, loglog(freq,Xr(:,1),freq,Xr(:,2),freq,Xr(:,3),freq,Xr(:,4),freq,Xr(:,5))
xlim([minfreq maxfreq])
%xlim([1e+1 1e+12])
hold on
loglog(freq,Xrinv(:,1), "o")
hold on
loglog(freq,Xrinv(:,2), "o")
hold on
loglog(freq,Xrinv(:,3), "o")
hold on
loglog(freq,Xrinv(:,4), "o")
hold on
loglog(freq,Xrinv(:,5), "o")
xlabel('Tempo de relaxação (s)')
ylabel("X' (SI )")
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 13nm', 'd = 16nm', 'd = 19nm', 'd = 22nm', 'dinv = 10nm', 'dinv = 13nm', 'dinv = 16nm', 'dinv = 19nm', 'dinv = 22nm' ,Location='northeast')