clc;clear;format long;

% Criacao do vetor tk e fk
tk = logspace(-10,1,1000);      % vetor de tempo tk 
tamanhotk = length(tk);         % tamanho do vetor tk 
tkmin = min(tk);
tkmax = max(tk);
fk = linspace(3,3.5,100);
tamanhofk = length(fk);

% Constants
t0 = 1e-9;  % segundos
Ms = 2.39e+5;   % A/m
m0 = 1.26e-6;   % H/m
Hc = 3.66e+4;   % A/m
Hk = Hc * 2.09; % A/m
Kb = 1.38e-23;  % J/K
T = [5;10;20;30;300];    % Temperature (K)
KbT = Kb.*T;
K = 1.15e+4;    % J/m3

% Particle sizes
diam = [1e-10; 19e-9];
raio = diam/2;
volume = (4/3)*pi*(raio.^3);

% KV (J)
KV = K.*volume;

% B
B10 = KV(1)./KbT;
B19 = KV(2)./KbT;

% t
t10 = t0*exp(B10);
t19 = t0*exp(B19);

%-----------Diameter of 10nm--------------------%
% Applying pk values to tk for T = 300 K
porc5 = t10(5)/20;
cont = 0;
pk10T300 = zeros(1,tamanhotk);   
for i=1:tamanhotk
    a = abs(t10(5)-tk(i));
    if a <= porc5
        pk10T300(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk10T300(i) == 1
        pk10T300(i) = pico;
    end
end

% Applying pk values to tk for T = 30 K
porc5 = t10(4)/20;
cont = 0;
pk10T30 = zeros(1,tamanhotk); 
for i=1:tamanhotk
    a = abs(t10(4)-tk(i));
    if a <= porc5
        pk10T30(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk10T30(i) == 1
        pk10T30(i) = pico;
    end
end

% Applying pk values to tk for T = 20 K
porc5 = t10(3)/20;
cont = 0;
pk10T20 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t10(3)-tk(i));
    if a <= porc5
        pk10T20(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk10T20(i) == 1
        pk10T20(i) = pico;
    end
end

% Applying pk values to tk for T = 5 K
porc5 = t10(2)/20;
cont = 0;
pk10T10 = zeros(1,tamanhotk);   
for i=1:tamanhotk
    a = abs(t10(2)-tk(i));
    if a <= porc5
        pk10T10(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk10T10(i) == 1
        pk10T10(i) = pico;
    end
end

% Applying pk values to tk for T = 5 K
porc5 = t10(1)/20;
cont = 0;
pk10T5 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t10(1)-tk(i));
    if a <= porc5
        pk10T5(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk10T5(i) == 1
        pk10T5(i) = pico;
    end
end

%---------Diameter of 19nm----------------------------%
% Applying pk values to tk for T = 300 K
porc5 = t19(5)/20;
cont = 0;
pk19T300 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t19(5)-tk(i));
    if a <= porc5
        pk19T300(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk19T300(i) == 1
        pk19T300(i) = pico;
    end
end

% Applying pk values to tk for T = 30 K
porc5 = t19(4)/20;
cont = 0;
pk19T30 = zeros(1,tamanhotk);   
for i=1:tamanhotk
    a = abs(t19(4)-tk(i));
    if a <= porc5
        pk19T30(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk19T30(i) == 1
        pk19T30(i) = pico;
    end
end

% Applying pk values to tk for T = 20 K
porc5 = t19(3)/20;
cont = 0;
pk19T20 = zeros(1,tamanhotk);   
for i=1:tamanhotk
    a = abs(t19(3)-tk(i));
    if a <= porc5
        pk19T20(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk19T20(i) == 1
        pk19T20(i) = pico;
    end
end

% Applying pk values to tk for T = 10 K
porc5 = t19(2)/20;
cont = 0;
pk19T10 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t19(2)-tk(i));
    if a <= porc5
        pk19T10(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk19T10(i) == 1
        pk19T10(i) = pico;
    end
end

% Applying pk values to tk for T = 35 K
porc5 = t19(1)/20;
cont = 0;
pk19T5 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t19(1)-tk(i));
    if a <= porc5
        pk19T5(i) = 1;
        cont = cont + 1;
    end
end
pico = 1/cont;
for i=1:tamanhotk
    if pk19T5(i) == 1
        pk19T5(i) = pico;
    end
end

%--------------Comparative figure for both sizes for different temperatures------------------%
%T = 5K
figure, semilogx(tk,pk10T5,tk,pk19T5)
xlabel('Tempo de relaxação (s)')
xlim([tkmin tkmax])
ylabel('P')
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 19nm',Location='northeast')
grid on
title('T = 5K')

%T = 10K
figure, semilogx(tk,pk10T10,tk,pk19T10)
xlabel('Tempo de relaxação (s)')
xlim([tkmin tkmax])
ylabel('P')
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 19nm',Location='northeast')
grid on
title('T = 10K')

%T = 20K
figure, semilogx(tk,pk10T20,tk,pk19T20)
xlabel('Tempo de relaxação (s)')
xlim([tkmin tkmax])
ylabel('P')
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 19nm',Location='northeast')
grid on
title('T = 20K')

%T = 30K
figure, semilogx(tk,pk10T30,tk,pk19T30)
xlabel('Tempo de relaxação (s)')
xlim([tkmin tkmax])
ylabel('P')
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 19nm',Location='northeast')
grid on
title('T = 30K')

%T = 300K
figure, semilogx(tk,pk10T300,tk,pk19T300)
xlabel('Tempo de relaxação (s)')
xlim([tkmin tkmax])
ylabel('P')
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 19nm',Location='northeast')
grid on
title('T = 300K')