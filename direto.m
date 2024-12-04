clc;clear;format long;

%------------Data reading-------------%
data = importdata("dados.txt"); % file with frequency and relaxion time 
t = data(:,1); 
freq = data(:,2); 

data2 = importdata("freq.txt"); % file with ΔX e Xhf
Xhf = data2(:,1);  
deltaX = data2(:,2); 

data3 = importdata("sus.txt");  % file with X' and X''
Xr = data3(:,1); % X'
Xi = data3(:,2); % X''

%-------------Array of time and frequency-----------%
tk = logspace(-10,1,1000);    
tamanhotk = length(tk);     
tkmin = min(tk);
tkmax = max(tk);
fk = linspace(3,3.5,100);
tamanhofk = length(fk);

%----------Applying pk values to tk for d = 10nm------------------%
porc5 = t(1)/20;
cont10 = 0;
pk10 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(1)-tk(i));
    if a <= porc5
        pk10(i) = 1;
        cont10 = cont10 + 1;
    end
end
pico10 = 1/cont10;
for i=1:tamanhotk
    if pk10(i) == 1
        pk10(i) = pico10;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk10)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 10nm')

%----------Applying pk values to tk for d = 11nm------------------%
porc5 = t(2)/20;
cont11 = 0;
pk11 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(2)-tk(i));
    if a <= porc5
        pk11(i) = 1;
        cont11 = cont11 + 1;
    end
end
pico11 = 1/cont11;
for i=1:tamanhotk
    if pk11(i) == 1
        pk11(i) = pico11;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk11)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 11nm')

%----------Applying pk values to tk for d = 12nm------------------%
porc5 = t(3)/20;
cont12 = 0;
pk12 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(3)-tk(i));
    if a <= porc5
        pk12(i) = 1;
        cont12 = cont12 + 1;
    end
end
pico12 = 1/cont12;
for i=1:tamanhotk
    if pk12(i) == 1
        pk12(i) = pico12;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk12)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 12nm')

%----------Applying pk values to tk for d = 13nm------------------%
porc5 = t(4)/20;
cont13 = 0;
pk13 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(4)-tk(i));
    if a <= porc5
        pk13(i) = 1;
        cont13 = cont13 + 1;
    end
end
pico13 = 1/cont13;
for i=1:tamanhotk
    if pk13(i) == 1
        pk13(i) = pico10;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk13)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 13nm')

%----------Applying pk values to tk for d = 14nm------------------%
porc5 = t(5)/20;
cont14 = 0;
pk14 = zeros(1,tamanhotk);   
for i=1:tamanhotk
    a = abs(t(5)-tk(i));
    if a <= porc5
        pk14(i) = 1;
        cont14 = cont14 + 1;
    end
end
pico14 = 1/cont14;
for i=1:tamanhotk
    if pk14(i) == 1
        pk14(i) = pico14;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk14)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 14nm')

%----------Applying pk values to tk for d = 15nm------------------%
porc5 = t(6)/20;
cont15 = 0;
pk15 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(6)-tk(i));
    if a <= porc5
        pk15(i) = 1;
        cont15 = cont15 + 1;
    end
end
pico15 = 1/cont15;
for i=1:tamanhotk
    if pk15(i) == 1
        pk15(i) = pico15;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk15)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 15nm')

%----------Applying pk values to tk for d = 16nm------------------%
porc5 = t(7)/20;
cont16 = 0;
pk16 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(7)-tk(i));
    if a <= porc5
        pk16(i) = 1;
        cont16 = cont16 + 1;
    end
end
pico16 = 1/cont16;
for i=1:tamanhotk
    if pk16(i) == 1
        pk16(i) = pico16;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk16)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 16nm')

%----------Applying pk values to tk for d = 17nm------------------%
porc5 = t(8)/20;
cont17 = 0;
pk17 = zeros(1,tamanhotk);   
for i=1:tamanhotk
    a = abs(t(8)-tk(i));
    if a <= porc5
        pk17(i) = 1;
        cont17 = cont17 + 1;
    end
end
pico17 = 1/cont17;
for i=1:tamanhotk
    if pk17(i) == 1
        pk17(i) = pico17;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk17)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 17nm')

%----------Applying pk values to tk for d = 18nm------------------%
porc5 = t(9)/20;
cont18 = 0;
pk18 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(9)-tk(i));
    if a <= porc5
        pk18(i) = 1;
        cont18 = cont18 + 1;
    end
end
pico18 = 1/cont18;
for i=1:tamanhotk
    if pk18(i) == 1
        pk18(i) = pico18;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk18)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 18nm')

%----------Applying pk values to tk for d = 19nm------------------%
porc5 = t(10)/20;
cont19 = 0;
pk19 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(10)-tk(i));
    if a <= porc5
        pk19(i) = 1;
        cont19 = cont19 + 1;
    end
end
pico19 = 1/cont19;
for i=1:tamanhotk
    if pk19(i) == 1
        pk19(i) = pico19;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk19)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 19nm')

%----------Applying pk values to tk for d = 20nm------------------%
porc5 = t(11)/20;
cont20 = 0;
pk20 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(11)-tk(i));
    if a <= porc5
        pk20(i) = 1;
        cont20 = cont20 + 1;
    end
end
pico20 = 1/cont20;
for i=1:tamanhotk
    if pk20(i) == 1
        pk20(i) = pico20;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk20)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 20nm')

%----------Applying pk values to tk for d = 21nm------------------%
porc5 = t(12)/20;
cont21 = 0;
pk21 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(12)-tk(i));
    if a <= porc5
        pk21(i) = 1;
        cont21 = cont21 + 1;
    end
end
pico21 = 1/cont21;
for i=1:tamanhotk
    if pk21(i) == 1
        pk21(i) = pico21;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk21)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 21nm')

%----------Applying pk values to tk for d = 22nm------------------%
porc5 = t(13)/20;
cont22 = 0;
pk22 = zeros(1,tamanhotk);    
for i=1:tamanhotk
    a = abs(t(13)-tk(i));
    if a <= porc5
        pk22(i) = 1;
        cont22 = cont22 + 1;
    end
end
pico22 = 1/cont22;
for i=1:tamanhotk
    if pk22(i) == 1
        pk22(i) = pico20;
    end
end
%----p distribution chart-----%
figure, semilogx(tk,pk22)
xlabel('log (Tempo de relaxação [s])')
xlim([tkmin tkmax])
ylabel('Pk')
title('Diâmetro = 22nm')

%--------p distribution chart for 10, 13, 16, 19, 22 nm----------%
figure, semilogx(tk,pk10,tk,pk13,tk,pk16,tk,pk19,tk,pk22)
xlabel('Tempo de relaxação (s)')
xlim([1e-9 1e-1])
ylim([0 0.4])
ylabel('P')
title('(A)')
fontsize(gca,12,"points")
legend('d = 10nm', 'd = 13nm', 'd = 16nm', 'd = 19nm', 'd = 22nm',Location='northwest')
grid on