clc;clear;close all;

%-----File reading and data selection-------%
SN=input('SHEET NAME  ');
DATA=xlsread('FDS_MPMS_cave.xlsx',SN);  
tester = zeros(29,5);
testeq = zeros(29,5);
contadorx = 1;
contadory = 1;
Temp=DATA(:,1);
frq=[1;3.2;10;32;100];
w=2*pi*frq;

Xr_1=DATA(:,2);Xq_1=DATA(:,3);
Xr_2=DATA(:,4);Xq_2=DATA(:,5);
Xr_3=DATA(:,6);Xq_3=DATA(:,7);
Xr_4=DATA(:,8);Xq_4=DATA(:,9);
Xr_5=DATA(:,10);Xq_5=DATA(:,11);
Xr=[Xr_1 Xr_2 Xr_3 Xr_4 Xr_5]; 
Xq=[Xq_1 Xq_2 Xq_3 Xq_4 Xq_5]; 

%-----variables initialization and initial values------%
p=5;
N=length(Temp);
nof=N;
ti=-3;
tf=1;
n=5*(tf-ti);
TK=logspace(ti,tf,n);
Xr_obs=zeros(5,10);
Xq_obs=zeros(5,10);
XH=zeros(nof,1);
XL=zeros(nof,1);
DX=zeros(nof,1);
LFE=zeros(nof,1);
Ft=zeros(nof,1);

for nu=1:nof
    Xr_obs=Xr(nu,1:5)';
    Xq_obs=Xq(nu,1:5)'; 
    
    for ii=1:p
        for jj=1:n
            m1=(w(ii)*TK(jj))^2;m2=w(ii)*TK(jj);
            Gr(ii,jj)=1/(1+m1);Gq(ii,jj)=m2/(1+m1);
        end
    end
    
    %-------------Determination of p, dx and Xhf------------%
    q1=lsqnonneg(Gq,Xq_obs); 
    dx=-sum(q1); 
    r=Xr_obs+(Gr*q1);
    xhf=mean(r);
    X_hf=ones(p,1)*xhf; 
    
    %-----------Determination of normalized data vector d--------%
    dr=(Xr_obs-X_hf)/dx; % Normalized real susceptibility
    di=-Xq_obs/dx; % Normalized imaginary susceptibility
    wt=norm(Xr_obs)/norm(Xq_obs);  % Weight
    d=[dr;wt*di]; % Normalized data vector
    
    %-----------Determination of sensitivity matriz G-----------%
    G=[Gr;wt*Gq]; 
    
    %-------Regularization and constrains----------%
    [U,S,V]=svd(G);CG=cond(S);L=diag(S);
    v=ones(n,1);
    A=diag(v);
    b=zeros(n,1);
    d1=[d;b]; % A1*x = d1
    lb=zeros(n,1);
    ub=ones(n,1);  % lb ? x ? ub. Set lb = [] and b = [] if no bounds exist.
    x0=zeros(n,1); % starting point to x0. 
    Aeq=ones(1,n); % sum(gk)==1

    %------Smoothing - L curve--------%
    MI=logspace(-4,1,25);
    mi=max(L)*MI;
    x=4;
    mi_opt=mi(x);
    A1=[G;mi_opt*A];
    pk=lsqlin(A1,d1,[],[],Aeq,1,lb,ub,[]); 
    res=d-(G*pk);
        
    %--------------Calculate X' and X''---------------%
    q=dx*pk;
    Xr_calc=X_hf-(Gr*q);
    Xq_calc=-Gq*q;

    XH(nu)=xhf;
    DX(nu)=dx;
    
    %-------------RMS-------------%
    d_XH(nu)=rms(XH); % RMS of Xhf
    d_DX(nu)=rms(Xr_obs); % RMS of dX
        
    %-----Figures------%
    figure
    subplot(2,2,[1 3])
    semilogx(TK,abs(pk),'-w');hold on;stairs(TK,abs(pk),'b');
    hold off;ylabel('p')
    xlabel('Relaxation Time (s)');title(num2str(Temp(nu)));
    fontsize(gca,12,"points")
    subplot(2,2,2)
    loglog(frq, Xr_obs,'.k');hold on;loglog(frq,Xr_calc,'r');hold off;
    xlabel('Frequency (Hz)');ylabel('\chi^, (S/m)');
    fontsize(gca,12,"points")
    subplot(2,2,4)
    loglog(frq,Xq_obs,'.k');hold on;loglog(frq,Xq_calc,'r');hold off;
    xlabel('Frequency (Hz)');ylabel('\chi^{,,} (S/m)');
    fontsize(gca,12,"points")

end
