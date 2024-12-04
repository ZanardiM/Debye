clc;clear;close all;

%-----File reading and data selection-------%
SN=input('SHEET NAME  '); 
DATA=xlsread('FDS_MPMS.xlsx',SN);      
T=DATA(:,1);F=DATA(:,2);
Xr_OBS=DATA(:,3);
dXr_OBS=DATA(:,4);
Xq_OBS=DATA(:,5);
dXq_OBS=DATA(:,6);

%-----variables initialization and initial values------%
p=5;
N=length(T);
nof=N/p;
ti=-3;
tf=1;
n=5*(tf-ti);
TK=logspace(ti,tf,n)';
Temp=zeros(nof,1);
XH=zeros(nof,1);
XL=zeros(nof,1);
DX=zeros(nof,1);
LFE=zeros(nof,1);
Ft=zeros(nof,1);

for x=0:nof-1
    Temp(x+1)=T((x*5)+1);
end

for nu=1:nof
    Xr_obs=Xr_OBS((p*nu)-(p-1):p*nu);  
    Xq_obs=Xq_OBS((p*nu)-(p-1):p*nu);  
    frq=F((p*nu)-(p-1):p*nu);  
    w=2*pi*frq;
    
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
    b=zeros(n,1);d1=[d;b]; % A1*x = d1
    lb=zeros(n,1);ub=ones(n,1); % lb ? x ? ub. Set lb = [] and b = [] if no bounds exist.
    x0=zeros(n,1); % starting point to x0. 
    Aeq=ones(1,n);% sum(gk)==1. 

    %------Smoothing - L curve--------%
    MI=logspace(-4,1,25);
    mi=max(L)*MI;
    x=4;
    mi_opt=mi(x);
    A1=[G;mi_opt*A];
    pk=lsqlin(A1,d1,[],[],Aeq,1,lb,ub,[]); 
    res=d-(G*pk);
    
    %------Calculate X' and X''------%
    q=dx*pk;
    Xr_calc=X_hf-(Gr*q);
    Xq_calc=-Gq*q;
     
    %---------RMS---------%
    d_XH(nu)=rms(XH); % RMS of Xhf
    d_DX(nu)=rms(Xr_obs); % RMS of dX
    
    %-----Figures-----%    
    figure
    subplot(2,2,[1 3])
    semilogx(TK,abs(pk),'-w');hold on;stairs(TK,abs(pk),'b');
    fontsize(gca,12,"points")
    hold off;ylabel('p')
    xlabel('Relaxation Time (s)');title(num2str(Temp(nu)));
    subplot(2,2,2)
    loglog(frq, Xr_obs,'.k');hold on;loglog(frq,Xr_calc,'r');hold off;
    fontsize(gca,12,"points")
    xlabel('Frequency (Hz)');ylabel('\chi^, (S/m)');
    subplot(2,2,4)
    loglog(frq,Xq_obs,'.k');hold on;loglog(frq,Xq_calc,'r');hold off;
    fontsize(gca,12,"points")
    xlabel('Frequency (Hz)');ylabel('\chi^{,,} (S/m)');
     
end
    
