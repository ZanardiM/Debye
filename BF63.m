clc;clear;close all;

%---------file reading and data selection-----------%
nof=input('number of files to be inverted   ');
Temp=zeros(nof,1);
T1=zeros(nof,1);
tau=logspace(-6,-2,80);
m=length(tau);
Z=zeros(nof,1);
XH=zeros(nof,1);
XL=zeros(nof,1);
DX=zeros(nof,1);
F=zeros(nof,1);
for nu=1:nof
    FN=input('FILE NAME  ');
    D=xlsread(FN);
    frq=D(:,2);
    Xr_obs=D(:,3);
    Xq_obs=D(:,4);
    n=length(frq); 
    T = mean(D(:,1)); % temperature mean
    Temp(nu)=T;
    T1(nu)=1/T;
    mu=0.0002;
    omg= 6.2832*frq;

%----------Creation and calculation of matrices---------%
    Gr=zeros(n,m);
    Gq=zeros(n,m);
    for j=1:m
        ot=tau(j)*omg;
        ot2=ot.*ot;
        pr=1./(1+ot2);
        pq=ot.*pr;
        Gr(:,j)=pr;
        Gq(:,j)=pq;
    end
    M=eye(m,m)*mu;
    o=zeros(m,1);
    lb=zeros(m,1);
    ub=ones(m,1);
    A=eye(m);
    b=ones(m,1);
    Aeq=ones(1,m);
    beq=1;
    G0=[Gq;M];

    %-------------Calculation of p, dx and Xsd-------------%
    p0=lsqlin(G0,[Xq_obs;o],[],[],[],[],lb,1000*ub);
    dX=norm(p0,1);
    p=p0/dX;
    r=Xr_obs-dX*Gr*p;
    Xsd=mean(r);
    x_sd=ones(n,1)*Xsd; 
    dXsd=std(r);
    Xrn=(Xr_obs-Xsd)/dX;Xqn=Xq_obs/dX;
    G=[Gr;Gq;M];
    d=[Xrn;Xqn;o];
    p=lsqlin(G,d,[],[],ones(1,m),1,lb,1*ub);
    Xr_calc=Xsd+(dX*Gr*p);Xq_calc=dX*Gq*p; % X' and X'' obtained through inversion   

    XH(nu)=Xsd;
    Xsp=dX+Xsd;
    DX(nu)=dX;
    XL(nu)=Xsp;
 
%-------p distribution, X' and X'' plotting------------%   
    figure
    subplot(2,2,[1 3])
    semilogx(tau,abs(p),'-w');hold on;stairs(tau,abs(p),'b');
    fontsize(gca,12,"points")
    hold off;ylabel('p')
    xlabel('Relaxation Time (s)');title(num2str(Temp(nu)));
    subplot(2,2,2)
    loglog(frq, Xr_obs,'.k');hold on;loglog(frq,Xr_calc,'r');hold off;
    fontsize(gca,12,"points")
    xlabel('Frequency (Hz)');ylabel('\chi^, (S/m)');
    subplot(2,2,4)
    loglog(frq,Xq_obs,'.k');hold on;loglog(frq,Xq_calc,'r');hold off;
    xlabel('Frequency (Hz)');ylabel('\chi^{,,} (S/m)');
    fontsize(gca,12,"points")
    plottools 'on'

end
