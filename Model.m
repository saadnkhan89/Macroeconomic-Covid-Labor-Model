%% Macroeconomic-COVID Model Simulation Without Government Policy
clc
clear all
warning('off')
%% Setup

% Number of samples
S=500;

%Number of time intervals
N=151;

% mimimum number of hours required (No governmnet policy)
hmin=0.125;

% mimimum number of hours required (governmnet policy)
%hmin=10^-8;

% lower and upper mortality bound (No government Policy)
mlow=50;
mhigh=2600;

% Hospital capacity
hoscap=0.05;

% penalty weights
mem=8;
weight=280/761;
w=weight./(1:mem);

% Maximum penalty
penmax=sqrt(1/(2*hmin));

% Initial population
Nint=8.4*10^6;

%% Initialization

% lag coefficient for initialization
lag=10;

% Setting matrices for variables

% Total population
tot=zeros(S,N+lag);
tot(1:S,1:lag+1)=Nint;

% Number of infections per time period
ifd=zeros(S,N+lag);
ifd(1:S,lag+1)=1;

% Cumulative infections
totifd=zeros(S,N+lag);
totifd(1:S,lag+1)=1;

% Recovery rate
rec=zeros(S,N+lag);
rec(1:S,lag+1)=1;

% Mortality rate
mor=zeros(S,N+lag);

% Labor supply 
opttim=zeros(S,N+lag);
opttim(1:S,1:lag+1)=0.5;

% Total Mortality
totdth=zeros(S,N+lag);

% Penalty at each time step
pen=zeros(S,N+lag);
pen(1:S,1:lag+1)=1;
% A=[tot;ifd;rec;mor;totdth;pen;opttim];

%% Simulating the model

for k=1:S
for i=lag+2:N+lag
    if opttim(k,i-1)>=0.4 && opttim(k,i-1)<=0.5
        ifd(k,i)=round(max(2,normrnd(3,0.5))*ifd(k,i-1));
    elseif opttim(k,i-1)>=0.2 && opttim(k,i-1)<0.4
        ifd(k,i)=round(max(1.5,normrnd(2,0.25))*ifd(k,i-1));
    else 
        ifd(k,i)=round(min(1,max(0.25,normrnd(0.5,0.25)))*ifd(k,i-1));
    end
    
    totifd(k,i)=totifd(k,i-1)+ifd(k,i);
        
   if ifd(k,i)/tot(k,i)<hoscap && ifd(k,i-1)/tot(k,i-1)<hoscap && ifd(k,i-2)/tot(k,i-2)<hoscap     
        mor(k,i)=round((0.05+normrnd(.025,.01))*ifd(k,i));
   else
       mor(k,i)=round((0.1+normrnd(.04,.015))*ifd(k,i));
   end
   
   rec(k,i)=ifd(k,i)-mor(k,i);
   
   tot(k,i)=tot(k,i-1)-mor(k,i);
   
   totdth(k,i)=totdth(k,i-1)+mor(k,i);
   
%    morord=mor(k,i-mem+1:i);
   morord=sort(mor(k,i-mem+1:i));
   
   if mor(k,i)<mlow && mor(k,i-1)<mlow && mor(k,i-2)<mlow && mor(k,i-3)<mlow &&... 
       mor(k,i-4)<mlow && mor(k,i-5)<mlow && mor(k,i-6)<mlow && mor(k,i-7)<mlow
       pen(k,i)=1;
   elseif mor(k,i)>mhigh && mor(k,i-1)>mhigh && mor(k,i-2)>mhigh && mor(k,i-3)>mhigh && ...
           mor(k,i-4)>mhigh && mor(k,i-5)>mhigh && mor(k,i-6)>mhigh && mor(k,i-7)>mhigh
       pen(k,i)=penmax;
   else pen(k,i)=max(1,sum(w.*((penmax-1)/(mhigh-mlow)*flip(morord)+1-(penmax-1)/(mhigh-mlow)*mlow)));
   end
   
   opttim(k,i)=1./(2*pen(k,i)^2);
    
end 
end


A=[mean(tot);mean(ifd);mean(totifd);mean(rec);mean(mor);mean(totdth);mean(pen);mean(opttim)];

%% Importing Data for Comparison

% Import USA COVID-19 dataset
T = readtable('usadata.xlsx');

% Select mortality rate from data 
dat_mor=T{:,7};
dat_mor=dat_mor(dat_mor>0);

% Select infection rate from data 
dat_inf=T{:,5};

%% Simulation Results

% Model predications

figure(1)
subplot(1,2,1)
bar([A(5,lag+1:end)' A(2,lag+1:end)'],'stacked')
lgd=legend('Mortality','Infected','interpreter','latex');
lgd.FontSize = 10;
xlabel('t','interpreter','latex','fontsize',14)
title('Per time Infection and Mortality','interpreter','latex','fontsize',14)
set(gca,'FontSize',14)
grid on
subplot(1,2,2)
hold on
plot(1:N,A(6,lag+1:end),'x-')
plot(1:N,A(3,lag+1:end),'o-')
hold off
grid on
title('Total Infection and Mortality','interpreter','latex','fontsize',14)
legend('Mortality','Infected','location','northwest','interpreter','latex')
xlabel('t','interpreter','latex')
set(gca,'FontSize',14)

figure(2)
subplot(1,2,1)
hold on
plot(1:N,A(end,lag+1:end),'s-')
hold off
grid on
xlabel('t','interpreter','latex')
legend('H_t','location','southeast','Interpreter','latex')
title('Labor Supply','interpreter','latex','fontsize',14)
set(gca,'FontSize',14)
subplot(1,2,2)
hold on
plot(1:N,sqrt(2*A(end,lag+1:end)),'x-')
hold off
grid on
xlabel('t','interpreter','latex')
legend('Y_t','location','southeast','Interpreter','latex')
title('Output','interpreter','latex','fontsize',14)
set(gca,'FontSize',14)

% Plot showing model and actual mortality rate

figure(3)
subplot(1,2,1)
bar(A(5,lag+1:end)')
lgd=legend('Mortality','interpreter','latex');
lgd.FontSize = 10;
xlabel('t','interpreter','latex','fontsize',14)
ylabel('$M_t$','interpreter','latex','fontsize',14)
title('Per time Model Mortality','interpreter','latex','fontsize',14)
set(gca,'FontSize',14)
grid on
subplot(1,2,2)
bar(dat_mor)
lgd=legend('Mortality','interpreter','latex');
lgd.FontSize = 10;
xlabel('t','interpreter','latex','fontsize',14)
ylabel('$M_t$','interpreter','latex','fontsize',14)
title('USA Daily Mortality ','interpreter','latex','fontsize',14)
set(gca,'FontSize',14)
grid on

 
% Plot showing model and actual infection rate

figure(4)
subplot(1,2,1)
bar(A(2,lag+1:end)')
lgd=legend('Infections','interpreter','latex');
lgd.FontSize = 10;
xlabel('t','interpreter','latex','fontsize',14)
ylabel('$I_t$','interpreter','latex','fontsize',14)
title('Per time Model Infections','interpreter','latex','fontsize',14)
set(gca,'FontSize',14)
grid on
subplot(1,2,2)
bar(dat_inf)
lgd=legend('Infections','interpreter','latex');
lgd.FontSize = 10;
xlabel('t','interpreter','latex','fontsize',14)
ylabel('$I_t$','interpreter','latex','fontsize',14)
title('USA Daily Infections','interpreter','latex','fontsize',14)
set(gca,'FontSize',14)
grid on

