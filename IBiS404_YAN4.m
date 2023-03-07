function  IBiS404_YAN4
clear
clc
clf
%%  Model
% This model includes Yan formation and degradation, the phosphorylation of Yan to Yanp via MAPK, and the
% degradation and dephosphorylation of YanP. Additionally, the model contains fast conversion of
% YanP to YanP-U (ubiquitinated Yan) along with the fast degradation of
% YanP-U. There is positive feedback where YanP represses the phosphatase that dephosphorylates it.
% The output in this case is the total level of Yan - all three isoforms.

% Concentrations of ubiquitin ligase and deubiquitinating enzyme are
% assumed to be constant.
%
% Variables:
%   Y               = non-phosphorylated Yan
%   Yp              = phosphorylated Yan
%   YpU             = ubiquitinated/phosphorylated Yan
%   Ytot            = total Yan (Y+Yp+YpU).
% Model parameters:
%   M               = Level of activated MAPK
%   Vy              = Production rate of non-phosphorylated Yan
%   d               = Degradation rate constant of non-phosphorylated Yan.
%   dp              = Degradation rate constant of phosphorylated YanP.
%   km              = Phosphorylation rate constant by activated MAPK
%   Vm              = Maximal phosphorylation rate of Yan by activated MAPK = km*M (Michaelis Menten Vmax).
%   Km              = Michaelis Menten Km constant for Yan phosphorylation.
%   Vp              = Maximal dephosphorylation rate of YanP (Michaelis Menten Vmax).
%   Kp              = Michaelis Menten Km constant for YanP dephosphorylation.
%   Kf              = Michaelis Menten Km constant for positive feedback.
%   Ku              = Michaelis Menten Km constant for  YanP ubiquitination.
%   Ke              = Michaelis Menten Km constant for YanP-U deubiquitination.
%   Ve              = Maximal deubiquitination rate of YanP-U (Michaelis Menten Vmax).
%   ke              = Rate constant for deubiquitination
%   ku              = Rate constant for ubiquitination
%   du              = Degradation rate constant of ubiquitinated YanP

close all; warning off MATLAB:divideByZero
clockstart = clock;

%% The input parameters:
Mvec  = [1e-4 5e-4 1e-3 5e-3 1e-2:1e-2:0.1 0.12:0.02:0.5 0.6:0.1:1.8 1.9:1e-2:2.6 2.7:0.1:4 4.2:0.2:5.8 5.9:0.05:6.3]; % MAPK vector

disp('STARTING THE MODELING')

[d, Km, Vp, Kp, dp, km, Kf, n, Ku, Ke, du, ku, Vy_vec] =  PositiveFeedback_parameters;

for i = 1:length(Vy_vec)
    Vy = Vy_vec(i);

    [Y, Yp, YpU, Ytot]  = DynamicSS3 (Mvec,Vy,d,dp,km,Km,Kp,Vp,Kf,n,Ku,Ke,du,ku);
    Ytot_mat(:,i) = Ytot;
end

PositiveFeedback_SS_figures_new(Vy,Vy_vec,d, Mvec, Ytot_mat);

disp('FINISHED MODELING')



clockend = clock; TotalRunTime = clockend-clockstart
return

%%
function [d, Km, Vp, Kp, dp, km, Kf, n, Ku, Ke, du, ku, Vy_vec] =  PositiveFeedback_parameters
d = 5e-5;
Km= 500;
Vp= 5e4 ; Kp= 500 ;
dp= 0.1;
km= 1;
Kf= 0.004  ;
n = 6;
Ku = 10;
Ke = 500;
du = 1000;
ku = 40;
%Vy = 1e-2;
Vy_vec = [1E-4 1E-3 1E-2 1E-1 1 10];
%ku_vec = [100 500 1000 5000 10000 50000];
%ku_vec = [0 20 40 60 80 100];
%ku_vec = 60*ones(1,6);
return
%%
function EQ = Equation_PF_Dynamics(t,var,Vy,d,dp,Km,Kp,Vm,Vp,Kf,n,Ku,Ke,Ve,du,ku)


Y  = var(1);
Yp = var(2);
YpU = var(3);

%%Original Code
EQ(1) = Vy - d*Y - Vm/(1+Km/Y) + (Vp/(1+Kp/Yp)) * (Kf^n/(Kf^n+Yp^n))    ;
EQ(2) = -1*dp*Yp + Vm/(1+Km/Y) - (Vp/(1+Kp/Yp)) * (Kf^n/(Kf^n+Yp^n)) - (ku*Yp)/(1+Ku/Yp) + Ve/(1+Ke/YpU)  ;
EQ(3) = -1*du*YpU - Ve/(1+Ke/YpU) + (ku*Yp)/(1+Ku/Yp);

%%Edit
% EQ(1) = Vy - d*Y - Vm/(1+Km/Y) + (Vp/Kp)*Yp * (Kf^n/(Kf^n+Yp^n))    ;
% EQ(2) = -1*dp*Yp + Vm/(1+Km/Y) - (Vp/Kp)*Yp * (Kf^n/(Kf^n+Yp^n)) - (ku*Yp)/(1+Ku/Yp) + Ve/(1+Ke/YpU)  ;
% EQ(3) = -1*du*YpU - Ve/(1+Ke/YpU) + (ku*Yp)/(1+Ku/Yp);

EQ = transpose(EQ);
return
%%
function  [Y, Yp, YpU, Ytot] = DynamicSS3 (Mvec,Vy,d,dp,km,Km,Kp,Vp,Kf,n,Ku,Ke,du,ku)
clock1=clock;

y0  = Vy/d;

for j=1:length(Mvec)
    Mtot= Mvec(j);
    Vm  = km*Mtot;
    Ve  = 10;

    options = odeset('RelTol',1e-5) ;
    [~,VAR] = ode15s(@Equation_PF_Dynamics, 0:10:2e4, [y0 0 0], options,Vy,d,dp,Km,Kp,Vm,Vp,Kf,n,Ku,Ke,Ve,du,ku);  % solver for EQ1, EQ2, and EQ3 from t=0 to t=20000, starting at y0 amount of Yan
    Y_vec  = VAR(:,1); %vector of Y levels at each time step
    Yp_vec = VAR(:,2); %vector of Yp levels at each time step
    YpU_vec = VAR(:,3); %vector of YpU levels at each time step
    Y(j) = Y_vec(end); % take only the final Y level in the vector to be the steady-state value for Y at the current MAPK concentration
    Yp(j) = Yp_vec(end); %take only the final Yp level in the vector to be the steady-state value for Yp at the current MAPK concentration
    YpU(j) = YpU_vec(end); %take only the final YpU level in the vector to be the steady-state value for YpU at the current MAPK concentration
    Ytot(j) = Y(j)+Yp(j)+YpU(j); %the total Yan level is the sum of the levels of all of its forms. Sum these values to get the total amount of Yan at the end of modeling time at the current MAPK concentration
end

clock2=clock;
Runtime_PositiveDynamics3 = clock2-clock1
return
%%
function PositiveFeedback_SS_figures_new(Vy,Vy_vec,d, Mvec, Ytot_mat)
y0    = Vy/d;
Mnorm = Mvec/max(Mvec) ; % vector of normalized MAPK concentrations
dgaus = sqrt(log(1./Mnorm))*2;       % Assume that MAPK activation is proportional to local EGF signaling. Assume EGF is secreted from a local source and
% diffuses following a Gaussian distribution. So the distribution of
% activated MAPK is also Gaussian.

%Figure
h=figure; set(0,'DefaultFigureWindowStyle','docked')%set(h,'units','normalized','position',[0 0.0741 1 0.838]);
H=plot(dgaus,Mnorm,'k',dgaus,Ytot_mat(:,1)/(Vy_vec(:,1)/d),'b-',dgaus,Ytot_mat(:,2)/(Vy_vec(:,2)/d),'r-',dgaus,Ytot_mat(:,3)/(Vy_vec(:,3)/d),'y-',dgaus,Ytot_mat(:,4)/(Vy_vec(:,4)/d),'g-',dgaus,Ytot_mat(:,5)/(Vy_vec(:,5)/d),'c-',dgaus,Ytot_mat(:,6)/(Vy_vec(:,6)/d),'m-'); %plots MAPK concentration and total Yan concentration (normalized)
set(H,'linewidth',5,'markersize',12);
h=gca;  set(h,'fontsize',30);
xlabel('CELL DIAMETERS FROM EGF SOURCE','color','k','fontsize',15,'fontweight','bold');
ylabel('FRACTIONAL YAN LEVEL / MAPK ACTIVITY','color','k','fontsize',15,'fontweight','bold');
title ('Positive Feedback + Ubiquitination','color','k','fontsize',15,'fontweight','bold');
% 
%labels =num2str(ku_vec.','Total Yan, ku = %d');

%legend('Activated MAPK','Total Yan, ku = 0','Total Yan, ku = 20','Total Yan, ku = 40','Total Yan, ku = 60','Total Yan, ku = 80','Total Yan, ku = 100') ;
legend('Activated MAPK','Total Yan, Vy = 1e-4','Total Yan, Vy = 1e-3','Total Yan, Vy = 1e-2','Total Yan, Vy =1e-1','Total Yan, Vy = 10','Total Yan, Vy = 100');
% legend(labels,'location','best')
%legend('Activated MAPK', labels)
legend_handle = legend ; set(legend_handle,'fontsize',12,'fontweight','bold');

return