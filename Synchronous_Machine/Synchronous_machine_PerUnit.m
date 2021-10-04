 %------------------------------------------------------------------------%
%START%
%------------------------------------------------------------------------%
% This program models a Synchronus Machine using Runge-Kutta 4th order Method %
%------------------------------------------------------------------------%

clear all;
clc;
close all;

%------------------------------------------------------------------------%
% Machine Constants %
%------------------------------------------------------------------------%
J=0.089;    %Machine Moment of Inertia in Kg.m^2
B=0;    %Co-efficient of Vicious Friction in S.I.
P=4;    %Number of Poles 
%------------------------------------------------------------------------%



%------------------------------------------------------------------------%
% Machine Dependent System Constants %
%------------------------------------------------------------------------%
Tm1=125;        %Constant Load Torque or Prime Mover Torque in N.m 
            %(Use when mechanical torque is constant)
Vms=311.127;        %Peak Line-Neutral Stator Voltage in V
V_f=0;        % Field Voltage in V
fs=50;         %Stator side Electrical Supply Frequency in Hz
Ws=2*pi*fs; %Stator Angular Velocity in rads/sec

%------------------------------------------------------------------------%




%% Number of turns:
%------------------------------------------------------------------------%
Ns =[];
Nkq1 =[];
Nkq2 =[];
Nfd =[];
Nkd =[];

%% List of Required Inductances:
%------------------------------------------------------------------------%
La =[];
Lb =[];
Lls =[];
Ll_kq1 =[];
Ll_kq2 =[];
Ll_fd =[];
Ll_kd =[];


Lmq = (3/2) * (La - Lb);
Lmd = (3/2) * (La + Lb);

Lskq1 = (Nkq1/Ns) * (2/3) * Lmq;
Lskq2 = (Nkq2/Ns) * (2/3) * Lmq;
Lsfd = (Nfd/Ns) * (2/3) * Lmd;
Lskd = (Nkd/Ns) * (2/3) * Lmd;
Lmkq1 = ((Nkq1/Ns)^2) * (2/3) * Lmq;
Lmkq2 = ((Nkq2/Ns)^2) * (2/3) * Lmq;
Lmfd = ((Nfd/Ns)^2) * (2/3) * Lmd;
Lmkd = ((Nkd/Ns)^2) * (2/3) * Lmd;
Lkq1kq2 = (Nkq2/Nkq1) * Lmkq1;
Lfdkd = (Nkd/Nfd) * Lmfd;

%% List of Required Reactances:
%------------------------------------------------------------------------%
Xls =[];
Xl_kq1 =[];
Xl_kq2 =[];
Xl_fd =[];
Xl_kd =[];


Xmq = [];
Xmd = [];


%% List of resistances:
%------------------------------------------------------------------------%
Rs = 
Rkq1 =
Rkq2 = 
Rfd =
Rkd = 

R_kq1 = (3/2) * ((Ns/Nkq1)^2) * Rkq1;  % r'kq1 = R_kq1
R_kq2 = (3/2) * ((Ns/Nkq2)^2) * Rkq2; % r'kq2 = R_kq2
R_fd = (3/2) * ((Ns/Nfd)^2) * Rfd; % r'fd = R_fd
R_kd = (3/2) * ((Ns/Nkd)^2) * Rkd; % r'kd = R_kd

%------------------------------------------------------------------------%
%% Initializations %
%------------------------------------------------------------------------%

% Stator Volatges in abc Frame in V (s-Stator;)%
%------------------------------------------------------------------------%
Vas(1,1)=0;
Vbs(1,1)=0;
Vcs(1,1)=0;


% Stator Voltages in qdo(Rotor) Frame in V %
%------------------------------------------------------------------------%
Vqsr(1,1)=0;
Vdsr(1,1)=0;
Vos(1,1)=0;

% Rotor Voltages in qdo(Rotor) Frame in V %
%------------------------------------------------------------------------%
V_kq1r(1,1)=0;
V_kq2r(1,1)=0;
V_fdr(1,1) = Vf ;
V_kdr(1,1)=0;


% Stator Currents in abc Frame in A (s-Stator)%
%------------------------------------------------------------------------%
Ias(1,1)=0;
Ibs(1,1)=0;
Ics(1,1)=0;

% Stator Currents in qdo Frame in A (s-Stator)%
%------------------------------------------------------------------------%
Iqsr(1,1)=0;
Idsr(1,1)=0;
Ios(1,1)=0;

% Rotor Currents in qdo(Rotor) Frame in A%
%------------------------------------------------------------------------%
I_kq1r(1,1)=0;
I_kq2r(1,1)=0;
I_fdr(1,1)=0;
I_kdr(1,1)=0;


% Stator Fluxes in qdo(Rotor) Frame (s-Stator)%
%------------------------------------------------------------------------%
Lamqsr(1,1)=0;
Lamdsr(1,1)=0;
Lamos(1,1)=0;

% Rotor Fluxes in qdo(Rotor) Frame (s-Stator)%
%------------------------------------------------------------------------%
Lam_kq1r(1,1)=0;
Lam_kq2r(1,1)=0;
Lam_fdr(1,1)=0;
Lam_kdr=0;

% Stator Fluxes Per Second in qdo(Rotor) Frame (s-Stator)%
%------------------------------------------------------------------------%
Siqsr(1,1) = Wb * Lamqsr(1,1);
Sidsr(1,1) = Wb * Lamdsr(1,1);
Sios(1,1) = Wb * Lamos(1,1);

% Rotor Fluxes Per Second in qdo(Rotor) Frame (s-Stator)%
%------------------------------------------------------------------------%
Si_kq1r(1,1) = Wb * Lam_kq1r(1,1);
Si_kq2r(1,1) = Wb * Lam_kq2r(1,1);
Si_fdr(1,1) = Wb * Lam_fdr(1,1);
Si_kdr(1,1) = Wb * Lam_kdr(1,1);

%------------------------------------------------------------------------%
% Electrical Torque(Te) in N.m, 
% Mechanical Load/Input Torque(Tm) in N.m, 
% Rotor Speed (wr-Electrical,wm-Mechanical)in rads/sec 
% Rotor Angle (ThetaR) in radians

Te(1,1)=0;
Tm(1,1)=0;% Use when Mechanical torque is varying
wr(1,1)=0;
wm(1,1)=0;
ThetaR(1,1)=0;

%------------------------------------------------------------------------%
% Loop Counter (i), Time (t) in secs, Time Step (h) in secs  %
%------------------------------------------------------------------------%
i=0;
t(1,1)=0;
h= 0.001;
%------------------------------------------------------------------------%

%------------------------------------------------------------------------%
% Starting with Machine Simulation 'For Loop' %
%------------------------------------------------------------------------%
for j=0:h:10
    
    i=i+1; % Incrementing Loop Counter
    t(1,i)=j; % Time Vector intialization
    
    % Stator Voltages in abc Frame for 'i'th count %
    %--------------------------------------------------------------------%
    Vas(1,i)=Vms*sin(Ws*t(1,i));
    Va1=Vas(1,i);
    Vbs(1,i)=Vms*sin((Ws*t(1,i))-(2*pi/3));
    Vb1=Vbs(1,i);
    Vcs(1,i)=Vms*sin((Ws*t(1,i))-(4*pi/3));
    Vc1=Vcs(1,i);
    
    % Calculation of Angles for dqo(Rotor) Transformation Matrix %
    %--------------------------------------------------------------------%
    Theta1 = Theta1 + (wr(1,i) * h); % For Stator Transformation into Rotor frame 
    
    
    % Calculation of Stator dqo Voltages for 'i'th count %
    %--------------------------------------------------------------------%
    c1=abc_dqo(Va1,Vb1,Vc1,Theta1);
    Vqsr(1,i)=c1(1,1);
    Vdsr(1,i)=c1(2,1);
    Vos(1,i)=c1(3,1);
    
    % Rotor Voltages in qdo(Rotor) Frame for 'i'th count %
    %--------------------------------------------------------------------%
    V_kq1r(1,i)=0;
    V_kq2r(1,i)=0;
    V_fdr(1,i)= Vf;
    V_kdr(1,i)=0;
       

    
    % Mechanical Load/Input Torque Function %
    %--------------------------------------------------------------------%
    %Tm(1,i)=function(t,wm);
    
    % Mechanical Load/Input Torque Constant %
    %--------------------------------------------------------------------%
    Tm(1,i)=Tm1(1,1);
    
    % Stator Flux Calculations %
    %--------------------------------------------------------------------% 
    % k = Slope for Lamqsr % 
    % m = Slope for Lamdsr %
    % n = Slope for Lamos %
    
    %% Calculations for k1, m1, n1:

k1(1,i) = Wb * (Vqsr(1,i)+(Rs*Iqsr(1,i))-(wr(1,i)*Sidsr(1,i)));
m1(1,i) = Wb * (Vdsr(1,i)+(Rs*Idsr(1,i))+(wr(1,i)*Siqsr(1,i)));
n1(1,i) = Wb * (Vos(1,i)+(Rs*Ios(1,i)));

%% Calculations for k2, m2, n2:

t2(1,i) = t(1,i) + (h/2);

Siqsr2(1,i) = Siqsr(1,i) + (k1(1,i)*(h/2));
Sidsr2(1,i) = Sidsr(1,i) + (m1(1,i)*(h/2));
Sios2(1,i) = Sios(1,i) + (n1(1,i)*(h/2));

 k2(1,i) = Wb * (Vqsr(1,i)+(Rs*Iqsr(1,i))-(wr(1,i)*Sidsr2(1,i)));
 m2(1,i) = Wb * (Vdsr(1,i)+(Rs*Idsr(1,i))+(wr(1,i)*Siqsr2(1,i)));
 n2(1,i) = Wb * (Vos(1,i)+(Rs*Ios(1,i)));
 
 
%% Calculations for k3, m3, n3:
 
t3(1,i) = t(1,i) + (h/2);
 
Siqsr3(1,i) = Siqsr2(1,i) + (k2(1,i)*(h/2));
Sidsr3(1,i) = Sidsr2(1,i) + (m2(1,i)*(h/2));
Sios3(1,i) = Sios2(1,i) + (n2(1,i)*(h/2));

 k3(1,i) = Wb * (Vqsr(1,i)+(Rs*Iqsr(1,i))-(wr(1,i)*Sidsr3(1,i)));
 m3(1,i) = Wb * (Vdsr(1,i)+(Rs*Idsr(1,i))+(wr(1,i)*Siqsr3(1,i)));
 n3(1,i) = Wb * (Vos(1,i)+(Rs*Ios(1,i)));
 
%% Calculations for k4, m4, n4:
t4(1,i) = t(1,i) +  h;

Siqsr4(1,i) = Siqsr3(1,i) + (k3(1,i)*(h/2));
Sidsr4(1,i) = Sidsr3(1,i) + (m3(1,i)*(h/2));
Sios4(1,i) = Sios3(1,i) + (n3(1,i)*(h/2));
 
  k4(1,i) = Wb * (Vqsr(1,i)+(Rs*Iqsr(1,i))-(wr(1,i)*Sidsr3(1,i)));
  m4(1,i) = Wb * (Vdsr(1,i)+(Rs*Idsr(1,i))+(wr(1,i)*Siqsr3(1,i)));
  n4(1,i) = Wb * (Vos(1,i)+(Rs*Ios(1,i)));
  
%% Calculations for the next value of Lamqsr, Lamdsr, Lamos:

Siqsr(1,i+1) = Siqsr(1,i) + (h/6) * (k1(1,i) + 2*k2(1,i) + 2*k3(1,i) + k4(1,i));

Sidsr(1,i+1) = Sidsr(1,i) + (h/6) * (m1(1,i) + 2*m2(1,i) + 2*m3(1,i) + m4(1,i));

Sios(1,i+1) = Sios(1,i) + (h/6) * (n1(1,i) + 2*n2(1,i) + 2*n3(1,i) + n4(1,i));
    %--------------------------------------------------------------------%
    
    %% Rotor Flux Calculations %
    %--------------------------------------------------------------------%
    % z = Si_kq1r%
    % x = Si_kq2r %
    % c = Si_fdr %
    % s = Si_kdr %
    %% Calculations for x1, c1, n1, s1:

z1(1,i) = Wb * (V_kq1r(1,i) - (R_kq1 * I_kq1r(1,i)));
x1(1,i) = Wb * (V_kq2r(1,i) - (R_kq2 * I_kq2r(1,i)));
c1(1,i) = Wb * (V_fdr(1,i) - (R_fd * I_fdr(1,i)));
s1(1,i) = Wb * (V_kdr(1,i) - (R_kd * I_kdr(1,i)));

%% Calculations for z2, x2, c2, s2:

t2(1,i) = t(1,i) + (h/2);

Si_kq1r2(1,i) = Si_kq1r(1,i) + (z1(1,i) * (h/2));
Si_kq2r2(1,i) = Si_kq2r(1,i) + (x1(1,i) * (h/2));
Si_fdr2(1,i) = Si_fdr(1,i) + (c1(1,i) * (h/2));
Si_kdr2(1,i) = Si_kdr(1,i) + (s1(1,i) * (h/2));

 z2(1,i) = Wb * (V_kq1r(1,i) - (R_kq1 * I_kq1r(1,i)));
 x2(1,i) = Wb * (V_kq2r(1,i) - (R_kq2 * I_kq2r(1,i)));
 c2(1,i) = Wb * (V_fdr(1,i) - (R_fd * I_fdr(1,i)));
 s2(1,i) = Wb * (V_kdr(1,i) - (R_kd * I_kdr(1,i)));
 
 
%% Calculations for z3, x3, c3, s4:
 
t3(1,i) = t(1,i) + (h/2);

Si_kq1r3(1,i) = Si_kq1r2(1,i) + (z2(1,i) * (h/2));
Si_kq2r3(1,i) = Si_kq2r2(1,i) + (x2(1,i) * (h/2));
Si_fdr3(1,i) = Si_fdr2(1,i) + (c2(1,i) * (h/2));
Si_kdr3(1,i) = Si_kdr2(1,i) + (s2(1,i) * (h/2));


 z3(1,i) = Wb * (V_kq1r(1,i) - (R_kq1 * I_kq1r(1,i)));
 x3(1,i) = Wb * (V_kq2r(1,i) - (R_kq2 * I_kq2r(1,i)));
 c3(1,i) = Wb * (V_fdr(1,i) - (R_fd * I_fdr(1,i)));
 s3(1,i) = Wb * (V_kdr(1,i) - (R_kd * I_kdr(1,i)));
 
%% Calculations for z4, x4, c4, s4:
t4(1,i) = t(1,i) +  h;

Si_kq1r4(1,i) = Si_kq1r3(1,i) + (z3(1,i) * (h/2));
Si_kq2r4(1,i) = Si_kq2r3(1,i) + (x3(1,i) * (h/2));
Si_fdr4(1,i) = Si_fdr3(1,i) + (c3(1,i) * (h/2));
Si_kdr4(1,i) = Si_kdr3(1,i) + (s3(1,i) * (h/2));

 z4(1,i) = Wb * (V_kq1r(1,i) - (R_kq1 * I_kq1r(1,i)));
 x4(1,i) = Wb * (V_kq2r(1,i) - (R_kq2 * I_kq2r(1,i)));
 c4(1,i) = Wb * (V_fdr(1,i) - (R_fd * I_fdr(1,i)));
 s4(1,i) = Wb * (V_kdr(1,i) - (R_kd * I_kdr(1,i)));
 
  
%% Calculations for the next value of Lamqr, Lamdr, Lamor:

Si_kq1r(1,i+1) = Si_kq1r(1,i) + (h/6) * (z1(1,i) + 2*z2(1,i) + 2*z3(1,i) + z4(1,i));

Si_kq2r(1,i+1) = Si_kq2r(1,i) + (h/6) * (x1(1,i) + 2*x2(1,i) + 2*x3(1,i) + x4(1,i));

Si_fdr(1,i+1) = Si_fdr(1,i) + (h/6) * (c1(1,i) + 2*c2(1,i) + 2*c3(1,i) + c4(1,i));

Si_kdr(1,i+1) = Si_kdr(1,i) + (h/6) * (s1(1,i) + 2*s2(1,i) + 2*s3(1,i) + s4(1,i));
%--------------------------------------------------------------------%  
    
    % Rotor Speed Calculations (Mechanical Degrees) %
    %--------------------------------------------------------------------%
    %S17=(1/J)*(Te(1,i)); % For Free Motoring action without B
    %S17=(1/J)*(Tm(1,i)-Te(1,i)); % For Generating action without B
    %S17=(1/J)*(Tm(1,i)-Te(1,i)-(B*wm(1,i))); % For Generating action with B
    %Above equations can be used for Motoring action when Tm is constant
    
    S17=(1/J)*(Te(1,i)-Tm(1,i)); % For Motoring action without B
    %S17=(1/J)*(Te(1,i)-Tm(1,i)-(B*wm(1,i))); % For Motoring action with B
    %Above equations can be used for Motoring action when Tm is Varying
    
    
    
    % Calculations for u1(u1 = slope of wm w.r.t. time):
     u1(1,i) = (1/J)*(Te(1,i)-Tm(1,i)); % Here for motoring without B, the formula for slope of 'wm' doesn't consist of wm
     % Formula of slope of 'wm' for motoring with B will consist of wm term
     % as seen above
    
    % Calculations for u2:
    t2(1,i) = t(1,i) + (h/2);
    
    wm2(1,i) = wm(1,i) + (u1(1,i)*(h/2));
    
     u2(1,i) = (1/J)*(Te(1,i)-Tm(1,i));
    
    % Calculations for u3:
    t3(1,i) = t(1,i) + (h/2);
    
    wm3(1,i) = wm2(1,i) + (u2(1,i)*(h/2));
    
     u3(1,i) = (1/J)*(Te(1,i)-Tm(1,i));
    
    % Calculations for u4:
    t4(1,i) = t(1,i) +  h;
    
    wm4(1,i) = wm3(1,i) + (u3(1,i)*(h/2));
    
     u4(1,i) = (1/J)*(Te(1,i)-Tm(1,i));
     
    % Calculation for next value of wm:
    wm(1,i+1) = wm(1,i) + (h/6) * (u1(1,i) + 2*u2(1,i) + 2*u3(1,i) + u4(1,i));

    
    
     
    %--------------------------------------------------------------------%
    
    % Rotor Speed Calculations (Electrical Degrees) %
    %--------------------------------------------------------------------%
    wr(1,i+1)=(P/2)*(wm(1,i+1));
    
    
    
   
    %-----------------------------------------------------------------------%
    %% Linear Solver for calculating Stator and Rotor Currents in qdo Frame:
    
    A = [-(Xls+Xmq) 0 0 Xmq Xmq 0 0; 0 -(Xls+Xmd) 0 0 0 Xmd Xmd; 0 0 -Xls 0 0 0 0; -Xmq 0 0 (Xl_kq1+Xmq) Xmq 0 0; -Xmq 0 0 Xmq (Xl_kq2+Xmq) 0 0; 0 -Xmd 0 0 0 (XLl_fd+Xmd) Xmd; 0 -Xmd 0 0 0 Xmd (XLl_kd+Xmd)];
    B = [Siqsr(1,i+1); Sidsr(1,i+1); Sios(1,i+1); Si_kq1r(1,i+1); Si_kq2r(1,i+1); Si_fdr(1,i+1); Si_kdr(1,i+1)];
    X = linsolve(A, B);
    
    
    %-----------------------------------------------------------------------%
    % Calculation of Angles for dqo Transformation Matrix %
    %--------------------------------------------------------------------%
    Theta3 = Theta1 + (wr(1,i+1) * h); % For Stator Transformation 
   
    
    % Calculation of Electrical Torque (Te)in N.m at 'i+1'th count%
    %--------------------------------------------------------------------%
    Te(1,i+1)=(3/2) * (P/2) * ((Lamdsr(1,i+1)*Iqsr(1,i+1)) - (Lamqsr(1,i+1)*Idsr(1,i+1)));
    
    % Calculation of Angles for dqo Transformation Matrix %
    %--------------------------------------------------------------------%
    Theta3=Wf*(t(1,i)+h); % For Stator Transformation 
    %Wrr1=Wf-wr(1,i+1); % Relative Angular Velocity between dqo frame and Rotor
    Theta4=Theta3-ThetaR(1,i+1); % For Rotor Transformation
    
    % Calculation of Stator abc Currents for 'i+1'th count %
    %--------------------------------------------------------------------%
    c3=dqo_abc(Iqsr,Idsr,Ios,Theta3);
    Ias(1,i+1)=c3(1,1);
    Ibs(1,i+1)=c3(2,1);
    Ics(1,i+1)=c3(3,1);
       

    
    
end
%------------------------------------------------------------------------%
% New Time Vector For Calculated Quantities (one element incremented) %
%------------------------------------------------------------------------%
t1=t;
t1(1,i+1)=t1(1,i)+h;

%------------------------------------------------------------------------%
% New Tm Vector For making elements equal to calculated quantities (one element incremented)  %
%------------------------------------------------------------------------%
Tm2=Tm;
Tm2(1,i+1)=Tm2(1,i);

%------------------------------------------------------------------------%
% Plotting Stator abc Voltages %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t,Vas)
title('Stator abc Voltages');
ylabel('Vas (V)');
subplot(3,1,2)
plot(t,Vbs)
ylabel('Vbs (V)');
subplot(3,1,3)
plot(t,Vcs)
ylabel('Vcs (V)');
xlabel('Time (sec)');

%------------------------------------------------------------------------%
% Plotting Rotor abc Voltages %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t,Var)
title('Rotor abc Voltages');
ylabel('Var (V)');
subplot(3,1,2)
plot(t,Vbr)
ylabel('Vbr (V)');
subplot(3,1,3)
plot(t,Vcr)
ylabel('Vcr (V)');
xlabel('Time (sec)');

%------------------------------------------------------------------------%
% Plotting Stator qdo Voltages %
%------------------------------------------------------------------------% 
figure
subplot(3,1,1)
plot(t,Vqs)
title('Stator qdo Voltages');
ylabel('Vqs (V)');
subplot(3,1,2)
plot(t,Vds)
ylabel('Vds (V)');
subplot(3,1,3)
plot(t,Vos)
ylabel('Vos (V)');
xlabel('Time (sec)');

%------------------------------------------------------------------------%
% Plotting Rotor qdo Voltages %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t,Vqr)
title('Rotor dqo Voltages');
ylabel('Vqr (V)');
subplot(3,1,2)
plot(t,Vdr)
ylabel('Vdr (V)');
subplot(3,1,3)
plot(t,Vor)
ylabel('Vor (V)');
xlabel('Time (sec)');
%------------------------------------------------------------------------%
% Plotting Stator qdo Fluxes %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t1,Lamqs)
title('Stator qdo Fluxes');
ylabel('Lamqs');
subplot(3,1,2)
plot(t1,Lamds)
ylabel('Lamds');
subplot(3,1,3)
plot(t1,Lamos)
ylabel('Lamos (V)');
xlabel('Time (sec)');

%------------------------------------------------------------------------%
% Plotting Rotor qdo Fluxes %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t1,Lamqr)
title('Rotor qdo Fluxes');
ylabel('Lamqr');
subplot(3,1,2)
plot(t1,Lamdr)
ylabel('Lamdr');
subplot(3,1,3)
plot(t1,Lamor)
ylabel('Lamor (V)');
xlabel('Time (sec)');


%------------------------------------------------------------------------%
% Plotting Stator qdo Currents %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t1,Iqs)
title('Stator qdo Currents');
ylabel('Iqs (A)');
subplot(3,1,2)
plot(t1,Ids)
ylabel('Ids (A)');
subplot(3,1,3)
plot(t1,Ios)
ylabel('Ios (A)');
xlabel('Time (sec)');    
   
%------------------------------------------------------------------------%
% Plotting Rotor dqo Currents %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t1,Iqr)
title('Rotor dqo Currents');
ylabel('Iqr (A)');
subplot(3,1,2)
plot(t1,Idr)
ylabel('Idr (A)');
subplot(3,1,3)
plot(t1,Ior)
ylabel('Ior (A)');
xlabel('Time (sec)');

%------------------------------------------------------------------------%
% Plotting Stator abc Currents %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t1,Ias)
title('Stator abc Currents');
ylabel('Ias (A)');
subplot(3,1,2)
plot(t1,Ibs)
ylabel('Ibs (A)');
subplot(3,1,3)
plot(t1,Ics)
ylabel('Ics (A)');
xlabel('Time (sec)');  

%------------------------------------------------------------------------%
% Plotting Rotor abc Currents %
%------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(t1,Iar)
title('Rotor abc Currents');
ylabel('Iar (A)');
subplot(3,1,2)
plot(t1,Ibr)
ylabel('Ibr (A)');
subplot(3,1,3)
plot(t1,Icr)
ylabel('Icr (A)');
xlabel('Time (sec)');    

%------------------------------------------------------------------------%
% Plotting Rotor Speed in Radian/sec and Revs/min %
%------------------------------------------------------------------------%
wm1= (30/pi)*(wm) ; % Calculating Rotor Speed in Rev/min 

figure
subplot(2,1,1)
plot(t1,wm)
title('Rotor Speed');
ylabel('wm (Rads/sec)');
subplot(2,1,2)
plot(t1,wm1)
ylabel('wm (Revs/min)');
xlabel('Time (sec)');

%------------------------------------------------------------------------%
% Plotting Mechanical and Electrical Torque %
%------------------------------------------------------------------------%
figure
subplot(2,1,1)
plot(t1,Te)
title('Elcetrical and Mechanical Torques');
ylabel('Te(N.m)');
subplot(2,1,2)
plot(t,Tm)
ylabel('Tm (N.m)');
xlabel('Time (sec)');

%------------------------------------------------------------------------%
% Plotting Mechanical and Electrical Torque vs Rotor Speed (rpm) %
%------------------------------------------------------------------------%
figure
subplot(2,1,1)
plot(wm1,Te)
title('Electrical and Mechanical Torques vs Rotor Speed');
ylabel('Te(N.m)');
subplot(2,1,2)
plot(wm1,Tm2)
ylabel('Tm (N.m)');
xlabel('Rotor Speed (rpm)');

%------------------------------------------------------------------------%
% END %
%------------------------------------------------------------------------%*