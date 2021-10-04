 %------------------------------------------------------------------------%
%START%
%------------------------------------------------------------------------%
% This program models a Doubly - Fed Induction Machine using Runge-Kutta 4th order Method %
%------------------------------------------------------------------------%

clear all;
clc;

%------------------------------------------------------------------------%
% Machine Constants %
%------------------------------------------------------------------------%
Rs=0.435;  %Stator Resistance in Ohms/Phase
Rr=0.816;  %Rotor Resistance referred to Stator in Ohms/Phase
Ls=0.0024; %Stator Self-Inductance in Henry/Phase 
Lr=0.0024;   %Rotor Self-Inductance referred to Stator in Henry/Phase 
Lms=0.005596;   %Stator Magnetizing-Inductance in Henry/Phase  
J=0.089;    %Machine Moment of Inertia in Kg.m^2
B=0;    %Co-efficient of Vicious Friction in S.I.
P=4;    %Number of Poles 
%------------------------------------------------------------------------%

%------------------------------------------------------------------------%
% Calculated Machine Constants %
%------------------------------------------------------------------------%
Lm=(3/2)*Lms;
ls=Ls+Lm;
lr=Lr+Lm;
le=((ls*lr)/Lm)-Lm;

%------------------------------------------------------------------------%
% Machine Reactances:
%------------------------------------------------------------------------%
Xls = [];
X_lr = [];
Xm = [];

%------------------------------------------------------------------------%
% Machine Dependent System Constants %
%------------------------------------------------------------------------%
Tm1=125;        %Constant Load Torque or Prime Mover Torque in N.m 
            %(Use when mechanical torque is constant)
Vms=311.127;        %Peak Line-Neutral Stator Voltage in V
Vmr=0;        %Peak Line-Neutral Rotor Voltage in V
fs=50;         %Stator side Electrical Supply Frequency in Hz
Ws=2*pi*fs; %Stator Angular Velocity in rads/sec
fr=0;         %Rotor side Electrical Supply Frequency in Hz
Wr=2*pi*fr; %Rotor Angular Velocity in rads/sec
%------------------------------------------------------------------------%

%------------------------------------------------------------------------%
% qdo Reference Frame rotation Frequency in HZ %
%------------------------------------------------------------------------%
f=0;
Wf=2*pi*f; %qdo Frame Angular Velocity in rads/sec
%------------------------------------------------------------------------%

%------------------------------------------------------------------------%
% Initializations %
%------------------------------------------------------------------------%

% Stator and Rotor Volatges in abc Frame in V (s-Stator;r-Rotor)%
%------------------------------------------------------------------------%
Vas(1,1)=0;
Vbs(1,1)=0;
Vcs(1,1)=0;
Var(1,1)=0;
Vbr(1,1)=0;
Vcr(1,1)=0;

% Stator and Rotor Base Volatge(RMS Voltage) in abc Frame in V:
%------------------------------------------------------------------------%
Vb_abc = [];

% Stator and Rotor Volatges in qdo Frame in V (s-Stator;r-Rotor)%
%------------------------------------------------------------------------%
Vqs(1,1)=0;
Vds(1,1)=0;
Vos(1,1)=0;
Vqr(1,1)=0;
Vdr(1,1)=0;
Vor(1,1)=0;

% Stator and Rotor (RMS voltage) Base Volatge in qdo Frame in V :
%------------------------------------------------------------------------%
Vb_qdo = sqrt(2) * Vq_abc;

% Stator and Rotor Currents in abc Frame in A (s-Stator;r-Rotor)%
%------------------------------------------------------------------------%
Ias(1,1)=0;
Ibs(1,1)=0;
Ics(1,1)=0;
Iar(1,1)=0;
Ibr(1,1)=0;
Icr(1,1)=0;

% Stator and Rotor Base current in abc Frame in A:
%------------------------------------------------------------------------%
Ib_abc = [];

% Stator and Rotor Currents in qdo Frame in A (s-Stator;r-Rotor)%
%------------------------------------------------------------------------%
Iqs(1,1)=0;
Ids(1,1)=0;
Ios(1,1)=0;
Iqr(1,1)=0;
Idr(1,1)=0;
Ior(1,1)=0;

% Stator and Rotor Base Current(RMS current) in qdo Frame in V:
%------------------------------------------------------------------------%
Ib_qdo = sqrt(2) * Ib_abc;

% Base Impedance in qdo Frame:
%------------------------------------------------------------------------%
Zb_qdo = (Vb_qdo/Ib_qdo);

% Base electrical Angular velocity in rad/sec:
%------------------------------------------------------------------------%
Wb = [];

% Stator and Rotor Fluxes in qdo Frame (s-Stator;r-Rotor)%
%------------------------------------------------------------------------%
Lamqs(1,1)=0;
Lamds(1,1)=0;
Lamos(1,1)=0;
Lamqr(1,1)=0;
Lamdr(1,1)=0;
Lamor(1,1)=0;

% Stator and Rotor Flux Per Second in qdo Frame (s-Stator;r-Rotor)%
%------------------------------------------------------------------------%
Siqs(1,1) = Lamqs(1,1)*Wb
Sids(1,1) = Lamds(1,1)*Wb;
Sios(1,1) = Lamos(1,1)*Wb;
Siqr(1,1) = Lamqr(1,1)*Wb;
Sidr(1,1) = Lamdr(1,1)*Wb;
Sior(1,1) = Lamor(1,1)*Wb;

% Electrical Torque(Te) in N.m, 
% Mechanical Load/Input Torque(Tm) in N.m, 
% Rotor Speed (wr-Electrical,wm-Mechanical)in rads/sec 
% Rotor Angle (ThetaR) in radians
%------------------------------------------------------------------------%
Te(1,1)=0;
Tm(1,1)=0;% Use when Mechanical torque is varying
wr(1,1)=0;
wm(1,1)=0;
ThetaR(1,1)=0;

% Loop Counter (i), Time (t) in secs, Time Step (h) in secs  %
%------------------------------------------------------------------------%
i=0;
t(1,1)=0;
h= 0.001;
%------------------------------------------------------------------------%

%------------------------------------------------------------------------%
% Starting with Machine Simulation For Loop %
%------------------------------------------------------------------------%
for j=0:h:10
    
    i=i+1; % Incrementing Loop Counter
    t(1,i)=j; % Time Vector intialization
    
    % Stator abc Voltages for 'i'th count %
    %--------------------------------------------------------------------%
    Vas(1,i)=Vms*sin(Ws*t(1,i));
    Va1=Vas(1,i);
    Vbs(1,i)=Vms*sin((Ws*t(1,i))-(2*pi/3));
    Vb1=Vbs(1,i);
    Vcs(1,i)=Vms*sin((Ws*t(1,i))-(4*pi/3));
    Vc1=Vcs(1,i);
    
    % Rotor abc Voltages for Squirrel Cage Machine for 'i'th count %
    %--------------------------------------------------------------------%
    Var(1,i)=0;
    Va2=Var(1,i);
    Vbr(1,i)=0;
    Vb2=Vbr(1,i);
    Vcr(1,i)=0;
    Vc2=Vcr(1,i);
    
    % Rotor abc Voltages for Wound Rotor Machine for 'i'th count%
    %--------------------------------------------------------------------%
%     Var(1,i)=Vms*sin(Wr*t(1,i));
%     Va2=Var(1,i);
%     Vbr(1,i)=Vms*sin((Wr*t(1,i))-(2*pi/3));
%     Vb2=Vbr(1,i);
%     Vcr(1,i)=Vms*sin((Wr*t(1,i))-(4*pi/3);
%     Vc2=Vcr(1,i);

    % Calculation of Angles for dqo Transformation Matrix %
    %--------------------------------------------------------------------%
    Theta1=Wf*t(1,i); % For Stator Transformation 
    Wrr=Wf-wr(1,i); % Relative Angular Velocity between dqo frame and Rotor
    Theta2=Theta1-ThetaR(1,i); % For Rotor Transformation
    
    % Calculation of Stator dqo Voltages for 'i'th count %
    %--------------------------------------------------------------------%
    c1=abc_dqo(Va1,Vb1,Vc1,Theta1);
    Vqs(1,i)=c1(1,1);
    Vds(1,i)=c1(2,1);
    Vos(1,i)=c1(3,1);
       
    % Calculation of Rotor dqo Voltages for 'i'th count  %
    %--------------------------------------------------------------------%
    c2=abc_dqo(Va2,Vb2,Vc2,Theta2);
    Vqr(1,i)=c2(1,1);
    Vdr(1,i)=c2(2,1);
    Vor(1,i)=c2(3,1);
    
    % Mechanical Load/Input Torque Function %
    %--------------------------------------------------------------------%
    %Tm(1,i)=function(t,wm);
    
    % Mechanical Load/Input Torque Constant %
    %--------------------------------------------------------------------%
    Tm(1,i)=Tm1(1,1);
    
    %--------------------------------------------------------------------%
    % Runge Kutta Method Calculations  %
    %--------------------------------------------------------------------%
    
    % Stator Flux Per Second Calculations %
    %--------------------------------------------------------------------% 
    % k = Slope for Siqs % 
    % m = Slope for Sids %
    % n = Slope for Sios %
    
    %% Calculations for k1, m1, n1:

k1(1,i) = Wb * (Vqs(1,i)-(Rs*Iqs(1,i))-((Wf/Wb)*Sids(1,i)));
m1(1,i) = Wb * (Vds(1,i)-(Rs*Ids(1,i))+((Wf/Wb)*Siqs(1,i)));
n1(1,i) = Wb * (Vos(1,i)-(Rs*Ios(1,i)));

%% Calculations for k2, m2, n2:

t2(1,i) = t(1,i) + (h/2);

Siqs2(1,i) = Siqs(1,i) + (k1(1,i)*(h/2));
Sids2(1,i) = Sids(1,i) + (m1(1,i)*(h/2));
Sios2(1,i) = Sios(1,i) + (n1(1,i)*(h/2));

 k2(1,i) =  Wb * (Vqs(1,i)-(Rs*Iqs(1,i))-((Wf/Wb)*Sids2(1,i)));
 m2(1,i) =  Wb * (Vds(1,i)-(Rs*Ids(1,i))+((Wf/Wb)*Siqs2(1,i)));
 n2(1,i) =  Wb * (Vos(1,i)-(Rs*Ios(1,i)));
 
 
%% Calculations for k3, m3, n3:
 
t3(1,i) = t(1,i) + (h/2);
 
Siqs3(1,i) = Siqs2(1,i) + (k2(1,i)*(h/2));
Sids3(1,i) = Sids2(1,i) + (m2(1,i)*(h/2));
Sios3(1,i) = Sios2(1,i) + (n2(1,i)*(h/2));


 k3(1,i) = Wb * (Vqs(1,i)-(Rs*Iqs(1,i))-((Wf/Wb)*Sids3(1,i)));
 m3(1,i) = Wb * (Vds(1,i)-(Rs*Ids(1,i))+((Wf/Wb)*Siqs3(1,i)));
 n3(1,i) = Wb * (Vos(1,i)-(Rs*Ios(1,i)));
 
%% Calculations for k4, m4, n4:
t4(1,i) = t(1,i) +  h;

Siqs4(1,i) = Siqs3(1,i) + (k3(1,i)*(h/2));
Sids4(1,i) = Sids3(1,i) + (m3(1,i)*(h/2));
Sios4(1,i) = Sios3(1,i) + (n3(1,i)*(h/2));

  k4(1,i) = Wb * (Vqs(1,i)-(Rs*Iqs(1,i))-((Wf/Wb)*Sids4(1,i)));
  m4(1,i) = Wb * (Vds(1,i)-(Rs*Ids(1,i))+((Wf/Wb)*Siqs4(1,i)));
  n4(1,i) = Wb * (Vos(1,i)-(Rs*Ios(1,i)));
 
  
%% Calculations for the next value of Siqs, Sids, Sios:

Siqs(1,i+1) = Siqs(1,i) + (h/6) * (k1(1,i) + 2*k2(1,i) + 2*k3(1,i) + k4(1,i));

Sids(1,i+1) = Sids(1,i) + (h/6) * (m1(1,i) + 2*m2(1,i) + 2*m3(1,i) + m4(1,i));

Sios(1,i+1) = Sios(1,i) + (h/6) * (n1(1,i) + 2*n2(1,i) + 2*n3(1,i) + n4(1,i));
    %--------------------------------------------------------------------%
    
    %% Rotor Flux Per Second Calculations %
    %--------------------------------------------------------------------%
    % z = Siqr%
    % x = Sidr %
    % c = Sior %
    %% Calculations for x1, c1, n1:

z1(1,i) = Wb * (Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Sidr(1,i)));
x1(1,i) = Wb * (Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Siqr(1,i)));
c1(1,i) = Wb * (Vor(1,i)-(Rr*Ior(1,i)));
%% Calculations for z2, x2, c2:

t2(1,i) = t(1,i) + (h/2);

Siqr2(1,i) = Siqr(1,i) + (z1(1,i)*(h/2));
Sidr2(1,i) = Sidr(1,i) + (x1(1,i)*(h/2));
Sior2(1,i) = Sior(1,i) + (c1(1,i)*(h/2));

 z2(1,i) = Wb * (Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Sidr2(1,i)));
 x2(1,i) = Wb * (Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Siqr2(1,i)));
 c2(1,i) = Wb * (Vor(1,i)-(Rr*Ior(1,i)));
 
 
%% Calculations for z3, x3, c3:
 
t3(1,i) = t(1,i) + (h/2);
 
Siqr3(1,i) = Siqr2(1,i) + (z2(1,i)*(h/2));
Sidr3(1,i) = Sidr2(1,i) + (x2(1,i)*(h/2));
Sior3(1,i) = Sior2(1,i) + (c2(1,i)*(h/2));

 z3(1,i) = Wb * (Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Sidr3(1,i)));
 x3(1,i) = Wb * (Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Siqr3(1,i)));
 c3(1,i) = Wb * (Vor(1,i)-(Rr*Ior(1,i)));
 
%% Calculations for z4, x4, c4:
t4(1,i) = t(1,i) +  h;

Siqr4(1,i) = Siqr3(1,i) + (k3(1,i)*(h/2));
Sidr4(1,i) = Sidr3(1,i) + (m3(1,i)*(h/2));
Sior4(1,i) = Sior3(1,i) + (n3(1,i)*(h/2));
 
  z4(1,i) = Wb * (Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Sidr3(1,i)));
  x4(1,i) = Wb * (Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Siqr3(1,i)));
  c4(1,i) = Wb * (Vor(1,i)-(Rr*Ior(1,i)));
  
%% Calculations for the next value of Siqr, Sidr, Sior:

Siqr(1,i+1) = Siqr(1,i) + (h/6) * (z1(1,i) + 2*z2(1,i) + 2*z3(1,i) + z4(1,i));

Sidr(1,i+1) = Sidr(1,i) + (h/6) * (x1(1,i) + 2*x2(1,i) + 2*x3(1,i) + x4(1,i));

Sior(1,i+1) = Sior(1,i) + (h/6) * (c1(1,i) + 2*c2(1,i) + 2*c3(1,i) + c4(1,i));
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
    
    
    %--------------------------------------------------------------------%
    
    % ThetaR Calculation %
    %--------------------------------------------------------------------%
    S18=wr(1,i);
   
    % Calculation for y1(y = slope of ThetaR):
    y1(1,i) = wr(1,i);
    
    % Calculation for y2:
    t2(1,i) = t(1,i) + (h/2);
    
    ThetaR2(1,i) = ThetaR(1,i) + (y1(1,i)*(h/2));
    
     y2(1,i) = wr(1,i);
     
    % Calculation for y3:
    t3(1,i) = t(1,i) + (h/2);
    
    ThetaR3(1,i) = ThetaR2(1,i) + (y2(1,i)*(h/2));
    
     y3(1,i) = wr(1,i);
     
    % Calculation for y4:
     t4(1,i) = t(1,i) +  h;
     
     ThetaR4(1,i) = ThetaR3(1,i) + (y3(1,i)*(h/2));
     
      y4(1,i) = wr(1,i);
      
    % Calculation for the next value of ThetaR:
    
    ThetaR(1,i+1) = ThetaR(1,i) + (h/6) * (y1(1,i) + 2*y2(1,i) + 2*y3(1,i) + y4(1,i));
    
    
    %--------------------------------------------------------------------%
     % Calculation of Currents using linsolve%
    %--------------------------------------------------------------------%
    A = [(Xls+xXm) 0 0 Xm 0 0; 0 (Xls+Xm) 0 0 Xm 0; 0 0 Xls 0 0 0; Xm 0 0 (X_lr+Xm) 0 0; 0 Xm 0 0 (X_lr+Xm) 0; 0 0 0 0 0 X_lr];
    
    B = [Siqs(1,i+1); Sids(1,i+1); Sios(1,i+1); Siqr(1,i+1); Sidr(1,i+1); Sior(1,i+1)];
    
    X = [Iqs; Ids; Ios; Iqr; Idr; Ior];
    
    X = linsolve(A,B);
    
    
    % Calculation of Electrical Torque (Te)in N.m at 'i+1'th count%
    %--------------------------------------------------------------------%
    Te(1,i+1)=(3/4)*P*Lm*((Iqs(1,i+1)*Idr(1,i+1))-(Ids(1,i+1)*Iqr(1,i+1)));
    
    % Calculation of Angles for dqo Transformation Matrix %
    %--------------------------------------------------------------------%
    Theta3=Wf*(t(1,i)+h); % For Stator Transformation 
    %Wrr1=Wf-wr(1,i+1); % Relative Angular Velocity between dqo frame and Rotor
    Theta4=Theta3-ThetaR(1,i+1); % For Rotor Transformation
    
    % Calculation of Stator abc Currents for 'i+1'th count %
    %--------------------------------------------------------------------%
    c3=dqo_abc(Iq1,Id1,Io1,Theta3);
    Ias(1,i+1)=c3(1,1);
    Ibs(1,i+1)=c3(2,1);
    Ics(1,i+1)=c3(3,1);
       
    % Calculation of Rotor abc Currents for 'i+1'th count  %
    %--------------------------------------------------------------------%
    c4=dqo_abc(Iq2,Id2,Io2,Theta4);
    Iar(1,i+1)=c4(1,1);
    Ibr(1,i+1)=c4(2,1);
    Icr(1,i+1)=c4(3,1);
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
plot(t1,Siqr)
title('Rotor qdo Fluxes');
ylabel('Lamqr');
subplot(3,1,2)
plot(t1,Sidr)
ylabel('Lamdr');
subplot(3,1,3)
plot(t1,Sior)
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
%------------------------------------------------------------------------%




    
    



