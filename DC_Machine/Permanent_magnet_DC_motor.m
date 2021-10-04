clear all;
clc;
close all;
%% inputs to DC motor %%
Ra = 7;
Laa = 120*(10^-3);
Kv = 1.41 * 10^-2;
Wr_NL = 351.1;
Ia_NL = 0.15;
J = 1.06*10^-6;
Tl = 3.53 * 10^-3;
Va = 6;

%% Calculation for Bm frm No-Load Te and Wr:
Te_NL = Kv * Ia_NL;
Bm = Te_NL/Wr_NL;

%% Initial conditions:
Ia1(1,1) = 0;
Wr1(1,1) = 0;
Te(1,1) = 0;
h = 0.001;
t1 = 0:h:20;

for i = 1:length(t1)
    
    %% calculations for finding k1 and m1:
    % calculation for k1:
    k1(1,i) =  -(Ra/Laa)*Ia1(1,i) - (Kv/Laa)*Wr1(1,i) + (1/Laa)*Va;
    
    % Finding electrical torque(Te) for current iteration:
    Te(1,i) = Kv * Ia1(1,i);
    
    % calculation for m1:
    m1(1,i) = -(Bm/J) * Wr1(1,i) + (Te(1,i)/J) - (1/J) * Tl;
    
    %% calculations for finding k2 and m2:
    
    t2(1,i) = t1(1,i) + (h/2);
    Ia2(1,i) = Ia1(1,i) + (k1(1,i)*(h/2));
    Wr2(1,i) = Wr1(1,i) + (m1(1,i)*(h/2));
    
    % calculations for k2:
    k2(1,i) = -(Ra/Laa)*Ia2(1,i) - (Kv/Laa)*Wr2(1,i) + (1/Laa)*Va ;
    
    %calculations for m2:
    m2(1,i) = -(Bm/J) * Wr2(1,i) + (Te(1,i)/J) - (1/J) * Tl;
    
    %% calculations for finding k3 and m3:
    
    t3(1,i) = t1(1,i) + (h/2);
    Ia3(1,i) = Ia1(1,i) + (k2(1,i)*(h/2));
    Wr3(1,i) = Wr1(1,i) + (m2(1,i)*(h/2));
    
    % calculations for k3:
    k3(1,i) = -(Ra/Laa)*Ia3(1,i) - (Kv/Laa)*Wr3(1,i) + (1/Laa)*Va;
    
    %calculations for m3:
    m3(1,i) =  -(Bm/J) * Wr3(1,i) + (Te(1,i)/J) - (1/J) * Tl;
    
    %% calculations for finding k4 and m4:
    t4(1,i) = t1(1,i) +  h;
    Ia4(1,i) = Ia1(1,i) + (k3(1,i)*h);
    Wr4(1,i) = Wr1(1,i) + (m3(1,i)*h);
    
    % calculations for k4:
    k4(1,i) = -(Ra/Laa)*Ia4(1,i) - (Kv/Laa)*Wr4(1,i) + (1/Laa)*Va;
    
    %calculations for m4:
    m4(1,i) = -(Bm/J) * Wr4(1,i) + (Te(1,i)/J) - (1/J) * Tl;
    
    
    %% calculation for next value of Armature current(Ia1):
    Ia1(1,i+1) = Ia1(1,i) + (h/6) * (k1(1,i) + 2*k2(1,i) + 2*k3(1,i) + k4(1,i));
    
   
    
    %% calculations for next value of Rotor-speed(Wr1):
    Wr1(1,i+1) = Wr1(1,i) + (h/6) * (m1(1,i) + 2*m2(1,i) + 2*m3(1,i) + m4(1,i));
   
       
end

%% plotting Ia1:
figure(1)
Ia1new = Ia1(1,(1:length(Ia1)-1));
plot(t1,Ia1new,'b');

%% plotting Wr1:
figure(2)
Wr1new = Wr1(1,(1:length(Wr1)-1));
plot(t1,Wr1new,'r');



