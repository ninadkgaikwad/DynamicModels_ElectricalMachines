clear all;
clc;
close all;

%% inputs to DC generator:
Ra = 7;
Laa = 120*(10^-3);
Kv = 1.41 * 10^-2;
J = 1.06*10^-6;
Tl = 100 * 3.53 * 10^-3;
Bm = 6.04 * 10^-6;
Rl = 1;

%% Output voltage across DC generator:
Va = 6;

%% Initial conditions:
Ia1(1,1) = 0;
Wr1(1,1) = 0;

h = 0.001;
t1 = 0:h:40;


for i = 1:length(t1);
    %% Calculation for Armature voltage(Va):
   %% Va(1,i) = Ia1(1,i) * Rl; %%
    
    %% Calculation for output power(Pa1):
    Pa1(1,i) = Va * Ia1(1,i);
    
    %% Calculations for m1 and k1:
    m1(1,i) = -(Bm/J) * Wr1(1,i) + (Kv/J)*Ia1(1,i) + (1/J) * Tl;
    k1(1,i) =  -(Ra/Laa)*Ia1(1,i) - (Kv/Laa)*Wr1(1,i) + Va*(1/Laa);
    
    %% calculations for m2 and k2:
    t2(1,i) = t1(1,i) + (h/2);
    Wr2(1,i) = Wr1(1,i) + (m1(1,i)*(h/2));
    Ia2(1,i) = Ia1(1,i) + (k1(1,i)*(h/2));
    
     m2(1,i) = -(Bm/J) * Wr2(1,i) + (Kv/J)*Ia2(1,i) + (1/J) * Tl;
    
     k2(1,i) =  -(Ra/Laa)*Ia2(1,i) - (Kv/Laa)*Wr2(1,i) + Va*(1/Laa);
    
    %% calculations for m3 and k3:
    t3(1,i) = t1(1,i) + (h/2);
    Wr3(1,i) = Wr1(1,i) + (m2(1,i)*(h/2));
    Ia3(1,i) = Ia1(1,i) + (k2(1,i)*(h/2));
    
     m3(1,i) = -(Bm/J) * Wr3(1,i) + (Kv/J)*Ia3(1,i) + (1/J) * Tl;
     
     k3(1,i) =  -(Ra/Laa)*Ia3(1,i) - (Kv/Laa)*Wr3(1,i) + Va*(1/Laa);
     
    %% calculations for m4 and k4:
    t4(1,i) = t1(1,i) +  h;
    Wr4(1,i) = Wr1(1,i) + (m3(1,i)*h);
    Ia4(1,i) = Ia1(1,i) + (k3(1,i)*h);
    
    m4(1,i) = -(Bm/J) * Wr4(1,i) + (Kv/J)*Ia4(1,i) + (1/J) * Tl;
    
    k4(1,i) =  -(Ra/Laa)*Ia4(1,i) - (Kv/Laa)*Wr4(1,i) + Va*(1/Laa);
    
    %% Calculations for next value of rotor speed and armature current:
    Wr1(1,i+1) = Wr1(1,i) + (h/6) * (m1(1,i) + 2*m2(1,i) + 2*m3(1,i) + m4(1,i));
    
    Ia1(1,i+1) = Ia1(1,i) + (h/6) * (k1(1,i) + 2*k2(1,i) + 2*k3(1,i) + k4(1,i));
    
    
end

%% plotting Rotor speed(Wr1new):
    figure(1)
    Wr1new = Wr1(1,(1:length(Wr1)-1));
    plot(t1,Wr1new,'r');
    
%% plotting Armature current(Ia1new):
    figure(2)
    Ia1new = (Ia1(1,(1:length(Ia1)-1)));
    plot(t1,Ia1new,'b');

%% plotting Ia1new aganst Wr1new:
    figure(3)
    plot(Wr1new,Ia1new,'c');

%% plotting output Power(Pa1) w.r.t time t1:
figure(4)
plot(t1,Pa1,'g');


