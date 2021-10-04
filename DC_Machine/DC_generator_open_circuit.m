clear all;
clc;
close all;

%% inputs to DC generator:
Laa = 120*(10^-3);
Kv = 1.41 * 10^-2;
J = 1.06*10^-6;
Tl = 2*3.53 * 10^-3;
Bm = 6.04 * 10^-6;

%% initial conditions:
Wr1(1,1) = 0;
Va(1,1) = 0;
h = 0.001;
t1 = 0:h:20;

for i = 1:length(t1)
    
    % calculation for m1:
      m1(1,i) = -(Bm/J) * Wr1(1,i) + (1/J) * Tl;
     
    % calculation for m2:
    t2(1,i) = t1(1,i) + (h/2);
    Wr2(1,i) = Wr1(1,i) + (m1(1,i)*(h/2));
     m2(1,i) = -(Bm/J) * Wr2(1,i) + (1/J) * Tl;
    
    % calculation for m3:
    t3(1,i) = t1(1,i) + (h/2);
    Wr3(1,i) = Wr1(1,i) + (m2(1,i)*(h/2)); 
     m3(1,i) = -(Bm/J) * Wr3(1,i) + (1/J) * Tl;
     
    % calculation for m4:
    t4(1,i) = t1(1,i) +  h;
    Wr4(1,i) = Wr1(1,i) + (m3(1,i)*h);
     m4(1,i) = -(Bm/J) * Wr4(1,i) + (1/J) * Tl;
     
    %% calculation for next value of Rotor speed(Wr1):
    Wr1(1,i+1) = Wr1(1,i) + (h/6) * (m1(1,i) + 2*m2(1,i) + 2*m3(1,i) + m4(1,i));
    
    %% calculation for next value of Armature voltage(Va):
    Va(1,i+1) = Kv * Wr1(1,i+1);
end


%% plotting Wr1new:
figure(1)
Wr1new = Wr1(1,(1:length(Wr1)-1));
plot(t1,Wr1new,'r');

%% plotting Vanew:
figure(2)
Vanew = Va(1,(1:length(Va)-1));
plot(t1,Vanew,'b');

%%plotting Vanew and Wr1new:
figure(3)
plot(Wr1new,Vanew,'c');