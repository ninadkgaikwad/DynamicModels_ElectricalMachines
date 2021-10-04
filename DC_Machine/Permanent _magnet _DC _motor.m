%% inputs to DC motor %%
Ra = 7;
Laa = 120*(10^-3);
Kv = 1.41 * 10^-2;
Wr_NL = 351.1;
Ia_NL = 0.15;
J = 1.06*10^-6;
Tl = 3.53 * 10^-3;
Va = 6;

%% Calculation of Bm %%
Te_NL = Kv * Ia_NL;
Bm = Te_NL/Wr_NL;

%% Initial conditions %%
Ia1(1,1) = 0;
Wr1(1,1) = 0;
h = 0.1;
t1 = 0:h:20;

%% calculation of Ia %%
for i = 1:length(t1)
    
    % calculation for k1:
    k1(1,i) =  -(Ra/Laa)*Ia1(1,i)-(1/Laa)*Va;
    
    %calculations for k2:
    t2(1,i) = t1(1,i) + (h/2);
    Ia2(1,i) = Ia1(1,i) + (k1(1,i)*(h/2));
    k2(1,i) = -(Ra/Laa)*Ia2(1,i)-(1/Laa)*Va ;
    
    %calculations for k3:
    t3(1,i) = t1(1,i) + (h/2);
    Ia3(1,i) = Ia1(1,i) + (k2(1,i)*(h/2));
    k3(1,i) =  -(Ra/Laa)*Ia3(1,i)-(1/Laa)*Va;
    
    %calculations for k4:
    t4(1,i) = t1(1,i) +  h;
    Ia4(1,i) = Ia1(1,i) + (k3(1,i)*h);
    k4(1,i) =  -(Ra/Laa)*Ia4(1,i)-(1/Laa)*Va;
    
    %calculation for new y1:
    Ia1(1,i+1) = Ia1(1,i) + (h/6) * (k1(1,i) + 2*k2(1,i) + 2*k3(1,i) + k4(1,i));
    
end

Ia1new = Ia1(1,(1:length(Ia1)-1));
plot(t1,Ia1new,'r');
hold on
%start ode45 calculations
interval = [0 20];
y_initial = 0;
[x,y] = ode45(@myode45func_1, interval, y_initial);
plot(x,y,'c*');
hold off;

