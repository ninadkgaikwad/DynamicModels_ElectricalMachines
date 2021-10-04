    S11=Vqs(1,i)-(Rs*Iqs(1,i))-(Wf*Lamds(1,i));
    Lamqs(1,i+1)=Lamqs(1,i)+(h*S11);
    
    S12=Vds(1,i)-(Rs*Ids(1,i))+(Wf*Lamqs(1,i));
    Lamds(1,i+1)=Lamds(1,i)+(h*S12);
    
    S13=Vos(1,i)-(Rs*Ios(1,i));
    Lamos(1,i+1)=Lamos(1,i)+(h*S13);
 
    
    
       % Stator Flux Calculations %
    %--------------------------------------------------------------------% 
    
    % k = Slope for Lamqs % 
    % m = Slope for Lamds %
    % n = Slope for Lamos %
    
%% initial conditions:
     h = 0.001;
     t1 = 0:h:40;
     
%% Calculations for k1, m1, n1:

k1 = Vqs(1,i)-(Rs*Iqs(1,i))-(Wf*Lamds(1,i));
m1 = Vds(1,i)-(Rs*Ids(1,i))+(Wf*Lamqs(1,i));
n1 = Vos(1,i)-(Rs*Ios(1,i));

%% Calculations for k2, m2, n2:

t2(1,i) = t1(1,i) + (h/2);

Lamqs2(1,i) = Lamqs(1,i) + (k1(1,i)*(h/2));
Lamds2(1,i) = Lamds(1,i) + (m1(1,i)*(h/2));
Lamos2(1,i) = Lamos(1,i) + (n1(1,i)*(h/2));

 k2 = Vqs(1,i)-(Rs*Iqs(1,i))-(Wf*Lamds2(1,i));
 m2 = Vds(1,i)-(Rs*Ids(1,i))+(Wf*Lamqs2(1,i));
 n2 = Vos(1,i)-(Rs*Ios(1,i));
 
 
%% Calculations for k3, m3, n3:
 
t3(1,i) = t1(1,i) + (h/2);
 
Lamqs3(1,i) = Lamqs2(1,i) + (k2(1,i)*(h/2));
Lamds3(1,i) = Lamds2(1,i) + (m2(1,i)*(h/2));
Lamos3(1,i) = Lamos2(1,i) + (n2(1,i)*(h/2));

 k3 = Vqs(1,i)-(Rs*Iqs(1,i))-(Wf*Lamds3(1,i));
 m3 = Vds(1,i)-(Rs*Ids(1,i))+(Wf*Lamqs3(1,i));
 n3 = Vos(1,i)-(Rs*Ios(1,i));
 
%% Calculations for k4, m4, n4:
t4(1,i) = t1(1,i) +  h;

Lamqs4(1,i) = Lamqs3(1,i) + (k3(1,i)*(h/2));
Lamds4(1,i) = Lamds3(1,i) + (m3(1,i)*(h/2));
Lamos4(1,i) = Lamos3(1,i) + (n3(1,i)*(h/2));
 
  k4 = Vqs(1,i)-(Rs*Iqs(1,i))-(Wf*Lamds4(1,i));
  m4 = Vds(1,i)-(Rs*Ids(1,i))+(Wf*Lamqs4(1,i));
  n4 = Vos(1,i)-(Rs*Ios(1,i));
  
%% Calculations for the next value of Lamqs, Lamds, Lamos:

Lamqs(1,i+1) = Lamqs(1,i) + (h/6) * (k1(1,i) + 2*k2(1,i) + 2*k3(1,i) + k4(1,i));

Lamds(1,i+1) = Lamds(1,i) + (h/6) * (m1(1,i) + 2*m2(1,i) + 2*m3(1,i) + m4(1,i));

Lamos(1,i+1) = Lamos(1,i) + (h/6) * (n1(1,i) + 2*n2(1,i) + 2*n3(1,i) + n4(1,i));



 % Rotor Flux Calculations %
    %--------------------------------------------------------------------%
 % z = Lamqr %
 % x = Lamdr %
 % c = Lamor %
 
 %% initial conditions:
     h = 0.001;
     t1 = 0:h:40;
     
%% Calculations for x1, c1, n1:

z1 = Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Lamdr(1,i));
x1 = Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Lamqr(1,i));
c1 = Vor(1,i)-(Rr*Ior(1,i));

%% Calculations for z2, x2, c2:

t2(1,i) = t1(1,i) + (h/2);

Lamqr2(1,i) = Lamqr(1,i) + (z1(1,i)*(h/2));
Lamdr2(1,i) = Lamdr(1,i) + (x1(1,i)*(h/2));
Lamor2(1,i) = Lamor(1,i) + (c1(1,i)*(h/2));

 z2 = Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Lamdr2(1,i));
 x2 = Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Lamqr2(1,i));
 c2 = Vor(1,i)-(Rr*Ior(1,i));
 
 
%% Calculations for z3, x3, c3:
 
t3(1,i) = t1(1,i) + (h/2);
 
Lamqr3(1,i) = Lamqr2(1,i) + (z2(1,i)*(h/2));
Lamdr3(1,i) = Lamdr2(1,i) + (x2(1,i)*(h/2));
Lamor3(1,i) = Lamor2(1,i) + (c2(1,i)*(h/2));

 z3 = Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Lamdr3(1,i));
 x3 = Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Lamqr3(1,i));
 c3 = Vor(1,i)-(Rr*Ior(1,i));
 
%% Calculations for z4, x4, c4:
t4(1,i) = t1(1,i) +  h;

Lamqr4(1,i) = Lamqr3(1,i) + (k3(1,i)*(h/2));
Lamdr4(1,i) = Lamdr3(1,i) + (m3(1,i)*(h/2));
Lamor4(1,i) = Lamor3(1,i) + (n3(1,i)*(h/2));
 
  z4 = Vqr(1,i)-(Rr*Iqr(1,i))-(Wrr*Lamdr4(1,i))
  x4 = Vdr(1,i)-(Rr*Idr(1,i))+(Wrr*Lamqr4(1,i));
  c4 = Vor(1,i)-(Rr*Ior(1,i));
  
%% Calculations for the next value of Lamqr, Lamdr, Lamor:

Lamqr(1,i+1) = Lamqr(1,i) + (h/6) * (z1(1,i) + 2*z2(1,i) + 2*z3(1,i) + z4(1,i));

Lamdr(1,i+1) = Lamdr(1,i) + (h/6) * (x1(1,i) + 2*x2(1,i) + 2*x3(1,i) + x4(1,i));

Lamor(1,i+1) = Lamor(1,i) + (h/6) * (c1(1,i) + 2*c2(1,i) + 2*c3(1,i) + c4(1,i));
