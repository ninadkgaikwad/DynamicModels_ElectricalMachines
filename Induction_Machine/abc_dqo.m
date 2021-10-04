function [ a1 ] = abc_dqo( a,b,c,theta )
%This function takes in abc frame quantities and converts them to dqo 
%quantities   
a2=[a;b;c];
a3=(2/3)*[cos(theta) cos(theta-(2*pi/3)) cos(theta+(2*pi/3));
    sin(theta) sin(theta-(2*pi/3)) sin(theta+(2*pi/3));
    0.5 0.5 0.5];
a1=a3*a2;

end

