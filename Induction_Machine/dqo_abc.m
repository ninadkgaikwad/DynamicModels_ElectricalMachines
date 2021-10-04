function [ a1 ] = dqo_abc( q,d,o, theta )

%This function takes in abc frame quantities and converts them to dqo 
%quantities   
a2=[q;d;o];
a3=[cos(theta) sin(theta) 1;
    cos(theta-(2*pi/3)) sin(theta-(2*pi/3)) 1;
    cos(theta+(2*pi/3)) sin(theta+(2*pi/3)) 1];
a1=a3*a2;


end

