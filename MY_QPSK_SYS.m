function [syp] = MY_QPSK_SYS(I,Q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


I=(I*2)-1;
Q=(Q*2)-1;
syp = (1/2)(I+ j*Q);


end

