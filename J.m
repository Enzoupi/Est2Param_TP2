function [Y] = J(F)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
n=length(F)+1;
h=1/n;
J=h/2*sum(direct(F)-Uobs); % A faire

end

