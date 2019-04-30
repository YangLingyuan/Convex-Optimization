function [ICR_next, solved]=CO_toolbox(P,q,A,b)
options = optimset('MaxIter',500);
[x,fval,exitflag] = quadprog(P,q,A,b);
ICR_next=x';
solved = exitflag;
