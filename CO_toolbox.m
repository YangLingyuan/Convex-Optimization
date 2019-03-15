P = [1 0; 0 1]; 
Xref=[2,2];
q=[-2*Xref(1),-2*Xref(2)];
r=Xref(1)^2+Xref(2)^2;
%f = [-2; -6];
A1=[1,-1];
b1=2;
A2=[-1,1];
b2=2;
A3=[-1,-1];
b3=2;
A4=[1,1];
b4=2;
A=[A1;A2;A3;A4];
b=[b1;b2;b3;b4];


[X,Y]=meshgrid(-10:0.1:10);
F=P(1,1)*X.^2+P(1,2)*X.*Y+P(2,1)*X.*Y+P(2,2)*Y.^2 + q(1)*X+q(2)*Y + r;
Y_c1=-A1(1)/A1(2)*X(1,:)-sign(A1(2))*b1;
Y_c2=-A2(1)/A2(2)*X(1,:)-sign(A2(2))*b2;
Y_c3=-A3(1)/A3(2)*X(1,:)-sign(A3(2))*b3;
Y_c4=-A4(1)/A4(2)*X(1,:)-sign(A4(2))*b4;

figure(1);
contour(X,Y,F,50);
xlabel('x');
ylabel('y');
hold on;
plot(X(1,:),Y_c1,X(1,:),Y_c2,X(1,:),Y_c3,X(1,:),Y_c4);
axis([-10 10 -10 10]);
hold on;
[x,fval,exitflag] = quadprog(P,q,A,b);
plot(x(1),x(2),'r*');