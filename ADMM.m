clc
clear;
[X,Y]=meshgrid(-10:0.1:10);

%QP problem specification
Xref=[1,1];
P=[1,0;0,1];
q=[-2*Xref(1),-2*Xref(2)]';
r=Xref(1)^2+Xref(2)^2;
%inequality constraint
A1=[1,-1];
b1=1;
A2=[-1,1];
b2=-1;
A3=[-1,-1];
b3=1;
A4=[1,1];
b4=-1;
A=[A1;A2;A3;A4];
b=[b1;b2;b3;b4];
%updating variables
%x_init=[-3,5]';
x_new=[-3,5]';
u_new=[0,0]';
z_new=[0,0]';
Max_Iterate=1000;
rho=1;
ABSTOL= 1e-4;
RELTOL= 1e-2;

%ploting related
history_x1=zeros(1000);
history_x2=zeros(1000);
F=P(1,1)*X.^2+P(1,2)*X.*Y+P(2,1)*X.*Y+P(2,2)*Y.^2 + q(1)*X+q(2)*Y + r;
f=zeros(1000);

Y_c1=-A1(1)/A1(2)*X(1,:)+b1;
Y_c2=-A2(1)/A2(2)*X(1,:)+b2;
Y_c3=-A3(1)/A3(2)*X(1,:)+b3;
Y_c4=-A4(1)/A4(2)*X(1,:)+b4;

figure(1);
contour(X,Y,F,50);
xlabel('x');
ylabel('y');
hold on;
plot(X(1,:),Y_c1,X(1,:),Y_c2,X(1,:),Y_c3,X(1,:),Y_c4);
axis([-10 10 -10 10]);
hold on;

for k=1:Max_Iterate
    history_x1(k)=x_new(1);
    history_x2(k)=x_new(2);
    f(k)=x_new'*P*x_new + q'*x_new + r;
    %x-update
    if k>1
        x_new=R\(R'\(rho*(z_new-u_new)-q));
    else
        R=chol(P+rho*eye(size(P,1)));
        x_new=R\(R'\(rho*(z_new-u_new)-q));
    end
    
    %z-update
    z_old=z_new;
    z_new=func_projection(x_new+u_new,A,b);
    
    %u-update
    u_new=u_new+x_new-z_new;
    
    %cost function
    
    plot(history_x1(k),history_x2(k),'r*');
    hold on;
    if(k>1)
        plot(history_x1(k-1:k),history_x2(k-1:k),'r');
        hold on;
    end
    
    %termination check
    r_norm(k)=norm(x_new-z_new);
    s_norm(k)=norm(-rho*(z_new-z_old));
    eps_pri(k) = sqrt(2)*ABSTOL + RELTOL*max(norm(x_new), norm(-z_new));
    eps_dual(k)= sqrt(2)*ABSTOL + RELTOL*norm(rho*u_new);
    if(r_norm(k)<eps_pri(k) && s_norm(k)<eps_dual(k))
        break;
    end
end
