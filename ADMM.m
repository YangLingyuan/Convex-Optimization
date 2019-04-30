clc
clear;
load A2;
ICR_ref(1)=-20;
ICR_ref(2)=20;
XP=-20:0.1:20;

r=ICR_ref(1)^2+ICR_ref(2)^2;
hx1 = 0.263;
hx2 = 0.263; 
hx3 = -0.263;
hx4 = -0.263;

hy1 = -0.263;
hy2 = 0.263;
hy3 = 0.263;
hy4 = -0.263;

hxi =[hx1, hx2, hx3, hx4];
hyi =[hy1, hy2, hy3, hy4];
h_i=[hxi(1),hyi(1);hxi(2),hyi(2);hxi(3),hyi(3);hxi(4),hyi(4)];


% A(:,1)=A(:,1)./100;


figure(1);
% contour(X,Y,F,50);
% xlabel('x');
% ylabel('y');
hold on;
 for i=1:8
    if(A(i,2)>0)
        YP=-A(i,1)*XP + b(i);
    else
        YP=A(i,1)*XP - b(i);
    end
%     %             YP=-cot(Beta(i))*XP + cot(Beta(i))*h_i(i,1)+h_i(i,2);
     plot(XP,YP);
     hold on;
 end
hold on;
scatter(h_i(1,1),h_i(1,2));
hold on;
scatter(h_i(2,1),h_i(2,2));
hold on;
scatter(h_i(3,1),h_i(3,2));
hold on;
scatter(h_i(4,1),h_i(4,2));
hold on;
scatter(ICR_ref(1),ICR_ref(2),"filled");
hold on;

P = [1,0;0,1];
q = [-ICR_ref(1),-ICR_ref(2)]';


x_new=ICR_ref';
u_new=zeros(size(A,1),1);
z_new=zeros(size(A,1),1);

Max_Iterate=500;
rho=1;
ABSTOL= sqrt(2)*1e-4;
RELTOL= 1e-2;



figure(1);

% N = -(P + rho*(A'*A));
% N2 = N\(rho*A');
% N1 = N\(q - rho*A'*b);

for k=1:Max_Iterate
    
    %x-update
    x_new=-inv(P+rho*A'*A)*(q+rho*A'  * ( z_new + u_new - b ) );
%     x_new=N1 + N2*(z_new + u_new);
    %z-update
    z_old=z_new;
    z_new = max(0, -A*x_new - u_new + b);
    %z_new=func_projection(x_new+u_new,A,b);
    
    %u-update
    u_new=u_new+(A*x_new - b + z_new);
    
    %cost function
    %f(k+1)=x_new'*P*x_new + q'*x_new + r;
%     history_x1(k+1)=x_new(1);
%     history_x2(k+1)=x_new(2);
%     
%     plot(history_x1(k+1),history_x2(k+1),'r*');
%     hold on;
%     plot(history_x1(k:k+1),history_x2(k:k+1),'r');
%     hold on;
    
    %termination check
    primal_residual=norm(A*x_new - b + z_new);
    dual_residual=norm(-rho*(z_new-z_old));
    eps_pri= ABSTOL + RELTOL*max(norm(x_new), norm(-z_new));
    eps_dual= ABSTOL + RELTOL*norm(rho*u_new);
    if(primal_residual<eps_pri && dual_residual<eps_dual)
        break;
    end
end
scatter(x_new(1),x_new(2));
