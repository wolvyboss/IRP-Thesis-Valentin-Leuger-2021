function [dx] = TWO_BP_Pos_Vel(t,x,Mu_ast,Rot_Speed_norm_ast)
r_norm = sqrt(x(1)^2+x(2)^2+x(3)^2);
dx = zeros(6,1);
dx(1)=x(4);%dx/dt
dx(2)=x(5);%dy/dt
dx(3)=x(6);%dz/dt
dx(4)=(-Mu_ast*x(1)/(r_norm^3))+2*Rot_Speed_norm_ast...
    *x(5)+(Rot_Speed_norm_ast^2)*x(1);%dvx/dt
dx(5)=(-Mu_ast*x(2)/(r_norm^3))-2*Rot_Speed_norm_ast...
   *x(4)+(Rot_Speed_norm_ast^2)*x(2);%dvy/dt
dx(6)=-Mu_ast*x(3)/(r_norm^3);%dvz/dt

end

