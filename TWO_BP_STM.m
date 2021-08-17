function [dx] = TWO_BP_STM(t,x,Mu_ast,Rot_Speed_norm_ast)
%the first 6 elements of x are the pos and vel, and the 36 remaining ones
%are the one of the STM. First we will compute the first element of dx
%(which are reltaed to the pos and vel) and then the compononent of the
%derivaitive of hte STM
r_norm = sqrt(x(1)^2+x(2)^2+x(3)^2);
dx = zeros(6,1);
dx(1)=x(4);%dx/dt
dx(2)=x(5);%dy/dt
dx(3)=x(6);%dz/dt
dx(4)=(-Mu_ast*x(1)/(r_norm^3))+2*Rot_Speed_norm_ast...
    *x(5)+(Rot_Speed_norm_ast^2)*x(1);%dvx/dt
dx(5)=(-Mu_ast*x(2)/(r_norm)^3)-(2*Rot_Speed_norm_ast...
   *x(4))+(Rot_Speed_norm_ast^2)*x(2);%dvy/dt
dx(6)=-Mu_ast*x(3)/(r_norm)^3;%dvz/dt

%Now we will compute the derivaitive of the STM, first we need to
%extract it from the state vector
Phi = reshape(x(7:end),6,6);%We re-create the STM as a matrix
    
%Then we will create the Jacobian Matrix for the 2BP. (cf the Extra
%Exercise 2 of the ATAD Class)
A_21_11 =(-Mu_ast/r_norm^3)+(3*Mu_ast*x(1)^2/(r_norm^5))+Rot_Speed_norm_ast^2;
A_21_12 = 3*Mu_ast*x(1)*x(2)/(r_norm)^5;
A_21_13 = 3*Mu_ast*x(1)*x(3)/r_norm^5;
A_21_21 = 3*Mu_ast*x(1)*x(2)/r_norm^5;
A_21_22 =(-Mu_ast/r_norm^3)+(3*Mu_ast*x(2)^2/(r_norm^5))+Rot_Speed_norm_ast^2;
A_21_23 = 3*Mu_ast*x(2)*x(3)/r_norm^5;
A_21_31 = 3*Mu_ast*x(1)*x(3)/r_norm^5;
A_21_32 = 3*Mu_ast*x(2)*x(3)/r_norm^5;
A_21_33 = (-Mu_ast/r_norm^3)+(3*Mu_ast*x(3)^2/r_norm^5);

A_21 = [A_21_11 A_21_12 A_21_13;
    A_21_21 A_21_22 A_21_23;
    A_21_31 A_21_32 A_21_33];
    
A_22 = [0 2*Rot_Speed_norm_ast 0;
    -2*Rot_Speed_norm_ast 0 0;
    0 0 0];

%we can now create the jacobian
A = [zeros(3) eye(3); 
    A_21 A_22];
d_Phi = A*Phi;
dx = [dx; reshape(d_Phi,36,1)]; 
    
    
    
end




    



