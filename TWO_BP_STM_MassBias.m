function [dx] = TWO_BP_STM(t,x,Rot_Speed_norm_ast)
%the first 6 elements of x are the pos and vel, and the 36 remaining ones
%are the one of the STM. First we will compute the first element of dx
%(which are reltaed to the pos and vel) and then the compononent of the
%derivaitive of hte STM
r_norm = sqrt(x(1)^2+x(2)^2+x(3)^2);
dx = zeros(7,1);
dx(1)=x(4);%dx/dt
dx(2)=x(5);%dy/dt
dx(3)=x(6);%dz/dt
dx(4)=((-x(7)*x(1))/(r_norm^3))+2*Rot_Speed_norm_ast...
    *x(5)+(Rot_Speed_norm_ast^2)*x(1);%dvx/dt
dx(5)=((-x(7)*x(2))/(r_norm^3))-(2*Rot_Speed_norm_ast...
   *x(4))+(Rot_Speed_norm_ast^2)*x(2);%dvy/dt
dx(6)=((-x(7)*x(3))/(r_norm^3));%dvz/dt
dx(7)=0; %derivative w.r.t time of the gravity constant of the asteroid constant of the asteroid


%Now we will compute the derivaitive of the STM, first we need to
%extract it from the state vector
Phi = reshape(x(8:end),7,7);%We re-create the STM as a 7x7 matrix
    
%Then we will create the Jacobian Matrix for the 2BP. (cf the Extra
%Exercise 2 of the ATAD Class)
A_21_11 =(-x(7)/(r_norm^3))+((3*x(7)*x(1)^2)/(r_norm^5))+Rot_Speed_norm_ast^2;
A_21_12 = (3*x(7)*x(1)*x(2))/(r_norm^5);
A_21_13 = (3*x(7)*x(1)*x(3))/(r_norm^5);
A_21_21 = (3*x(7)*x(1)*x(2))/(r_norm^5);
A_21_22 =(-x(7)/(r_norm^3))+((3*x(7)*x(2)^2)/(r_norm^5))+Rot_Speed_norm_ast^2;
A_21_23 = (3*x(7)*x(2)*x(3))/(r_norm^5);
A_21_31 = (3*x(7)*x(1)*x(3))/(r_norm^5);
A_21_32 = (3*x(7)*x(2)*x(3))/(r_norm^5);
A_21_33 =(-x(7)/(r_norm^3))+((3*x(7)*x(3)^2)/(r_norm^5));


A_21 = [A_21_11 A_21_12 A_21_13;
    A_21_21 A_21_22 A_21_23;
    A_21_31 A_21_32 A_21_33];
    
A_22 = [0 2*Rot_Speed_norm_ast 0;
    -2*Rot_Speed_norm_ast 0 0;
    0 0 0];
A_3x1=[-x(1)/(r_norm^3);
    -x(2)/(r_norm^3);
    -x(3)/(r_norm^3)];

    

%we can now create the jacobian
A = [zeros(3) eye(3) zeros(3,1); 
    A_21 A_22 A_3x1;
    zeros(1,3) zeros(1,3) zeros(1,1)];

d_Phi = A*Phi;
dx = [dx; reshape(d_Phi,49,1)]; 
    
    
    
end




    





