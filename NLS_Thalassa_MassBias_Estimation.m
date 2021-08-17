%%=========================================================================
%                        Cranfield University
%            NonLinear Least Square calculator for s/c orbiting around an 
%                         asteroid for the 2-Body-Problem
%             
%%=========================================================================
%                   
%

clear all 
close all
clc

%% Defining the scenario of our s/c

% The aim of this section is to define the parameters of the asteroid
Prot=5*3600;    %[s] Rotationnal Period of the asteroid
Rho_ast = 2000;                                %[kg/m^3]Density of the asteroid
R_ast =250;                                    %[m] Radius of the asteroid
Rot_Speed_norm_ast =2*pi/Prot;                 %[rad/s] Rotationnal Speed of the asteroid 
Volume_ast = 4/3*pi*R_ast^3;                   %[m^3] Volume of the asteroid
G = 6.674E-11;                                 %[m^3⋅kg^−1⋅s^−2] Gravitationnal constante
Mu_ast = Volume_ast*Rho_ast*G;                 %standard gravitational parameter of the asteroid
Ast_Rot_Axis = [0 0 1]  ;                      %Orientation of the rotation axis of the asteroid
Rot_Speed_Vect_ast = Ast_Rot_Axis*Rot_Speed_norm_ast; %rotationnal speed of %the asteroif as a vector
                              %[m] the mean deviation of the position from the measurement

%% Defining the initial condition of our s/c

%The aim of this section is to define the initial position and velocity of
%the s/c with respect to the asteroid

Alt_SC = 3*R_ast; % The alitutde of the spacecraft
R0_SC = [Alt_SC, 0, 0]; %Position Vector
V_SC = sqrt(Mu_ast/norm(Alt_SC)); %Velocity Vector in inertial ref frame for a circular orbit around the asteroid
V_SC_ECI = [0 V_SC 0]-cross(Rot_Speed_Vect_ast, R0_SC);%Vel wrt the asteroid
X0_SC_ECI = [R0_SC V_SC_ECI Mu_ast];%Initial State Vector for the nominal trajectory

%% Generate The Nominal Trajectory

%The aim of this section is to create the trajectory that the s/c will
%follow for real. The s/c will have to retrieve this trajectory with the
%use of the NLS

SimulationTime_Hours = 2*pi*sqrt((Alt_SC)^3/Mu_ast)/3600;      %[h]The amount of time we will propagate our orbit
DT = SimulationTime_Hours*3600; %[s] That same time expresed in [s]
tspan = [0,DT];                 % the peiod of time where we integrate the trajectory

% In order to integrate the trajcetory, we must use the equation of motion
% of the 2 Body Problem around a spherical asteroid. We will use a function
% that will compute the derivative of the state-vector and that we will
% then integrate through ode45. The function that provides dx is called
% 'TWO_BP_Pos_Vel'
options = odeset('AbsTol',1e-6,'RelTol',1e-9);
[Tnominal,SV_Nominal] = ode45(@(t,x)TWO_BP_Pos_Vel_MassBias(t,x,Rot_Speed_norm_ast),tspan,X0_SC_ECI,options);

%% Setting up the measurement process

%The aim of this section is to prepare the measurement process (which will
%be an imitation of what the sensors really output. The strategy is to take
%some points along the nominal trajcetory and to modify them a bit so they
%will look like outputs from sensors. The first step is done in this
%section and its about "pre-extracting" points from the nominal trajectory
%So we will extract an observation every Sample_H

Sample_H=0.5; %[h] sampling time in hour
tspan_H=[0, SimulationTime_Hours];
Obs_time_H=tspan_H(1):Sample_H:tspan_H(2);%The step between each obsbervations
tspan_sec = Obs_time_H*3600;
%Now we exract the points from the nominal trajcetory every Sample_H
options = odeset('AbsTol',1e-6,'RelTol',1e-9);
[~,SV_Nominal_Obs] = ode45(@(t,x)TWO_BP_Pos_Vel_MassBias(t,x,Rot_Speed_norm_ast),tspan_sec,X0_SC_ECI,options);
 
%% Creation of the measurements

%The aim of this section is to take the datas extracted in the previous
%section, and to alter them so they will look like a real measurement with
%an error for example

[nb_Points, ~]=size(SV_Nominal_Obs); %Number of observations that we want

%Initialisation of the matrices used
 
Qi = zeros(6,6,nb_Points); %Covariance matrix associated the the measuring process 
Phi = zeros(7,7,nb_Points); %STM
SV_Obs = zeros (nb_Points, 6);  
 
SigmaPos = 10;  
qqx = (SigmaPos/3)^2; %Error in position for the X axis
qqy = (SigmaPos/3)^2; %Error in position for the Y axis
qqz = (SigmaPos/3)^2; %Error in position for the Z axis

Q_0 = diag([qqx,qqy,qqz]);%qqvx,qqvy,qqvz]); %Cov Matrix 
cQ = diag(Q_0);
%Creation of the obsrvations, we will corrupt the points extracted 
 for i = 1:nb_Points
    SV_Obs(i,1:6) = CorruptX(SV_Nominal_Obs(i,1:6),cQ(1), cQ(2), cQ(3),0,0,0);
 end
 
 
%% Generate the first guess of the trajectory
 
%The aim of this section is to provide the s/c a first guess of the
%trajectory he's following so that he will be able to start the first loop
%of the NLS process. To do so, we need a first guess of the trajectory. Therefore
%will corrupt the nominal trajectory to obtain something close enough of the real trajectory
%that could be used as a first guess (the real traj is indeed the nominal
%one)
 %Definition of the error in the covariance matix
SigmaVel_FirstGuess = 0.1*V_SC; %Standard deviation for the velocity
SigmaPos_FirstGuess = 50;
qqPos = (SigmaPos_FirstGuess/3)^2;  %Variance for the position
qqVel = (SigmaVel_FirstGuess/3)^2;  %Variance for the velocity

 
%Creation of (x0,y0,z0) and (vx0,vy0,vz0) of the first guess
Pos_FirstGuess = CorruptX(R0_SC, qqPos, qqPos, qqPos);
Vel_FirstGuess = CorruptX(V_SC_ECI, qqVel, qqVel, qqVel);
Mu_FirstGuess = Mu_ast + normrnd(0,0.5);
SV_FirstGuess = [Pos_FirstGuess, Vel_FirstGuess, Mu_FirstGuess];%Initial state vector of 
%First guess
 
Phi0 = eye(7); %initialisation of the STM
reshape_Phi0 = reshape(Phi0, 1, 49); %Reshape Phi0 as a row vector bc its simpler 
% for all the computation
 
%Now we will create the state vector that we want to propagate, and it is a
%wise choice to include the STM in the state vector
 X0_FirstGuess = [SV_FirstGuess,reshape_Phi0];
 options = odeset('AbsTol',1e-9,'RelTol',1e-9);
 
%We can Now propagate the initial state vector to obtain the first guess of
%referecen trajectory. The X_FirstGuess vector will contain every potision
%and terms of the STM for each time in [0,tspan_sec]. The X_FirstGuess_Traj
%vector represents only the (pos,vol) vector. It is therefore used to
%plot the trajector of the first guess (cyan line on the plot)
[~,X_FirstGuess] = ode45(@(t,x)TWO_BP_STM_MassBias(t,x,Rot_Speed_norm_ast),tspan_sec,X0_FirstGuess,options);
[~,X_FirstGuess_Traj] = ode45(@(t,x)TWO_BP_Pos_Vel_MassBias(t,x,Rot_Speed_norm_ast),Tnominal,SV_FirstGuess,options);

%% NLS Loop
 
%The Aim of this section is to put the first guess into the NLS loop in
%order to create our estimate by correcting our reference trajcetory at
%each loop of the NLS until a certain condition on the residuals is
%achieved
 
%First we do a loop of the NLS to obtain our first correction
[x0new, Delx_Vec, PDelx_Mat, RRi] = NLSCameraSimple_MassBias(X_FirstGuess,SV_Obs,tspan_sec,Q_0);
 
%Now we initialise the loop
N = 10; %Number of iteration
X0_update = zeros(N,7); %the updated initial state vetor that the NLSCameraSImple will output
Residual = zeros(N,1); 
Correction_applied  = zeros(N,7);
Cov_Mat_Correction_applied =zeros(7,7,N);


X0_update(1,:)=x0new;
Residual(1)=sum(sqrt(RRi(:,1).^2+RRi(:,2).^2+RRi(:,3).^2));
%Now we iterate the NLS N times
j=2;
condition = norm(Residual(j)-Residual(j-1));
while condition>0.01 %j = 2:N
     %First we define the initial state vector
     x0_newIter(1:7) = X0_update(j-1,:); %x0_newIter is used to propagate the...
     %trajectory with ode45. It will take the value of the previous loop.
     %The idea is to first propagate the new x0 with ode45 and with the
     %trajectory obtained, obtain a correction to it with the NLSCameraSimple
      %we can propagate
     x0_newIter(8:56) = reshape(Phi0,1,49); %Reshaphe the STM into a row vector
    
     %Now we intergrate the x0_newIter to obtain the trajectory
     options = odeset('AbsTol',1e-6,'RelTol',1e-9);
     [T,X_corrected] = ode45(@(t,x)TWO_BP_STM_MassBias(t,x,Rot_Speed_norm_ast),tspan_sec,x0_newIter,options);
     
     %Now we compute the correction that the above trajectory needs to be
     %closer to the nominal trajectory (which is usually unknown)
     [x0new, Delx_Vec, PDelx_Mat, RRi] = NLSCameraSimple_MassBias(X_corrected,SV_Obs,tspan_sec,Q_0);
     Correction_applied(j,:)=Delx_Vec;
     Cov_Mat_Correction_applied(:,:,j)=PDelx_Mat;
     X0_update(j,:)=x0new;
     Residual(j) = sum(sqrt(RRi(:,1).^2+RRi(:,2).^2+RRi(:,3).^2));
     condition = norm(Residual(j)-Residual(j-1));
     j=j+1;
end

%semilogy([1:10,Residual])
%% Ploting of the different trajectory
% Defining Asteroid
Rmean=250;
Rx = Rmean; Ry = Rmean; Rz = Rmean; %Target axes
Rho = 2000; %Target Density
Prot=5*3600;    %[s] Rotationnal Period of the asteroid
Om = 2*pi/Prot; % Target rotational speed (radians/s); TBRMarco:where 4.29 rev comes from?
RAx = [0 0 1];  % Orientation of the target rotation axis
Ancillary(1) = Rx; Ancillary(2) = Ry; Ancillary(3) = Rz; Ancillary(4) = Rho; Ancillary(5) = Om;
Ancillary(6:8) = RAx;
%First we are going to plot the nominal trajectory, so the real trajectory
%the s/c will be following, we will plot this in Black

%We plot the asteroid first
run('plot_asteroid.m')

%Now the nominal trajectory
hold on 
plot3(SV_Nominal(:,1),SV_Nominal(:,2),SV_Nominal(:,3),'k') %plot full of the nominal trajectory
    plot3(SV_Nominal(1,1),SV_Nominal(1,2),SV_Nominal(1,3),'MarkerFaceColor',[1 0 0],'Marker','square',...
        'LineStyle','none',...
        'Color',[0 0 0]);%plot a marker for first and last point of the nominal traj
    plot3(SV_Nominal(end,1),SV_Nominal(end,2),SV_Nominal(end,3),'MarkerFaceColor',[1 0 0],'Marker','o',...
        'LineStyle','none',...
        'Color',[0 0 0]);

%Now we plot the observations points
for k = 1:nb_Points
    plot3(SV_Obs(k,1),SV_Obs(k,2),SV_Obs(k,3),'MarkerFaceColor',[0 1 1],'Marker','diamond',...
            'LineStyle','none',...
            'Color',[0 0 0]);

end
 
%Now we plot the first guess so we can know how far we can be from the
%nominal trajcetory in cyan
plot3(X_FirstGuess_Traj(:,1),X_FirstGuess_Traj(:,2),X_FirstGuess_Traj(:,3),'c')
%plot full trajectory of the first guess
 
%Now we will plot every guess of the loop in green
for j=1:N
     x0new=X0_update(j,:);
     [~,X_refined] = ode45(@(t,x)TWO_BP_Pos_Vel_MassBias(t,x,Rot_Speed_norm_ast),Tnominal,x0new,options);
     plot3(X_refined(:,1),X_refined(:,2),X_refined(:,3),'g') %plot every guess
end
%ax.Clipping = 'off';    % turn clipping off
