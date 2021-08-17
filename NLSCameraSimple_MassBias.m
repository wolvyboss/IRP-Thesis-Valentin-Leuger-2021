    function [x0new, Delx_Vec, PDelx_Mat, RRi] = NLSCameraSimple(Xref, Xobs, T, Qi)
%Non-linear Least Squares
%Inputs:
%Xref0 is the matrix of reference states
%Xobs the matrix of state observations
%T is the matrix containing the timesteps
%Qi is the matrix containing the covariance of the

%Outputs: 
%x0new is the initial state vector of the new trajectory
%Delx_Vec is the correction applied to the input trajectory
%Pdelx_Mat is the covariance matrix of the correction applied
%RRi is the residual of the corrected trajectory

%% Matrices allocation:
% !!! Given a matrix S, Si represents that matrix for the timestep i.
%Hi is the derivative w.r.t. the state of the observation relation
Hi = zeros(3,7, length(T));
%Ti represents the observation matrix Hi*Phi
Ti = zeros(3,7, length(T));
%State Transition Matrix
PhiT = zeros(7,7, length(T));
%Residuals
RRi = zeros(length(T),3);
%Running sum of
SUM1_Mat = zeros(7,7,length(T));
%Running sum of
SUM2_Vec = zeros(7,length(T));
%Tensor containing the inverse of the covariance matrix
QiInv = zeros(3,3,length(T));
%The tensor containing the covariances is considered given

%This process is the one explained in page 74 of Modern Orbit Determination
for tt = 1:length(T)
    PhiT(:,:,tt) = reshape(Xref(tt,8:end),7,7);%we put the stm as a matrix
    QiInv(:,:,tt) = inv(Qi);
    RRi(tt,:) = Xobs(tt,1:3)' - Xref(tt,1:3)'; %Obtain the residuals
    Hi(1:3,1:3, tt) = eye(3);
    Ti(:,:,tt) = Hi(:,:, tt)*PhiT(:,:,tt);
     if tt > 1 %TBR: why you don't use the first observable?
        SUM1_Mat(:,:, tt) = SUM1_Mat(:,:,tt-1) + (transpose(Ti(:,:,tt)))*QiInv(:,:,tt)*Ti(:,:,tt);
        SUM2_Vec(:,tt) = SUM2_Vec(:, tt-1) + (transpose(Ti(:,:,tt)))*QiInv(:,:,tt)*RRi(tt,:)';
     elseif tt==1
        SUM1_Mat(:,:, tt) = (transpose(Ti(:,:,tt)))*QiInv(:,:,tt)*Ti(:,:,tt);
        SUM2_Vec(:,tt) = (transpose(Ti(:,:,tt)))*QiInv(:,:,tt)*RRi(tt,:)'; 
     end
   % end
end
% RRi;
PDelx_Mat = inv(SUM1_Mat(:,:,end)); %The Covariance Matrix of the estimate
Delx_Vec = PDelx_Mat*SUM2_Vec(:,end); %The correction to apply at the ref traj
x0new = Xref(1,1:7); 
x0new(1:7) = x0new(1:7) + Delx_Vec';%The new initial vector with the correction
end

