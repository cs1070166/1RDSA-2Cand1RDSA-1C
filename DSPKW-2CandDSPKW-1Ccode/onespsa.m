%J.C. Spall, Jan. 2000
%This code implements SPSA with constraints for theta to lie in 
%a specified hypercube (i.e., component-wise constraints).  Allows for multiple replications   
%for purposes of statistical evaluation based on knowledge of true (noise-free) loss value
%(set replications=1 if user only wants one run).

% Parameters:
% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss
% numSimulation -> this is the simulation budget that impacts the number of 2SPSA iterations
% replications -> number of independent simulations
% theta_0 -> initial point
%
function onespsa(p, sigma, type, numSimulations, replications, theta_0)
global oneSPSA_Results
% alpha=1;
alpha = 0.602;  % should be such that 2*alpha-2*gamma>1, for the gamma below alpha> 0.601)
gamma =0.101;
a=1;
c=1.15;
% c=1.9;      				%chosen by standard guidelines
A=1000;
%A=50;
errtheta=0;
lossfinal=0;            %variable for cum. loss values
theta_lo=-2.048*ones(p,1);   %lower bounds on theta  
theta_hi=2.047*ones(p,1);    %upper bounds on theta 

lossfinaleval='loss_myexample';  %choice of loss function for final perf. evaluations (noise-free)                            % evaluation
loss='loss_myexample_noise';     %loss function used in algorithm operations (with noise)
rand('seed',31415297)
randn('seed',3111113)

Ltheta0 = feval(lossfinaleval, p, theta_0, type);
if(Ltheta0==0) 
    Ltheta0=1; 
end

thetaStar = getOptima(p, type);

% lossesAllReplications = zeros(1,p);
% nmseAllReplications = zeros(1,p);

 lossesAllReplications = zeros(1,replications);
 nmseAllReplications = zeros(1,replications);

mseTheta0=(theta_0-thetaStar)'*(theta_0-thetaStar);
if(mseTheta0==0) 
    mseTheta0=1; 
end

% outer loop for replications
for i=1:replications
  theta=theta_0;
  
% inner loop runs 1SPSA for numIterations (each iter=2 function evals)
  for k=0:numSimulations/2-1
    ak = a/(k+1+A)^alpha;
    ck = c/(k+1)^gamma;
    % ck=0.01;
    delta = 2*round(rand(p,1))-1;
    thetaplus = theta + ck*delta;
    thetaminus = theta - ck*delta;
    yplus=feval(loss, p, thetaplus, sigma, type);
    yminus=feval(loss, p, thetaminus, sigma, type);
    ghat = (yplus - yminus)./(2*ck*delta);
    
    theta=theta-ak*ghat;
  
    % Project theta onto a bounded set, component-wise
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
    oneSPSA_Results(k+1,i)=(theta-thetaStar)'*(theta-thetaStar)/mseTheta0;
    
  end  
  
  lossvalue=feval(lossfinaleval, p, theta, type);
  lossfinal=lossfinal+lossvalue;
  errtheta=errtheta+(theta-thetaStar)'*(theta-thetaStar); 
  lossesAllReplications(1, i) = lossvalue/Ltheta0;
  nmseAllReplications(1, i) = (theta-thetaStar)'*(theta-thetaStar)/mseTheta0;
end
    

% Display results: normalized loss and mean square error

% str = sprintf('Normalized loss: %e +- %e, Normalised MSE: %e +- %e', lossfinal/replications/Ltheta0, std(lossesAllReplications)/(replications^.5), errtheta/replications/mseTheta0, std(nmseAllReplications)/(replications^.5));
str = sprintf('Normalized loss: %e +- %e, Normalised MSE: %e +- %e', lossfinal/replications/Ltheta0, std(lossesAllReplications), errtheta/replications/mseTheta0, std(nmseAllReplications));
disp(str);
disp(mat2str(theta,4));
