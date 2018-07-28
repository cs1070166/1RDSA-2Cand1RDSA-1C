p=10;
sigma=0.01;
type=2; % quadratic loss is 1
replications = 100;
theta_0=ones(p,1);
numSimulations = 10000;

global oneSPSA_Results
global Hadamard_Results
global Circulant_Results

oneSPSA_Results = zeros(numSimulations/2,replications);
Hadamard_Results= zeros(numSimulations/2,replications);
Circulant_Results = zeros(numSimulations/2,replications);


str=sprintf('sigma: %f dim: %d sim: %d type: %d numSimulations: %d replications: %d, theta_0: %s', sigma, p, 2, type, numSimulations, replications, mat2str(theta_0));
disp(str);

disp('---------------1SPSA---------------');
onespsa(p, sigma, type, numSimulations, replications, theta_0);

disp('---------------Hadamardpert---------------');
hadamardpert(p, sigma, type, numSimulations, replications, theta_0);
 

disp('---------------Circulantpert---------------');
circulantpert(p, sigma, type, numSimulations, replications, theta_0);


oneSPSA_mean=mean(oneSPSA_Results,2);
Hadamard_mean=mean(Hadamard_Results,2);
Circulant_mean=mean(Circulant_Results,2);

oneSPSA_std=std(oneSPSA_Results,0,2);
Hadamard_std=std(Hadamard_Results,0,2);
Circulant_std=std(Circulant_Results,0,2);
 
plot(log10(oneSPSA_mean),'r-');
hold on;
plot(log10(Hadamard_mean),'b-');
hold on;
plot(log10(Circulant_mean),'g-');
hold off;
xlabel('Number of Iterations','fontweight','bold','fontsize',20);
ylabel('log10(Normalised Mean Square Error)','fontweight','bold','fontsize',20);
title('Comparision of the 1RDSA-C with 1SPSA and 1SPSA-H for 2 simulaton methods','fontweight','bold','fontsize',20);
legend('1SPSA-2R','1SPSA-2H','1RDSA-2C');


% 
% 
% 
% global oneSPSA_Results_1sim
% global Hadamard_Results_1sim
% global Circulant_Results_1sim
% 
% oneSPSA_Results_1sim = zeros(numSimulations,replications);
% Hadamard_Results_1sim = zeros(numSimulations,replications);
% Circulant_Results_1sim = zeros(numSimulations,replications);
% 
% str=sprintf('sigma: %f dim: %d sim: %d type: %d numSimulations: %d replications: %d, theta_0: %s', sigma, p, 1, type, numSimulations, replications, mat2str(theta_0));
% disp(str);
% 
% 
% disp('---------------1SPSA_1s---------------');
% onespsa_1s(p, sigma, type, numSimulations, replications, theta_0);
% 
% disp('---------------Hadamardpert_1s---------------');
% hadamardpert_1s(p, sigma, type, numSimulations, replications, theta_0);
%  
% 
% disp('---------------Circulantpert_1s---------------');
% circulantpert_1s(p, sigma, type, numSimulations, replications, theta_0);
% 
% 
% oneSPSA_mean_1sim=mean(oneSPSA_Results_1sim,2);
% Hadamard_mean_1sim=mean(Hadamard_Results_1sim,2);
% Circulant_mean_1sim=mean(Circulant_Results_1sim,2);
% 
% oneSPSA_std_1sim=std(oneSPSA_Results_1sim,0,2);
% Hadamard_std_1sim=std(Hadamard_Results_1sim,0,2);
% Circulant_std_1sim=std(Circulant_Results_1sim,0,2);
% 
% 
% %  I=1:length(oneSPSA_mean_1sim);
% %  jbfill(I,(log10(oneSPSA_mean_1sim+oneSPSA_std_1sim))',(log10(oneSPSA_mean_1sim-oneSPSA_std_1sim))',[1 0 0],[1 0 0],0,0.2);
% %  hold on;
% %  jbfill(I,(log10(Hadamard_mean_1sim+Hadamard_std_1sim))',(log10(Hadamard_mean_1sim-Hadamard_std_1sim))',[0 1 1],[0 0 1],0,0.2);
% %  hold on;
% %  jbfill(I,(log10(Circulant_mean_1sim+Circulant_std_1sim))',(log10(Circulant_mean_1sim-Circulant_std_1sim))',[1 1 0],[0 1 0],0,0.2);
% %  hold on
%  
%  
%  plot(log10(oneSPSA_mean_1sim),'r-');
%  hold on;
%  plot(log10(Hadamard_mean_1sim),'b-');
%  hold on;
%  plot(log10(Circulant_mean_1sim),'g-');
%  hold off;
% 
% xlabel('Number of Iterations','fontweight','bold','fontsize',20);
% ylabel('log10(Normalised Mean Square Error)','fontweight','bold','fontsize',20);
% title('Comparision of the 1RDSA-1C with 1SPSA-1R and 1SPSA-1H','fontweight','bold','fontsize',20);
% legend('1SPSA-1R','1SPSA-1H','1RDSA-1C');
% 

