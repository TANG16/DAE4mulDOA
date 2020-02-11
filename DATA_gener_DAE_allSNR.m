%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Georgios K. Papageorgiou
% Date: 08/02/2020
% Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast Direction-of-arrival
% Estimation of Multiple Targets Using Deep Learning and Sparse Arrays,"
% IEEE International Conference on Acoustics, Speech and Signal Processing
% (ICASSP), Barcelona, May 4-8 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training DATA Generator for the scenarios introduced in [1].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

% Vectorization operator
vec = @(X) X(:);

% Parameters
pars.T = 16; % number of snapshots (samples)
Nsim = 5e+5; % number of simulations (Monte Carlo runs)
pars.fc = 2e9; % this parameter can be ignored - it is used for the evaluations of MISC array's properties
pars.c = physconst('LightSpeed'); 
pars.lambda = pars.c/pars.fc; 
pars.d = pars.lambda/2; 
SNR_dB_vec = 0:5:30; % the SNR values
SOURCE.interval = 60;
SOURCE.K = 6; % number of sources
SOURCE.power = ones(1,SOURCE.K).^2; % power of the sources/targets

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MISC.N = 6; % number of the MISC array sensors
MISC.P = 2*floor(MISC.N/4)+2; % P>= 4 should be satisfied by definition
MISC.A = [1, MISC.P-3, MISC.P*ones(1,MISC.N-MISC.P), 2*ones(1,(MISC.P-4)/2),3,2*ones(1,(MISC.P-4)/2)];
MISC.S = zeros(1,MISC.N);
for i=2:MISC.N
    MISC.S(i) = MISC.S(i-1) + MISC.A(i-1);
end
MISC.DASw = zeros(MISC.N,1);
MISC.DASw(MISC.S+1) = 1;
MISC.size = (MISC.S(end)+1)*pars.d*100; % in CM
MISC.weight = conv(MISC.DASw,flip(MISC.DASw));

% figure();
% stem(MISC.weight);
% title(['Weight Vector of a ', num2str(MISC.N),' element MISC array'], 'interpreter','latex');

% The steering vector for the MISC array for d_min = lambda/2
MISC_steer_vec = @(x) exp(1j*(2*pi*pars.d/pars.lambda)*sin(deg2rad(x))*MISC.S).';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The consecutive intervals
ang_d = linspace(-SOURCE.interval,SOURCE.interval,SOURCE.K+1);
dis = min(diff(ang_d));
% The guard band that is used as exclusion zone
safe_d = 5;

for n = 1:Nsim
% Select a random value of SNR and calculate the noise variance
SNR_dB = datasample(SNR_dB_vec,1); 
noise_power = min(SOURCE.power)*10^(-SNR_dB/10);
% The array manifold matrix
A_misc = zeros(MISC.N, SOURCE.K);
for k=1:SOURCE.K
    if k==1 % the 5 degree exclusion zone is used on the right side of the first interval
        SOURCE.angles(k, n) = ang_d(k) + (dis-safe_d)*rand;
    elseif k==SOURCE.K % the 5 degree exclusion zone is used on the left side of the last interval
        SOURCE.angles(k, n) = ang_d(k) + safe_d + (dis-safe_d)*rand;
    else % the 5 degree exclusion zone is used on the both sides of the interval
       SOURCE.angles(k, n) = ang_d(k) + safe_d + (dis-2*safe_d)*rand;
    end
    A_misc(:,k) = MISC_steer_vec(SOURCE.angles(k,n)); % this is the array manifold vector at direction \theta_k
end
ang_deg = SOURCE.angles(:,n)';

% The signal plus noise 
S = (randn(SOURCE.K,pars.T)+1j*randn(SOURCE.K,pars.T))/sqrt(2); 
X = A_misc*S;
Eta = sqrt(noise_power)*(randn(MISC.N,pars.T)+1j*randn(MISC.N,pars.T))/sqrt(2);
Y = X + Eta;

% The sample covariance matrix
R_xx_sampled = Y*Y'/pars.T; 
% The true covariance matrix
R_xx_theor = A_misc*diag(SOURCE.power)*A_misc' + noise_power*eye(MISC.N);
% Checks how far these two are
R_xx_norm_dif = norm(R_xx_theor-R_xx_sampled,'fro');

% Real and Imaginary part of the sample matrix, SIZE: MISC.N^2
[re_R_sam, im_R_sam] = conv_cov2vec(R_xx_sampled);
vR.sam(:,n) = [re_R_sam; im_R_sam];

% Real and Imaginary part for the theoretical (true) covariance matrix R, SIZE: MISC.N^2
[re_R_the, im_R_the] = conv_cov2vec(R_xx_theor);
vR.the(:,n) = [re_R_the; im_R_the];

% The true angles (ground truth)
vR.angles(:,n) = ang_deg';
n
end
% Save the Training Data to a specified directory for each T
if pars.T==16                                     
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5','/sam', size(vR.sam));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5', '/sam', vR.sam);
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5','/theor',size(vR.the));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5', '/theor', vR.the);
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5','/angles',size(vR.angles));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5', '/angles', vR.angles);
h5disp('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5');
elseif pars.T==32                                     
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T32.h5','/sam', size(vR.sam));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T32.h5', '/sam', vR.sam);
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T32.h5','/theor',size(vR.the));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T32.h5', '/theor', vR.the);
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T32.h5','/angles',size(vR.angles));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T32.h5', '/angles', vR.angles);
h5disp('C:\...\DATA_R_re_im_500K_s6_allSNR_T32.h5');
elseif pars.T==64   
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T64.h5','/sam', size(vR.sam));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T64.h5', '/sam', vR.sam);
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T64.h5','/theor',size(vR.the));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T64.h5', '/theor', vR.the);
h5create('C:\...\DATA_R_re_im_500K_s6_allSNR_T64.h5','/angles',size(vR.angles));
h5write('C:\...\DATA_R_re_im_500K_s6_allSNR_T64.h5', '/angles', vR.angles);
h5disp('C:\...\DATA_R_re_im_500K_s6_allSNR_T64.h5');   
end
