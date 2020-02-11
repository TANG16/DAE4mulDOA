%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Georgios K. Papageorgiou
% Date: 08/02/2020
% Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast Direction-of-arrival
% Estimation of Multiple Targets Using Deep Learning and Sparse Arrays,"
% IEEE International Conference on Acoustics, Speech and Signal Processing
% (ICASSP), Barcelona, May 4-8 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the performance of the DAE for each SNR level and number of snapshots 
% T in terms of root-mean-squared-error (RMSE).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

% Load the testing DATA
h5disp('C:\...\Testing_DATA\Test_DATA_R_re_im_75K_s6_30dB_T16.h5');
z_sam = h5read('C:\...\Testing_DATA\Test_DATA_R_re_im_75K_s6_30dB_T16.h5', '/sam');
z_the = h5read('C:\...\Testing_DATA\Test_DATA_R_re_im_75K_s6_30dB_T16.h5', '/theor');
% This is the ground truth that is used for the evaluation
True_angles = h5read('C:\...\Testing_DATA\Test_DATA_R_re_im_75K_s6_30dB_T16.h5', '/angles');

% Predict the vectors r
[n, N_test] = size(z_sam);
Xtest_noisy = reshape(z_sam, [1 n 1 N_test]);
load('DAE_K6_allSNR_T16.mat'); % load the trained weights of the proposed DAE
r_est = predict(net, Xtest_noisy)';
SOURCE_K = size(True_angles,1); % the number of sources/targets
res = 0.1; % SS-MUSIC's resolution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MISC array
MISC.N = 6; % number of the MISC array's sensors
MISC.P = 2*floor(MISC.N/4)+2; % P>= 4 should be satisfied
MISC.A = [1, MISC.P-3, MISC.P*ones(1,MISC.N-MISC.P), 2*ones(1,(MISC.P-4)/2),3,2*ones(1,(MISC.P-4)/2)];
MISC.S = zeros(1,MISC.N);
for i=2:MISC.N
    MISC.S(i) = MISC.S(i-1) + MISC.A(i-1);
end
MISC.DASw = zeros(MISC.N,1);
MISC.DASw(MISC.S+1) = 1;
MISC.weights = conv(MISC.DASw,flip(MISC.DASw));

MISC.Nva = MISC.N^2/2+3*MISC.N-9;% for the case that MISC.N%2=0
L = (MISC.Nva-1)/2; % the smothing parameter (maximum number of sources that can be identified)

% figure();
% stem(MISC.weight);
% title(['Weight Vector of a ', num2str(MISC.N),' element MISC array'], 'interpreter','latex');

% Define the difference coarray of MISC corresponding to a virtual ULA
MISC.S_dif_misc = zeros(MISC.N,MISC.N);
for n=1:MISC.N
   MISC.S_dif_misc(:,n) = -MISC.S(n) + MISC.S;
end
MISC.Set_D_misc = unique(vec(MISC.S_dif_misc)).';
MISC.DOF_dif_misc = numel(MISC.Set_D_misc);
[c, ia] = unique(vec(MISC.S_dif_misc),'sorted'); % identify the unique rows of the difference co-array
% The extraction of repeated entries of the difference co-array manifold
% via the selection matrix
J = zeros(MISC.DOF_dif_misc,MISC.N^2);
for ii=1:MISC.DOF_dif_misc
   J(ii,ia(ii))=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
rmse_the = 0;
rmse_sam = 0;
rmse_est = 0;
% use parfor or simply for if parallel processing is not available
for n=1:N_test
   
   % The smoothed true covariance matrix
   Rx = conv2matcom(z_the(:,n));
   Rz = vec(Rx); 
   Rzu = J*Rz; 
   Rss_z  = spsmooth(Rzu*Rzu',L+1);
    
   % The smoothed sample covariance matrix
   Rx_sam = conv2matcom(z_sam(:,n));
   Rz_sam = vec(Rx_sam); 
   Rzu_sam = J*Rz_sam; 
   Rss_z_sam  = spsmooth(Rzu_sam*Rzu_sam',L+1);
   
   % The smoothed predicted covariance matrix
   R_est = conv2matcom(r_est(:,n));
   Rz_est = vec(R_est); 
   Rzu_est = J*Rz_est;
   Rss_z_est  = spsmooth(Rzu_est*Rzu_est',L+1);
   
    % Estimate the angles with MUSIC
   [doas_the, spec_the, specang_the] = musicdoa(Rss_z,SOURCE_K,'ScanAngles', -90:res:90);
   % Comment this line for T=16 since MUSIC fails to resolve the angles -
   % uncomment for the rest of the cases
   %    [doas_sam, spec_sam, specang_sam] = musicdoa(Rss_z_sam,SOURCE_K, 'ScanAngles', -90:res:90);
   [doas_est, spec_est, specang_est] = musicdoa(Rss_z_est,SOURCE_K, 'ScanAngles', -90:res:90);
   
   ang_the = sort(doas_the)';
   % Comment this line for T=16 since MUSIC fails to resolve the angles -
   % uncomment for the rest of the cases
   %    ang_sam = sort(doas_sam)';
   ang_est = sort(doas_est)'; 
   ang_gt = sort(True_angles(:,n));
   
   % RMSE calculation
   rmse_the = rmse_the + norm(ang_the - ang_gt)^2;
   % Comment this line for T=16 since MUSIC fails to resolve the angles -
   % uncomment for the rest of the cases
   %    rmse_sam = rmse_sam + norm(ang_sam - ang_gt)^2;
   rmse_est = rmse_est + norm(ang_est - ang_gt)^2;
   
%    figure(1);
%    plot(spec_the);
%    hold on;
%    plot(ang_the,'.');
%    hold off;
   
   n
end

% The average RMSE
RMSE_the = sqrt(rmse_the/SOURCE_K/N_test);
RMSE_sam = sqrt(rmse_sam/SOURCE_K/N_test);
RMSE_est = sqrt(rmse_est/SOURCE_K/N_test);

% A table with the results
RMSE = [ RMSE_the; RMSE_sam; RMSE_est];
VarNames = {'Theoretical'; 'Sampled'; 'Estimated'};
Tab = table(VarNames, RMSE)

Improv = (RMSE_sam - RMSE_est)*100/RMSE_sam

% Save the results per SNR and T (rename accordingly)
save('RESULTS_K6_30dB_T16_res01.mat','RMSE','Improv');