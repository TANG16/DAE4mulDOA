%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Georgios K. Papageorgiou
% Date: 08/02/2020
% Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast Direction-of-arrival
% Estimation of Multiple Targets Using Deep Learning and Sparse Arrays,"
% IEEE International Conference on Acoustics, Speech and Signal Processing
% (ICASSP), Barcelona, May 4-8 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results with the probabilities vs the SNR for T=16, 32 and 64 
% (scenario 1) in [1].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR_vec = 0:5:30; 
SVL = length(SNR_vec);
Prob_R_64 = zeros(1,SVL);
Prob_R_sam_64 = zeros(1,SVL);
Prob_R_est_64 = zeros(1,SVL);
Improv_est_64 = zeros(1,SVL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T = 64 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
load('RESULTS_K6_0dB_T64_res1_ongrid_prob.mat');
Prob_R_64(n) = Prob(1);
Prob_R_sam_64(n) = Prob(2);
Prob_R_est_64(n) = Prob(3);
Improv_est_64(n) = Improv;

n = 2;
load('RESULTS_K6_5dB_T64_res1_ongrid_prob.mat');
Prob_R_64(n) = Prob(1);
Prob_R_sam_64(n) = Prob(2);
Prob_R_est_64(n) = Prob(3);
Improv_est_64(n) = Improv;

n = 3;
load('RESULTS_K6_10dB_T64_res1_ongrid_prob.mat');
Prob_R_64(n) = Prob(1);
Prob_R_sam_64(n) = Prob(2);
Prob_R_est_64(n) = Prob(3);
Improv_est_64(n) = Improv;

n = 4;
load('RESULTS_K6_15dB_T64_res1_ongrid_prob.mat');
Prob_R_64(n) = Prob(1);
Prob_R_sam_64(n) = Prob(2);
Prob_R_est_64(n) = Prob(3);
Improv_est_64(n) = Improv;

n = 5;
load('RESULTS_K6_20dB_T64_res1_ongrid_prob.mat');
Prob_R_64(n) = Prob(1);
Prob_R_sam_64(n) = Prob(2);
Prob_R_est_64(n) = Prob(3);
Improv_est_64(n) = Improv;

n = 6;
load('RESULTS_K6_25dB_T64_res1_ongrid_prob.mat');
Prob_R_64(n) = Prob(1);
Prob_R_sam_64(n) = Prob(2);
Prob_R_est_64(n) = Prob(3);
Improv_est_64(n) = Improv;

n = 7;
load('RESULTS_K6_30dB_T64_res1_ongrid_prob.mat');
Prob_R_64(n) = Prob(1);
Prob_R_sam_64(n) = Prob(2);
Prob_R_est_64(n) = Prob(3);
Improv_est_64(n) = Improv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prob_R_32 = zeros(1,SVL);
Prob_R_sam_32 = zeros(1,SVL);
Prob_R_est_32 = zeros(1,SVL);
Improv_est_32 = zeros(1,SVL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T = 32 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
load('RESULTS_K6_0dB_T32_res1_ongrid_prob.mat');
Prob_R_32(n) = Prob(1);
Prob_R_sam_32(n) = Prob(2);
Prob_R_est_32(n) = Prob(3);
Improv_est_32(n) = Improv;

n = 2;
load('RESULTS_K6_5dB_T32_res1_ongrid_prob.mat');
Prob_R_32(n) = Prob(1);
Prob_R_sam_32(n) = Prob(2);
Prob_R_est_32(n) = Prob(3);
Improv_est_32(n) = Improv;

n = 3;
load('RESULTS_K6_10dB_T32_res1_ongrid_prob.mat');
Prob_R_32(n) = Prob(1);
Prob_R_sam_32(n) = Prob(2);
Prob_R_est_32(n) = Prob(3);
Improv_est_32(n) = Improv;

n = 4;
load('RESULTS_K6_15dB_T32_res1_ongrid_prob.mat');
Prob_R_32(n) = Prob(1);
Prob_R_sam_32(n) = Prob(2);
Prob_R_est_32(n) = Prob(3);
Improv_est_32(n) = Improv;

n = 5;
load('RESULTS_K6_20dB_T32_res1_ongrid_prob.mat');
Prob_R_32(n) = Prob(1);
Prob_R_sam_32(n) = Prob(2);
Prob_R_est_32(n) = Prob(3);
Improv_est_32(n) = Improv;

n = 6;
load('RESULTS_K6_25dB_T32_res1_ongrid_prob.mat');
Prob_R_32(n) = Prob(1);
Prob_R_sam_32(n) = Prob(2);
Prob_R_est_32(n) = Prob(3);
Improv_est_32(n) = Improv;

n = 7;
load('RESULTS_K6_30dB_T32_res1_ongrid_prob.mat');
Prob_R_32(n) = Prob(1);
Prob_R_sam_32(n) = Prob(2);
Prob_R_est_32(n) = Prob(3);
Improv_est_32(n) = Improv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prob_R_est_16 = zeros(1,SVL);
Improv_est_16 = zeros(1,SVL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T = 16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
load('RESULTS_K6_0dB_T16_res1_ongrid_prob.mat');
Prob_R_est_16(n) = Prob(3);
Improv_est_16(n) = Improv;

n = 2;
load('RESULTS_K6_5dB_T16_res1_ongrid_prob.mat');
Prob_R_est_16(n) = Prob(3);
Improv_est_16(n) = Improv;

n = 3;
load('RESULTS_K6_10dB_T16_res1_ongrid_prob.mat');
Prob_R_est_16(n) = Prob(3);
Improv_est_16(n) = Improv;

n = 4;
load('RESULTS_K6_15dB_T16_res1_ongrid_prob.mat');
Prob_R_est_16(n) = Prob(3);
Improv_est_16(n) = Improv;

n = 5;
load('RESULTS_K6_20dB_T16_res1_ongrid_prob.mat');
Prob_R_est_16(n) = Prob(3);
Improv_est_16(n) = Improv;

n = 6;
load('RESULTS_K6_25dB_T16_res1_ongrid_prob.mat');
Prob_R_est_16(n) = Prob(3);
Improv_est_16(n) = Improv;

n = 7;
load('RESULTS_K6_30dB_T16_res1_ongrid_prob.mat');
Prob_R_est_16(n) = Prob(3);
Improv_est_16(n) = Improv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Prob results

figure(1);
plot(SNR_vec, Prob_R_64, '-o','LineWidth', 1, 'MarkerSize', 5);
hold on;
plot(SNR_vec, Prob_R_sam_64, '-x','LineWidth', 1, 'MarkerSize', 6);
plot(SNR_vec, Prob_R_est_64, '-v','LineWidth', 1, 'MarkerSize', 5);
plot(SNR_vec, Prob_R_sam_32, '--x','LineWidth', 1, 'MarkerSize', 6);
plot(SNR_vec, Prob_R_est_32, '--v','LineWidth', 1, 'MarkerSize', 5);
plot(SNR_vec, Prob_R_est_16, '-.v','LineWidth', 1, 'MarkerSize', 5);
hold off;
xlabel('SNR [dB]', 'interpreter', 'latex');
ylabel('Probability [\%]', 'interpreter', 'latex');
legend('$\mathbf{R}_x$','$\widetilde{\mathbf{R}}_x,\ T=64$', '$\hat{\mathbf{R}}_x,\ T=64$',...
    '$\widetilde{\mathbf{R}}_x,\ T=32$', '$\hat{\mathbf{R}}_x,\ T=32$', '$\hat{\mathbf{R}}_x,\ T=16$',...
    'interpreter', 'latex');
%title('Resolution $1^\circ$','interpreter', 'latex');
%set(gca, 'YScale', 'log')
grid on;

figure(2);
plot(SNR_vec, Improv_est_64,'-v','LineWidth', 1, 'MarkerSize', 5);
hold on;
plot(SNR_vec, Improv_est_32,'-v','LineWidth', 1, 'MarkerSize', 5);
hold off;
xlabel('SNR [dB]', 'interpreter', 'latex');
ylabel('percentage [\%]', 'interpreter', 'latex');
legend('$T=64$','$T=32$', 'interpreter', 'latex');
grid on;

