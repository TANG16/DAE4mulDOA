%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Georgios K. Papageorgiou
% Date: 08/02/2020
% Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast Direction-of-arrival
% Estimation of Multiple Targets Using Deep Learning and Sparse Arrays,"
% IEEE International Conference on Acoustics, Speech and Signal Processing
% (ICASSP), Barcelona, May 4-8 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training of the Denoising Autoencoder (DAE) introduced in [1]. The
% networks parameters are saved so that prediction can be performed at the
% next stage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
% Loads the training data for T=16 snapshots (change the directory and the name of saved weights for other values of T)
h5disp('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5');
vR.sam = h5read('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5', '/sam');
vR.the = h5read('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5', '/theor');
vR.ang = h5read('C:\...\DATA_R_re_im_500K_s6_allSNR_T16.h5', '/angles');

% Split the training DATA
[n, N] = size(vR.sam);
% Indices
[trainInd,valInd] = dividerand(N, 0.8,0.2);

Ntrain = length(trainInd);
Nval = length(valInd);

% Split Data into training-validation sets
xTrain = vR.sam(:,trainInd);
yTrain = vR.the(:,trainInd);
xVal = vR.sam(:,valInd);
yVal = vR.the(:,valInd);

% Reshape the noisy/corrupted input
X4DTrain_noisy = reshape(xTrain, [1 n 1 Ntrain]);
X4DVal_noisy = reshape(xVal, [1 n 1 Nval]);

% Reshape the desired output
YTrain = yTrain';
YVal = yVal';

% Batch size and number of iterations
miniBatchSize = 500;
validationFrequency = floor(numel(trainInd)/miniBatchSize);

% The encoder - decoder architecture
encoder = [fullyConnectedLayer(1000) reluLayer ... 
           fullyConnectedLayer(800) reluLayer ...           
           fullyConnectedLayer(300) reluLayer ... 
           fullyConnectedLayer(100) reluLayer];

decoder = [fullyConnectedLayer(300) reluLayer ...           
           fullyConnectedLayer(800) reluLayer...
           fullyConnectedLayer(1000) reluLayer];

% The proposed DAE
layers = [imageInputLayer([1 n], 'Normalization', 'none')...
            encoder ... 
            decoder ... 
            fullyConnectedLayer(n) regressionLayer()];
        
% Training Options        
options = trainingOptions('adam', ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',4, ...
    'LearnRateDropFactor',0.5,...
    'MaxEpochs',50, ...
    'ValidationFrequency',validationFrequency,...
    'L2Regularization', 1e-4,... 
    'MiniBatchSize',miniBatchSize, ...
    'Plots','training-progress',...
    'ValidationData',{X4DVal_noisy,YVal});
% Training the DAE
net = trainNetwork(X4DTrain_noisy, YTrain, layers, options);

% Save the network weights and biases in net for T=16 (change the name for
% other values of T)
save DAE_K6_allSNR_T16 net 
