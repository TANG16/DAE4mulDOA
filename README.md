# DAE-for-multiple-DoA-estimation
MATLAB repository for the fast DoA estimation of multiple sources/targets using a Denoising Autoencoder (DAE).

System requirements: NVIDIA GPU (at least 2GB DDR3 memory), MATLAB 2019b.

Here is a short description of how to generate the training and testing data, train the neural network and evaluate its performance:

1. DATA_gener_DAE_allSNR.m : generates the training data and saves them to a specified directory.
2. DenoisingAutoencoder_allSNR.m : trains the DAE for each T after loading the data from the directory of step 1. The weights, i.e., net, are saved so that they can be later used for prediction. The trained weights are also provided in the repository. 
3. Testing_DATA_gener_ongrid_angles_per_SNR.m : generates the testing data for random angles on a discrete grid (integers) in the intervals specified in [1] (scenario 1).
4. Testing_DATA_gener_per_SNR.m : generates the testing data for random angles in the intervals specified in [1] (scenario 2 - continuous domain).
5. DAE_testing_per_SNR_ongrid_probabilities.m : evaluates the scheme over the testing set of step 3. Metric: probability of detection with SS-MUSIC (resolution equal to 1 degree).
6. DAE_testing_per_SNR.m : evaluates the scheme over the testing set of step 4. Metric: root-mean-squared-error (RMSE) with SS-MUSIC (resolution equal to 0.1 degree).
7. plot_Probabilities.m: plot the probability vs the SNR results after running all experiments of step 5 (all SNRs and T values).
8. plot_RMSE_res01.m: plot the RMSE vs the SNR results after running all experiments at step 6 (all SNRs and T values).

Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast direction-of-arrival estimation of multiple targets using deep learning and sparse arrays," IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Barcelona, May 4-8 2020.
