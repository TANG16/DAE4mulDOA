# DAE-for-multiple-DoA-estimation
Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast direction-of-arrival estimation of multiple targets using deep learning and sparse arrays," IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2020.

System requirements: NVIDIA GPU (at least 2GB DDR3 memory), MATLAB 2019b.

MATLAB repository for the fast DoA estimation of multiple sources/targets using a Denoising Autoencoder (DAE). Here is a short description of how to generate the training and testing data, train the neural network and evaluate its performance.

1. Generate the training data using DATA_gener_DAE_allSNR.m and save them to a specified directory.
2. Train the DAE using DenoisingAutoencoder_allSNR.m after loading the data from the directory of step 1.

