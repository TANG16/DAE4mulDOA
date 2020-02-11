# DAE-for-multiple-DoA-estimation
MATLAB repository for the fast DoA estimation of multiple sources/targets using a Denoising Autoencoder (DAE).

System requirements: NVIDIA GPU (at least 2GB DDR3 memory), MATLAB 2019b.

Here is a short description of how to generate the training and testing data, train the neural network and evaluate its performance:

1. DATA_gener_DAE_allSNR.m : generates the training data and saves them to a specified directory.

2. DenoisingAutoencoder_allSNR.m : trains the DAE for each T after loading the data from the directory of step 1. The weights, i.e., net, are saved so that they can be later used for prediction. The trained weights are also provided in the repository.

3. Testing data generator

4. evaluation of the scheme over the testing sets:





Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast direction-of-arrival estimation of multiple targets using deep learning and sparse arrays," IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2020.
