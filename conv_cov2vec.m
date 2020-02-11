%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Georgios K. Papageorgiou
% Date: 08/02/2020
% Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast Direction-of-arrival
% Estimation of Multiple Targets Using Deep Learning and Sparse Arrays,"
% IEEE International Conference on Acoustics, Speech and Signal Processing
% (ICASSP), Barcelona, May 4-8 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the operators vtr() and vti() as described in [1].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [re_R, im_R] = conv_cov2vec(R)
%%%%%%%%%%%%%%%%%%%%%%Input%%%%%%%%%%% 
% R: complex-valued Hermitian matrix
%%%%%%%%%%%%%%%%%%%%%Output%%%%%%%%%%%
% re_R: a vector of the upper triangular's real part of R 
% im_R: a vector of the strict (+1 diagonal) upper triangular's imaginary part of R 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
re_R1 = real(R);
re_R2 = triu(re_R1);
re_R2_in = triu(true(size(re_R2)));
re_R = re_R2(re_R2_in);

im_R1 = imag(R);
im_R2 = triu(im_R1,1);
im_R2_in = triu(true(size(im_R2)),1);
im_R = im_R2(im_R2_in);
end

