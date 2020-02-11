%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Georgios K. Papageorgiou
% Date: 08/02/2020
% Cite: [1]. G. K. Papageorgiou and M. Sellathurai, "Fast Direction-of-arrival
% Estimation of Multiple Targets Using Deep Learning and Sparse Arrays,"
% IEEE International Conference on Acoustics, Speech and Signal Processing
% (ICASSP), Barcelona, May 4-8 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is the operator uvt() as described in [1].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rx = conv2matcom(c)
%%%%%%%%%%%%%%%%%%%%%%Input%%%%%%%%%%% 
% c: real-valued vector
%%%%%%%%%%%%%%%%%%%%%Output%%%%%%%%%%%
% R_x: a Hermitian matrix according to the reverse operator of vt() 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nd  = sqrt(size(c,1));
    Ndr = Nd*(Nd+1)/2;
    re_R = zeros(Nd);
    re_R(triu(ones(Nd),0)==1) = c(1:Ndr);
    re_R = re_R' + re_R - diag(diag(re_R));
    im_R = zeros(Nd); %this is the upper triangular
    im_R(triu(ones(Nd),1)==1) = c(Ndr+1:end);
    im_R = im_R - im_R' ;
    Rx = re_R + 1j*im_R;
end