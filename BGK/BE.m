function BE = BE(z,nu)
%% Bose-Einstein function
% Implementation of the Bose-Einstein Statistical function
% in latex words:
%
% $$F_\nu(z)=\frac{1}{\Gamma(\nu)}\int_{0}^{\infty}\frac{x^{\nu-1}}{z^{-1}e^x-1} 
%   \approx \sum_{l=1}^{\infty}\frac{z^l}{l^\nu}$$
%
% As we can notice this function is of the order $\nu$ 

L   = 1:50; % up to l = 50 to ensure a fair accuaracy

    L = repmat(L',1,length(z));
    z = repmat(z,50,1);
    BE  = sum((z.^L) ./ (L.^nu));
return