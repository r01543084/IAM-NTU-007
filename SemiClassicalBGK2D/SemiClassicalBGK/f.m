function y = f(p34)	% Basic shock tube relation equation (10.51)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riemann Solver for solving shoc-tube problems
%
% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and
% distributed under the GNU General Public License: 
% http://www.gnu.org/copyleft/gpl.html
%
% Theory is given in Section 10.2 of
% PRINCIPLES OF COMPUTATIONAL FLUID DYNAMICS, by P. Wesseling
% Springer-Verlag, Berlin etc., 2001. ISBN 3-540-67853-0
% See http://dutita0.twi.tudelft.nl/nw/users/wesseling/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global PRL  CRL MACHLEFT  gamma

wortel = sqrt(2*gamma*(gamma-1+(gamma+1)*p34));
yy = (gamma-1)*CRL*(p34-1)/wortel;
yy = (1 + MACHLEFT*(gamma-1)/2-yy)^(2*gamma/(gamma-1));
y = yy/p34 - PRL;
