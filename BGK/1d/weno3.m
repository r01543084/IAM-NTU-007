%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               WENO3   5th order                                         %
%                      Atmosphere@ntu 2013.12.2                           %
%                                   thanks Tony                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uout] = weno3(u)
%%
%beta is the same!!
um  = circshift(u,[0 1]);
umm = circshift(u,[0 2]);
up  = circshift(u,[0 -1]);
upp = circshift(u,[0 -2]);

%bc
um(:,1)=um(:,2);
umm(:,2)=umm(:,3);umm(:,1)=umm(:,2);
up(:,end)=up(:,end-1);
upp(:,end-1)=upp(:,end-2);upp(:,end)=upp(:,end-1);


beta1 = 1/3*(4*umm.^2-19*umm.*um+25*um.^2+11*umm.*u-31*um.*u+10*u.^2);    
beta2 = 1/3*(4*um.^2-13*um.*u+13*u.^2+5*um.*up-13*u.*up+4*up.^2);           
beta3 = 1/3*(10*u.^2-31*u.*up+25*up.^2+11*u.*upp-19*up.*upp+4*upp.^2);

epsilon = 10^-6;%5th order = 10^-(order+1)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   right hand side  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 |         |                              
%------------------o---------o----|----o----|----o---------o---------------
%                 i-2       i-1   |    i    |   i+1       i+2
%                                 |         |                              
%                                         -                        
%                                        u
%                                         i+1/2                            
gammaR = [1/16 5/8 5/16];

wRt1 = (gammaR(1)./((epsilon+beta1).^2));%w are 'not' real w  
wRt2 = (gammaR(2)./((epsilon+beta2).^2));
wRt3 = (gammaR(3)./((epsilon+beta3).^2));

wR1 = wRt1./(wRt1+wRt2+wRt3);%w are real w
wR2 = wRt2./(wRt1+wRt2+wRt3);
wR3 = wRt3./(wRt1+wRt2+wRt3);

uoutR = (3/8*umm-5/4*um+15/8*u).*wR1+...
        (-1/8*um+3/4*u+3/8*up).*wR2+ ...
        (3/8*u+3/4*up-1/8*upp).*wR3;
uoutR = [uoutR(:,1) uoutR];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   right hand side  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 |         |                              
%------------------o---------o----|----o----|----o---------o---------------
%                 i-2       i-1   |    i    |   i+1       i+2
%                                 |         |                              
%                                  +                        
%                                 u
%                                  i-1/2                            
gammaL = [5/16 5/8 1/16];

wLt1 = (gammaL(1)./((epsilon+beta1).^2));%w are 'not' real w 
wLt2 = (gammaL(2)./((epsilon+beta2).^2));
wLt3 = (gammaL(3)./((epsilon+beta3).^2));

wL1 = wLt1./(wLt1+wLt2+wLt3);%w are real w
wL2 = wLt2./(wLt1+wLt2+wLt3);
wL3 = wLt3./(wLt1+wLt2+wLt3);

     
uoutL = (-1/8*umm+3/4*um+3/8*u).*wL1+ ...
       (3/8*um+3/4*u-1/8*up).*wL2+   ...
       (15/8*u-5/4*up+3/8*upp).*wL3;
uoutL = [uoutL uoutL(:,end)];

uout = [uoutR;
        uoutL];
