%flux calculator
function [LF] = LF_flux(type,a,u)
switch type
    case 1
        %linear
        flux = @(w) a*w;
        apha  = max(abs(a));%·L¤Àflux
    case 2
        %non-linear
        flux = @(w) w.*w;
        apha  = max(max(abs(u)));%·L¤Àflux
end
    f = flux(u);
    
	F = 0.5*((f(2,:)+f(1,:))-apha*(u(2,:)-u(1,:)));%LF flux
    LF = ( F(2:end)-F(1:end-1) );