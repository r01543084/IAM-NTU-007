%flux calculator
function [LF] = LF_flux(a,uR,uL)
    %分辨微觀速度(mirco_v)為正or負
   
    flux = @(w) [a a(:,end)].*w;
    apha  = max(a(:,1));%微分flux
        
    fR = flux(uR);fL = flux(uL);
    
	F = 0.5*((fL+fR)-apha*(uL-uR));%LF flux
    LF = ( F(:,2:end)-F(:,1:end-1) );
end