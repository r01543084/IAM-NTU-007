%flux calculator
function [LF] = LF_flux(a,uR,uL)
    %����L�[�t��(mirco_v)����or�t
   
    flux = @(w) [a a(:,end)].*w;
    apha  = max(a(:,1));%�L��flux
        
    fR = flux(uR);fL = flux(uL);
    
	F = 0.5*((fL+fR)-apha*(uL-uR));%LF flux
    LF = ( F(:,2:end)-F(:,1:end-1) );
end