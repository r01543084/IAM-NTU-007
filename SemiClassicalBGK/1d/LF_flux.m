%flux calculator
function [LF] = LF_flux(a,u,idv,nv)

    uR = u(1:nv,:); uL = u(nv+1:end,:);%���F�ϥ�RK�ҥH�NuR,uL�@�ֿ�J�A�b���@�i�}  
    a = a(idv);%�N�L�[�t�װ��i�}
    
    flux = @(w) [a a(:,end)].*w;%�]�wflux func.
    apha  = abs([a a(:,end)]);%�L��flux
    
    fR = flux(uR);fL = flux(uL);
    F = 0.5*((fL+fR)-apha.*(uL-uR));%LF flux
    
    LF = ( F(:,2:end)-F(:,1:end-1) );
end