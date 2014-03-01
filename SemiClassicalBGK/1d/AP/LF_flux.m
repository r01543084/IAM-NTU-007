%flux calculator
function [LF] = LF_flux(a,u,idv,nv)

    uR = u(1:nv,:); uL = u(nv+1:end,:);%為了使用RK所以將uR,uL一併輸入，在此作展開  
    a = a(idv);%將微觀速度做展開
    
    flux = @(w) [a a(:,end)].*w;%設定flux func.
    apha  = abs([a a(:,end)]);%微分flux
    
    fR = flux(uR);fL = flux(uL);
    F = 0.5*((fL+fR)-apha.*(uL-uR));%LF flux
    
    LF = ( F(:,2:end)-F(:,1:end-1) );
end