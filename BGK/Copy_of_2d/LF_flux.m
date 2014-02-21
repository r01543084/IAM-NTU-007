%flux calculator
function [LF] = LF_flux(a,u,nx)

    uR = u(:,:,1:nx,:); uL = u(:,:,nx+1:end,:);%���F�ϥ�RK�ҥH�NuR,uL�@�ֿ�J�A�b���@�i�}  
    
    flux = @(w) a*w;%�]�wflux func.
    apha  = abs(a);%�L��flux
    
    fR = flux(uR);fL = flux(uL);
    F = 0.5*((fL+fR)-apha*(uL-uR));%LF flux
    
    detector = size(u);
    if detector(1)-detector(2) == -1
        %x dir
        LF = ( F(:,2:end,:,:)-F(:,1:end-1,:,:) );
    elseif detector(1)-detector(2) == 1
        LF = ( F(2:end,:,:,:)-F(1:end-1,:,:,:) );
    else
        disp('���~404')
    end
end