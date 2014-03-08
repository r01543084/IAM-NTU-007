function []=savefunc(T,density,p,e,marco_ux,marco_uy,tstep,counter,ID)
        
        tt = [ID,'/T/T',num2str(counter),'.mat'];
        save(tt,'T','tstep');
        tt = [ID,'/density/density',num2str(counter),'.mat'];
        save(tt,'density','tstep');
        tt = [ID,'/p/p',num2str(counter),'.mat'];
        save(tt,'p','tstep');
        tt = [ID,'/e/e',num2str(counter),'.mat'];
        save(tt,'e','tstep');
        tt = [ID,'/Ux/Ux',num2str(counter),'.mat'];
        save(tt,'marco_ux','tstep');
        tt = [ID,'/Uy/Uy',num2str(counter),'.mat'];
        save(tt,'marco_uy','tstep');
end