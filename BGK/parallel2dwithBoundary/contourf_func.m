function contourf_func(x,y,tstep,T,density,p,e,marco_ux,marco_uy,lagtime,cline)

    %plot Temperature
    subplot(2,3,1); axis([x(1),x(end),y(1),y(end)]); contourf(x,y,T,cline)
    axis equal; title({['Temperature'];['time=',num2str(tstep)]});grid on
    %plot density
    subplot(2,3,2);  axis([x(1),x(end),y(1),y(end)]);contourf(x,y,density,cline)
    axis equal; title({['Density'];['time=',num2str(tstep)]});grid on
    %plot internal energy
    subplot(2,3,3);axis([x(1),x(end),y(1),y(end)]); contourf(x,y,e,cline) 
    axis equal; title({['Total Energy'];['time=',num2str(tstep)]});grid on
    %plot pressure
    subplot(2,3,4); axis([x(1),x(end),y(1),y(end)]); contourf(x,y,p,cline); 
    axis equal; title({['Pressure'];['time=',num2str(tstep)]});grid on
    %plot Ux
    subplot(2,3,5); axis([x(1),x(end),y(1),y(end)]);contourf(x,y,marco_ux,cline)
    axis equal; title({'Marcoscopic Ux';['time=',num2str(tstep)]});grid on
    %plot Uy
    subplot(2,3,6); axis([x(1),x(end),y(1),y(end)]);contourf(x,y,marco_uy,cline)
    axis equal; title({['Marcoscopic Uy'];['time=',num2str(tstep)]});grid on
    pause(lagtime)
end