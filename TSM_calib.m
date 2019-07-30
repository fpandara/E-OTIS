function [RMSE]=TSM_calib(beta,outloc,dx,dt,T,M,V,tin,RD,OD)
%% Stream characteristics
l=RD.Value(1);%Reachlength
q= RD.Value(3); %discharge (m3/s)

%% Calibration parameters

A=beta(1);
As=beta(2); 
D=beta(3); 
alpha=beta(4);

u=q/A; %average velocity
x = [0:dx:l];%x-values where you want the solution
t = [0:dt:T];%times for which you calculate solution (longer than the initial times)



%initialize concentration matrices
consmat=sparse(size(x,2),size(t,2))+RD.Value(7);
%%storage zone matrix initializations
cscons=consmat;

%%Injection concentrations%%
consin=M*1000/V;
consmat(1,1:tin)=(consin*(V/tin)+RD.Value(7)*q*1000)/(q*1000+V/tin);%conservative tracer input

% Calculate concentration at next time step
 for j=1:length(t)-1
        %Conservative tracer
        consdt(:,1)=0;
        consdt1(:,1)=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Advection-dispers1on Calculations
           consmat(end,j+1)=consmat(end-1,j);
           cscons(2:end-1,j+1)=(-1*alpha*A/As)*(cscons(2:end-1,j)-consmat(2:end-1,j))*dt+cscons(2:end-1,j)+consdt1(:,1)*dt/(24*3600);
           consmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*consmat(1:end-2,j)+(1-2*D*dt/(dx^2))*consmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*consmat(3:end,j)+consdt(:,1)*dt/(24*3600)+dt*alpha*(cscons(2:end-1,j)-consmat(2:end-1,j));
%  j
end


clf
plot(t./60,consmat(find(x==outloc),:),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
hold on;
plot(OD{:,1}./60,OD{:,2},'ks');
% set(h,'MarkerEdgeColor','none','MarkerFaceColor','k','markersize',8);
    xlabel('Time since injection (min)', 'Fontsize',15);ylabel('Tracer Concentration','Fontsize',15);
    legend('Modelled','Observed');
    ylim([min(OD{:,2}) max(OD{:,2})])
    set(gca, 'fontsize', 15);
%     title('Algal Concentration in storage zone');

shg
pause(0.1)
yhat=interp1(t./60,consmat(find(x==outloc),:),OD{:,1}./60);
RMSE = (sqrt(mean((((OD{:,2}) - (yhat))).^2))); 


