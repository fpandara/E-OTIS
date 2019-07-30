function TSM_runmodel(dx,dt,T,M,V,tin,A,As,D,alpha,kt,RD,figc,xrun,trun,filen,plott)
%% Stream characteristics
l=RD.Value(1);%Reachlength
q= RD.Value(3); %discharge (m3/s)

%% Calibration parameters

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
        if figc==1
            consdt(:,1)=0;
            consdt1(:,1)=0;
        elseif figc==5
            consdt(:,1)=kt;
            consdt1(:,1)=kt;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Advection-dispers1on Calculations
           consmat(end,j+1)=consmat(end-1,j);
           cscons(2:end-1,j+1)=(-1*alpha*A/As)*(cscons(2:end-1,j)-consmat(2:end-1,j))*dt+cscons(2:end-1,j)+consdt1(:,1)*dt/(24*3600);
           consmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*consmat(1:end-2,j)+(1-2*D*dt/(dx^2))*consmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*consmat(3:end,j)+consdt(:,1)*dt/(24*3600)+dt*alpha*(cscons(2:end-1,j)-consmat(2:end-1,j));
%  j
 end
figure(1)
if plott==1
    plot(x,consmat(:,find(t==trun)),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
    xlabel('Distance from injection (m)', 'Fontsize',15);ylabel('Tracer Concentration','Fontsize',15);
    ylim([min(consmat(:,find(t==trun))) max(consmat(:,find(t==trun)))])
    set(gca, 'fontsize', 15);
elseif plott==2
    plot(t./60,consmat(find(x==xrun),:),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
    xlabel('Time since injection (m)', 'Fontsize',15);ylabel('Tracer Concentration','Fontsize',15);
    ylim([min(consmat(find(x==xrun),:)) max(consmat(find(x==xrun),:))])
    set(gca, 'fontsize', 15);    
end
xlswrite(filen,consmat)




