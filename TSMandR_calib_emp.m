function [RMSE]=TSMandR_calib_emp(beta,outloc,dx,dt,T,M,V,tin,RD,OD,trt,datt)
%% Calibration parameters

A=beta(1);
As=beta(2); 
D=beta(3); 
alpha=beta(4);
Decr=beta(5);


l=RD.Value(1);%Reachlength
q= RD.Value(3); %discharge (m3/s)

u=q/A; %average velocity
% dx=5;dt=1; %time and distance steps (m and sec)
x = [0:dx:l];%x-values where you want the solution
t = [0:dt:T];%times for which you calculate solution (longer than the initial times)

%% Background concentrations
%initialize concentration matrices


%%%%%REACTIVE TRACER%%%%%%%%
%initialize concentration matrices
if trt==1
    reacmat=sparse(size(x,2),size(t,2))+RD.Value(12);
    csreac=reacmat;
    reacin=M*1000/V;
    reacmat(1,1:tin)=(reacin*(V/tin)+RD.Value(12)*q*1000)/(q*1000+V/tin);%reactive tracer input
elseif trt==2
    reacmat=sparse(size(x,2),size(t,2))+RD.Value(10);
    csreac=reacmat;
    reacin=M*1000/V;
    reacmat(1,1:tin)=(reacin*(V/tin)+RD.Value(10)*q*1000)/(q*1000+V/tin);%reactive tracer input
elseif trt==3
    reacmat=sparse(size(x,2),size(t,2))+RD.Value(13);
    csreac=reacmat;
    reacin=M*1000/V;
    reacmat(1,1:tin)=(reacin*(V/tin)+RD.Value(13)*q*1000)/(q*1000+V/tin);%reactive tracer input
elseif trt==4
    reacmat=sparse(size(x,2),size(t,2))+RD.Value(17);
    csreac=reacmat;
    reacin=M*1000/V;
    reacmat(1,1:tin)=(reacin*(V/tin)+RD.Value(17)*q*1000)/(q*1000+V/tin);%reactive tracer input
end

% Calculate concentration at next time step

% Calculate concentration at next time step
 for j=1:length(t)-1
          %% Reactions/Interactions
             reacdt(:,1)=Decr*-1*reacmat(2:end-1,j);
             reacdt1(:,1)=Decr*-1*reacmat(2:end-1,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Advection-dispers1on Calculations
           reacmat(end,j+1)=reacmat(end-1,j);
           csreac(2:end-1,j+1)=(-1*alpha*A/As)*(csreac(2:end-1,j)-reacmat(2:end-1,j))*dt+csreac(2:end-1,j)+reacdt1(:,1)*dt/(24*3600);
           reacmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*reacmat(1:end-2,j)+(1-2*D*dt/(dx^2))*reacmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*reacmat(3:end,j)+reacdt(:,1)*dt/(24*3600)+dt*alpha*(csreac(2:end-1,j)-reacmat(2:end-1,j));
 
 end

if datt==1
        OD=OD(:,5:6);
        OD(any(ismissing(OD),2),:)=[];
        clf
        plot(x(:),reacmat(:,find(t==outloc)),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1},OD{:,2},'ks');
        ylim([min(reacmat(:,find(t==outloc))) max(reacmat(:,find(t==outloc)))])
        yhat=interp1(x,reacmat(:,find(t==outloc)),OD{:,1});
else  
        OD=OD(:,3:4);
        OD(any(ismissing(OD),2),:)=[];
        clf
        plot(t./60,reacmat(find(x==outloc),:),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1}./60,OD{:,2},'ks');
        ylim([min(reacmat(find(x==outloc),:)) max(reacmat(find(x==outloc),:))])
        yhat=interp1(t./60,reacmat(find(x==outloc),:),OD{:,1}./60);
end
        xlabel('Time since injection (m)', 'Fontsize',15);ylabel('Tracer Concentration','Fontsize',15);
        set(gca, 'fontsize', 15);
        shg
        pause(0.1)
        
RMSE = (sqrt(mean((((OD{:,2}) - (yhat))).^2))); 
 


