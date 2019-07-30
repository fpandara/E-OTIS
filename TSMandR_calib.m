function [RMSE]=TSMandR_calib(beta,outloc,dx,dt,T,M,V,tin,RD,OD,PD,trt,datt)
%% Calibration parameters
% 
A=beta(1);
As=beta(2); 
D=beta(3); 
alpha=beta(4);
Alg=beta(5);

l=RD.Value(1);%Reachlength
rchdep= RD.Value(2); %flow depth
q= RD.Value(3); %discharge (m3/s)
dayl=RD.Value(4);%daylength
ra=RD.Value(5);%total solar radiation
wtmp= RD.Value(6); %stream temperature

%% Reaction Parameters %%
lambda0=PD{'lambda0',(2)}; lambda1=PD{'lambda1',(2)};
lambda2=PD{'lambda2',(2)}; ai0=PD{'ai0',(2)}; ai1=PD{'ai1',(2)};
ai2=PD{'ai2',(2)};ai3=PD{'ai3',(2)};ai4=PD{'ai4',(2)};
ai5=PD{'ai5',(2)};ai6=PD{'ai6',(2)};k_n=PD{'k_n',(2)};
k_p=PD{'k_p',(2)};tfact=PD{'tfact',(2)};k_l=PD{'k_l',(2)};
mumax=PD{'mumax',(2)};thgra=PD{'thgra',(2)};rhoq=PD{'rhoq',(2)};
thrho=PD{'thrho',(2)};rs1=PD{'rs1',(2)}; thrs1=PD{'thrs1',(2)};
rs2=PD{'rs2',(2)};thrs2=PD{'thrs2',(2)};rs3=PD{'rs3',(2)};
thrs3=PD{'thrs3',(2)};rs4=PD{'rs4',(2)};thrs4=PD{'thrs4',(2)};
rs5=PD{'rs5',(2)};thrs5=PD{'thrs5',(2)};rk1=PD{'rk1',(2)};
thrk1=PD{'thrk2',(2)};rk2=PD{'rk2',(2)};thrk2=PD{'thrk2',(2)};
rk3=PD{'rk3',(2)};thrk3=PD{'thrk3',(2)};rk4=PD{'rk4',(2)};
thrk4=PD{'thrk4',(2)};bc1=PD{'bc1',(2)};thbc1=PD{'thbc1',(2)};
bc2=PD{'bc2',(2)};thbc2=PD{'thbc2',(2)};bc3=PD{'bc3',(2)};
thbc3=PD{'thbc3',(2)};bc4=PD{'bc4',(2)};thbc4=PD{'thbc4',(2)};
k_nx=PD{'k_nx',(2)};

%%temperature dependent constants- calculation%%
thetarhoq=rhoq*thrho^(wtmp-20);
thetars1=rs1*thrs1^(wtmp-20);
thetars2=rs2*thrs2^(wtmp-20);
thetars3=rs3*thrs3^(wtmp-20);
thetars4=rs4*thrs4^(wtmp-20);
thetars5=rs5*thrs5^(wtmp-20);
thetark1=rk1*thrk1^(wtmp-20);
thetark2=rk2*thrk2^(wtmp-20);
thetark3=rk3*thrk3^(wtmp-20);
thetark4=rk4*thrk4^(wtmp-20);
thetabc3=bc3*thbc3^(wtmp-20);
thetabc4=bc4*thbc4^(wtmp-20);


u=q/A; %average velocity
% dx=5;dt=1; %time and distance steps (m and sec)
x = [0:dx:l];%x-values where you want the solution
t = [0:dt:T];%times for which you calculate solution (longer than the initial times)

%% Background concentrations
%initialize concentration matrices


%%%%%REACTIVE TRACER%%%%%%%%

%initialize concentration matrices
algmat=sparse(size(x,2),size(t,2))+Alg; 
domat=sparse(size(x,2),size(t,2))+RD.Value(9);
nh3nmat=sparse(size(x,2),size(t,2))+RD.Value(10);
no2nmat=sparse(size(x,2),size(t,2))+RD.Value(11);
no3nmat=sparse(size(x,2),size(t,2))+RD.Value(12);
disspmat=sparse(size(x,2),size(t,2))+RD.Value(13);
cbodmat=sparse(size(x,2),size(t,2))+RD.Value(14);
orgnmat=sparse(size(x,2),size(t,2))+RD.Value(15);
orgpmat=sparse(size(x,2),size(t,2))+RD.Value(16);

%%storage zone matrix initializations
cscbod=cbodmat;
csorgn=orgnmat;
csorgp=orgpmat;
csalg=algmat;
csdo=domat;
csno3n=no3nmat;
csno2n=no2nmat;
csnh3n=nh3nmat;
csdissp=disspmat;

%%Injection concentrations%%
if trt==1
    no3nin=M*1000/V;
    no3nmat(1,1:tin)=(no3nin*(V/tin)+RD.Value(12)*q*1000)/(q*1000+V/tin);%nitrogen tracer input
elseif trt==2
    nh3nin=M*1000/V;
    nh3nmat(1,1:tin)=(nh3nin*(V/tin)+RD.Value(10)*q*1000)/(q*1000+V/tin);%nitrogen tracer input
else
    disspin=M*1000/V;
    disspmat(1,1:tin)=(disspin*(V/tin)+RD.Value(13)*q*1000)/(q*1000+V/tin);%phosphate tracer input

end

% Calculate concentration at next time step

%%photosynthetically active light intensity
algi=ra*tfact/dayl;
%%Saturation concentration for dissolved oxygen
soxy=exp(-139.34410+(1.575701*10^5/(wtmp+273.15))-6.642308*10^7/((wtmp+273.15)^2)+1.243800*10^10/((wtmp+273.15)^3)-8.621949*10^11/((wtmp+273.15)^4));
% Calculate concentration at next time step
 for j=1:length(t)-1
          %% Reactions/Interactions
        %%%%%%%%%%%%%%%%%%%%%%%Algae calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Light extinction coefficient
        lambda(:,1)=lambda0+(lambda1*ai0*algmat(2:end-1,j))+lambda2*(ai0*algmat(2:end-1,j)).^(2/3);
        
        %%algal growth limitation factors
        fnn(:,1)=(no3nmat(2:end-1,j)+no2nmat(2:end-1,j)+nh3nmat(2:end-1,j))./((no3nmat(2:end-1,j)+no2nmat(2:end-1,j)+nh3nmat(2:end-1,j))+k_n);
        fpp(:,1)=disspmat(2:end-1,j)./(disspmat(2:end-1,j)+k_p);
               
        %%algal growth attenuation factor
        fl_1(:,1)=(1./(lambda(:,1)*rchdep)).*log((k_l+algi)./(k_l+algi*(exp(-1*lambda(:,1)*rchdep))));
        fll(:,1)=0.92*(dayl/24)*fl_1(:,1);
        
        %%local algal growth rate
        gra(:,1)=mumax*fll(:,1).*min(fnn(:,1), fpp(:,1));
        
        %%algal biomass concentration change due to transformations/interactions
        algdt(:,1)=(gra(:,1)*thgra^(wtmp-20)).*algmat(2:end-1,j)-thetarhoq*algmat(2:end-1,j)-thetars1/rchdep*algmat(2:end-1,j);

        %%%%%%%%%%%%%%%%%%%%%%%Oxygen calculations%%%%%%%%%%%%%%%%%%%%%%%
        
        %Nitrification rate correction factor for low oxygen
        cordo(:,1)=1-exp(-0.6*domat(2:end-1,j));
        bc1mod(:,1)=bc1*cordo(:,1);
        bc2mod(:,1)=bc2*cordo(:,1);
        
        %%Dissolved oxygen calculation
        dodt(:,1)=thetark2*(soxy-domat(2:end-1,j))+(ai3*(gra(:,1)*thgra^(wtmp-20))-ai4 *thetarhoq).*algmat(2:end-1,j)-thetark1*cbodmat(2:end-1,j)...
            -thetark4/(rchdep*1000)-ai5*(bc1mod(:,1)*thbc1^(wtmp-20)).*nh3nmat(2:end-1,j)-ai6*(bc2mod(:,1)*thbc2^(wtmp-20)).*no2nmat(2:end-1,j);
        
        %%cbod calculations
        cboddt(:,1)=-1*(thetark1*cbodmat(2:end-1,j)-thetark3*cbodmat(2:end-1,j));
        
        %%%%%%%%%%%%%%%%%%%%%%%Nitrogen calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%organic N concentration
        orgndt(:,1)= ai1*thetarhoq*algmat(2:end-1,j)-thetabc3*orgnmat(2:end-1,j)-thetars4*orgnmat(2:end-1,j);
       
        %%New Equation (May 27 2016)%%
        f1(:,1)=((nh3nmat(2:end-1,j).*no3nmat(2:end-1,j))./((k_nx+nh3nmat(2:end-1,j)).*(k_nx+nh3nmat(2:end-1,j))))+...
                ((nh3nmat(2:end-1,j)*k_nx)./((nh3nmat(2:end-1,j)+no3nmat(2:end-1,j)).*(k_nx+no3nmat(2:end-1,j))));
%         end
        %%Ammonia nitrogen concentration
        nh3ndt(:,1)=thetabc3*orgnmat(2:end-1,j)-(bc1mod(:,1)*thbc1^(wtmp-20)).*nh3nmat(2:end-1,j)+thetars3/(rchdep*1000)-f1(:,1)*ai1.*algmat(2:end-1,j).*(gra(:,1)*thgra^(wtmp-20));
        no2ndt(:,1)=(bc1mod(:,1)*thbc1^(wtmp-20)).*nh3nmat(2:end-1,j)-(bc2mod(:,1)*thbc2^(wtmp-20)).*no2nmat(2:end-1,j);

        
        %%Nitrate nitrogen concentration
        no3ndt(:,1)=(bc2mod(:,1)*thbc2^(wtmp-20)).*no2nmat(2:end-1,j)-(1-f1(:,1))*ai1.*algmat(2:end-1,j).*(gra(:,1)*thgra^(wtmp-20));
        %%%%%%%%%%%%%%%%%%%%%%%Phosphorus calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%Organic Phosphorus calculations
        orgpdt(:,1)=ai2*thetarhoq*algmat(2:end-1,j)-thetabc4*orgpmat(2:end-1,j)-thetars5*orgpmat(2:end-1,j);
        
        %%Dissolved Phosphorus concentration
        disspdt(:,1)=thetabc4*orgpmat(2:end-1,j)+(thetars2/(rchdep*1000))-ai2*(gra(:,1)*thgra^(wtmp-20)).*algmat(2:end-1,j);

        % storage zone calculations%%
        %%%%%%%%%%%%%%%%%%%%%%Algae calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Light extinction coefficient
        lambda(:,1)=lambda0+(lambda1*ai0*csalg(2:end-1,j))+lambda2*(ai0*csalg(2:end-1,j)).^(2/3);
        
        %%algal growth limitation factors
        fnn(:,1)=(csno3n(2:end-1,j)+csno2n(2:end-1,j)+csnh3n(2:end-1,j))./((csno3n(2:end-1,j)+csno2n(2:end-1,j)+csnh3n(2:end-1,j))+k_n);
        fpp(:,1)=csdissp(2:end-1,j)./(csdissp(2:end-1,j)+k_p);
               
        %%algal growth attenuation factor
        fl_1(:,1)=(1./(lambda(:,1)*rchdep)).*log((k_l+algi)./(k_l+algi*(exp(-1*lambda(:,1)*rchdep))));
        fll(:,1)=0.92*(dayl/24)*fl_1(:,1);
        
        %%local algal growth rate
        gra(:,1)=mumax*fll(:,1).*min(fnn(:,1), fpp(:,1));
        
        %%algal biomass concentration change due to transformations/interactions
        algdt1(:,1)=(gra(:,1)*thgra^(wtmp-20)).*csalg(2:end-1,j)-thetarhoq*csalg(2:end-1,j);%-thetars1/rchdep*csalg(2:end-1,j);

        %%%%%%%%%%%%%%%%%%%%%%%Oxygen calculations%%%%%%%%%%%%%%%%%%%%%%%
        
        %Nitrification rate correction factor for low oxygen
        cordo(:,1)=1-exp(-0.6*csdo(2:end-1,j));
        bc1mod(:,1)=bc1*cordo(:,1);
        bc2mod(:,1)=bc2*cordo(:,1);
        
        %%Dissolved oxygen calculation
        dodt1(:,1)=thetark2*(soxy-csdo(2:end-1,j))+(ai3*(gra(:,1)*thgra^(wtmp-20))-ai4 *thetarhoq).*csalg(2:end-1,j)-thetark1*cscbod(2:end-1,j)...
            -thetark4/(rchdep*1000)-ai5*(bc1mod(:,1)*thbc1^(wtmp-20)).*csnh3n(2:end-1,j)-ai6*(bc2mod(:,1)*thbc2^(wtmp-20)).*csno2n(2:end-1,j);
        
        %%cbod calculations
        cboddt1(:,1)=-1*(thetark1*cscbod(2:end-1,j)-thetark3*cscbod(2:end-1,j));
        
        %%%%%%%%%%%%%%%%%%%%%%%Nitrogen calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%organic N concentration
        orgndt1(:,1)= ai1*thetarhoq*csalg(2:end-1,j)-thetabc3*csorgn(2:end-1,j)-thetars4*csorgn(2:end-1,j);
       
        %%New Equation (May 27 2016)%%
            f1(:,1)=((csnh3n(2:end-1,j).*csno3n(2:end-1,j))./((k_nx+csnh3n(2:end-1,j)).*(k_nx+csnh3n(2:end-1,j))))+...
                ((csnh3n(2:end-1,j)*k_nx)./((csnh3n(2:end-1,j)+csno3n(2:end-1,j)).*(k_nx+csno3n(2:end-1,j))));
%         end
        %%Ammonia nitrogen concentration
        nh3ndt1(:,1)=thetabc3*csorgn(2:end-1,j)-(bc1mod(:,1)*thbc1^(wtmp-20)).*csnh3n(2:end-1,j)+thetars3/(rchdep*1000)-f1(:,1)*ai1.*csalg(2:end-1,j).*(gra(:,1)*thgra^(wtmp-20));
        no2ndt1(:,1)=(bc1mod(:,1)*thbc1^(wtmp-20)).*csnh3n(2:end-1,j)-(bc2mod(:,1)*thbc2^(wtmp-20)).*csno2n(2:end-1,j);

        
        %%Nitrate nitrogen concentration
        no3ndt1(:,1)=(bc2mod(:,1)*thbc2^(wtmp-20)).*csno2n(2:end-1,j)-(1-f1(:,1))*ai1.*csalg(2:end-1,j).*(gra(:,1)*thgra^(wtmp-20));
        %%%%%%%%%%%%%%%%%%%%%%%Phosphorus calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%Organic Phosphorus calculations
        orgpdt1(:,1)=ai2*thetarhoq*csalg(2:end-1,j)-thetabc4*csorgp(2:end-1,j)-thetars5*csorgp(2:end-1,j);
        
        %%Dissolved Phosphorus concentration
        disspdt1(:,1)=thetabc4*csorgp(2:end-1,j)+(thetars2/(rchdep*1000))-ai2*(gra(:,1)*thgra^(wtmp-20)).*csalg(2:end-1,j);
     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Advection-dispers1on Calculations
            algmat(end,j+1)=algmat(end-1,j);domat(end,j+1)=domat(end-1,j);cbodmat(end,j+1)=cbodmat(end-1,j);
            orgnmat(end,j+1)=orgnmat(end-1,j);no3nmat(end,j+1)=no3nmat(end-1,j);no2nmat(end,j+1)=no2nmat(end-1,j);
            nh3nmat(end,j+1)=nh3nmat(end-1,j);orgpmat(end,j+1)=orgpmat(end-1,j);disspmat(end,j+1)=disspmat(end-1,j);
            no3nmat(end,j+1)=no3nmat(end-1,j);
              
            csalg(2:end-1,j+1)=(-1*alpha*A/As)*(csalg(2:end-1,j)-algmat(2:end-1,j))*dt+csalg(2:end-1,j)+algdt1(:,1)*dt/(24*3600);
            csdo(2:end-1,j+1)=(-1*alpha*A/As)*(csdo(2:end-1,j)-domat(2:end-1,j))*dt+csdo(2:end-1,j)+dodt1(:,1)*dt/(24*3600);
            cscbod(2:end-1,j+1)=(-1*alpha*A/As)*(cscbod(2:end-1,j)-cbodmat(2:end-1,j))*dt+cscbod(2:end-1,j)+cboddt1(:,1)*dt/(24*3600);
            csorgn(2:end-1,j+1)=(-1*alpha*A/As)*(csorgn(2:end-1,j)-orgnmat(2:end-1,j))*dt+csorgn(2:end-1,j)+orgndt1(:,1)*dt/(24*3600);
            csno3n(2:end-1,j+1)=(-1*alpha*A/As)*(csno3n(2:end-1,j)-no3nmat(2:end-1,j))*dt+csno3n(2:end-1,j)+no3ndt1(:,1)*dt/(24*3600);
            csno2n(2:end-1,j+1)=(-1*alpha*A/As)*(csno2n(2:end-1,j)-no2nmat(2:end-1,j))*dt+csno2n(2:end-1,j)+no2ndt1(:,1)*dt/(24*3600);
            csnh3n(2:end-1,j+1)=(-1*alpha*A/As)*(csnh3n(2:end-1,j)-nh3nmat(2:end-1,j))*dt+csnh3n(2:end-1,j)+nh3ndt1(:,1)*dt/(24*3600);
            csorgp(2:end-1,j+1)=(-1*alpha*A/As)*(csorgp(2:end-1,j)-orgpmat(2:end-1,j))*dt+csorgp(2:end-1,j)+orgpdt1(:,1)*dt/(24*3600);
            csdissp(2:end-1,j+1)=(-1*alpha*A/As)*(csdissp(2:end-1,j)-disspmat(2:end-1,j))*dt+csdissp(2:end-1,j)+disspdt1(:,1)*dt/(24*3600);
            
            algmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*algmat(1:end-2,j)+(1-2*D*dt/(dx^2))*algmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*algmat(3:end,j)+algdt(:,1)*dt/(24*3600)+dt*alpha*(csalg(2:end-1,j)-algmat(2:end-1,j));
            domat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*domat(1:end-2,j)+(1-2*D*dt/(dx^2))*domat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*domat(3:end,j)+dodt(:,1)*dt/(24*3600)+dt*alpha*(csdo(2:end-1,j)-domat(2:end-1,j));
            cbodmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*cbodmat(1:end-2,j)+(1-2*D*dt/(dx^2))*cbodmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*cbodmat(3:end,j)+cboddt(:,1)*dt/(24*3600)+dt*alpha*(cscbod(2:end-1,j)-cbodmat(2:end-1,j));
            orgnmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*orgnmat(1:end-2,j)+(1-2*D*dt/(dx^2))*orgnmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*orgnmat(3:end,j)+orgndt(:,1)*dt/(24*3600)+dt*alpha*(csorgn(2:end-1,j)-orgnmat(2:end-1,j));
            no3nmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*no3nmat(1:end-2,j)+(1-2*D*dt/(dx^2))*no3nmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*no3nmat(3:end,j)+no3ndt(:,1)*dt/(24*3600)+dt*alpha*(csno3n(2:end-1,j)-no3nmat(2:end-1,j));
            no2nmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*no2nmat(1:end-2,j)+(1-2*D*dt/(dx^2))*no2nmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*no2nmat(3:end,j)+no2ndt(:,1)*dt/(24*3600)+dt*alpha*(csno2n(2:end-1,j)-no2nmat(2:end-1,j));
            nh3nmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*nh3nmat(1:end-2,j)+(1-2*D*dt/(dx^2))*nh3nmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*nh3nmat(3:end,j)+nh3ndt(:,1)*dt/(24*3600)+dt*alpha*(csnh3n(2:end-1,j)-nh3nmat(2:end-1,j));
            orgpmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*orgpmat(1:end-2,j)+(1-2*D*dt/(dx^2))*orgpmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*orgpmat(3:end,j)+orgpdt(:,1)*dt/(24*3600)+dt*alpha*(csorgp(2:end-1,j)-orgpmat(2:end-1,j));
            disspmat(2:end-1,j+1)=(u*dt/(2*dx)+D*dt/(dx^2))*disspmat(1:end-2,j)+(1-2*D*dt/(dx^2))*disspmat(2:end-1,j)+(-1*u*dt/(2*dx) + D*dt/(dx^2))*disspmat(3:end,j)+disspdt(:,1)*dt/(24*3600)+dt*alpha*(csdissp(2:end-1,j)-disspmat(2:end-1,j));

 end

if datt==1
    OD=OD(:,5:6);
    OD(any(ismissing(OD),2),:)=[];
    if trt==1
        clf
        plot(x(:),no3nmat(:,find(t==outloc)),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1},OD{:,2},'ks');
        ylim([min(no3nmat(:,find(t==outloc))) max(no3nmat(:,find(t==outloc)))])
        yhat=interp1(x,no3nmat(:,find(t==outloc)),OD{:,1});
    elseif trt==2
        clf
        plot(x(:),nh3nmat(:,find(t==outloc)),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1},OD{:,2},'ks');
        ylim([min(nh3nmat(:,find(t==outloc))) max(nh3nmat(:,find(t==outloc)))])
        yhat=interp1(x,nh3nmat(:,find(t==outloc)),OD{:,1});
    elseif trt==3
        clf
        plot(x(:),disspmat(:,find(t==outloc)),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1},OD{:,2},'ks');
        ylim([min(disspmat(:,find(t==outloc))) max(disspmat(:,find(t==outloc)))])
        yhat=interp1(x,disspmat(:,find(t==outloc)),OD{:,1});
    end
elseif datt==2  
    OD=OD(:,3:4);
    OD(any(ismissing(OD),2),:)=[];
    if trt==1
        clf
        plot(t./60,no3nmat(find(x==outloc),:),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1}./60,OD{:,2},'ks');
        ylim([min(no3nmat(find(x==outloc),:)) max(no3nmat(find(x==outloc),:))])
        yhat=interp1(t./60,no3nmat(find(x==outloc),:),OD{:,1}./60);
    elseif trt==2
        clf
        plot(t./60,nh3nmat(find(x==outloc),:),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1}./60,OD{:,2},'ks');
        ylim([min(nh3nmat(find(x==outloc),:)) max(nh3nmat(find(x==outloc),:))])
        yhat=interp1(t./60,nh3nmat(find(x==outloc),:),OD{:,1}./60);
    elseif trt==3
        clf
        plot(t./60,disspmat(find(x==outloc),:),'-o','color', [0.6 0.6 0.6],'MarkerSize',2,'LineWidth',5);
        hold on;
        plot(OD{:,1}./60,OD{:,2},'ks');
        ylim([min(disspmat(find(x==outloc),:)) max(disspmat(find(x==outloc),:))])
        yhat=interp1(t./60,disspmat(find(x==outloc),:),OD{:,1}./60);
    end
end
xlabel('Time since injection (m)', 'Fontsize',15);ylabel('Tracer Concentration','Fontsize',15);
set(gca, 'fontsize', 15);
shg
pause(0.1)
RMSE = (sqrt(mean((((OD{:,2}) - (yhat))).^2))); 
 


