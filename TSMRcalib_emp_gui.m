%code for simulating Kiel tracer tests
%created by Femeena P V  starting 11/10/2016
%initial code copied from model_kiel_11102016
%modified 12/06/2016
%data for Freienwill Phosphate

%optimized params after running for EC and dissP
%A =0.6722
%As=0.1548
%K=0.2817
%alpha=0.0052
%disspmat(1:3)=57.9481;
function params=TSMRcalib_emp_gui(cData,outloc,dx,dt,T,M,V,tin,MI,algC,alg1,alg2,trt,datt)
% clear all;clc;
A =cData(1,1);
As=cData(2,1);
D=cData(3,1);
alpha=cData(4,1);
beta0=[A As D alpha algC];
% rl=140;
% beta0=[disspinp];
% rs1=1;
% rhoq=0.3;
% rs2=0.05;
% lambda0=1;
% lambda1=0.03;
% ai0=50;
% k_n=0.02;
% k_p=0.025;
% k_l=0.75;
% beta0=[rs1 rhoq rs2 lambda0 lambda1 ai0 k_n k_p k_l];
RD=readtable('InputData.xlsx','Sheet','ReachData','ReadVariableNames',true);
OD=readtable('InputData.xlsx','Sheet','ObservedData','ReadVariableNames',true);
figure(1)
options=optimset('Display','final','TolFun',1e-2,'MaxFunEvals',MI,'MaxIter',MI);
LB=[cData(:,2);alg1];
UB=[cData(:,3);alg2];
params=fminsearchbnd(@TSMandR_calib_emp,beta0,LB,UB,options,outloc,dx,dt,T,M,V,tin,RD,OD,trt,datt);

% params=fminsearchbnd(@model,beta0,[61],[61],options);






