function [fs]=constants(jj)

% Rock density in kg/m3
fs.density = 2600;
pressures=[1013.25 898.75 794.95 701.08 616.40 540.20 471.81];
p=pressures(jj); %hPa atmospheric pressure [0m=1013.25hPa, 1000m=898.75hPa, 2000m=794.95hPa, 3000m=701.08hPa, 4000m=616.40hPa, 5000m=540.20hPa, 6000m=471.81hPa]


%Half lives
H10=1.387e6; %Chemeleff2010
fs.L10=log(2)/H10;
H26=0.705e6;
fs.L26=log(2)/H26;
H14=5730;
fs.L14=log(2)/H14;
fs.L21=0;

%Attenuation lengths from Margreth 2016
fs.att_l_spal=1600; %[kg/m2], Equal to 150g/cm2
%fs.att_l_fm=43200; %[kg/m2], Equal to 4320g/cm2
%fs.att_l_nm=15000; %[kg/m2], Equal to 1500g/cm2

p_scaling_spal=1*exp((1013.25-p)/242);

fs.P10_top_spal=p_scaling_spal*4.01e3; %SLHL Phillips 2016
fs.P26_top_spal=p_scaling_spal*27.93e3; %SLHL Phillips 2016
fs.P14_top_spal=p_scaling_spal*12.24e3; %SLHL Phillips 2016



%%
%Muons
mindepth = 0; % g/cm2
maxdepth = 7800; % g/cm2

%Muons 10Be
f_star=0.00191; %Model 1A, alpha=1;
Natoms = 2.006e22; %
sigma0 = 0.280e-30; % model 1A, alpha=1;

p_muons=p_rate_calc2(f_star,Natoms,sigma0,p,mindepth,maxdepth);

fs.P10_att_l_fm=p_muons.L(1)*10;
fs.P10_att_l_nm=p_muons.L(2)*10;
fs.P10_top_fm=p_muons.P(1)*1000;
fs.P10_top_nm=p_muons.P(2)*1000;

%Muons 14C
f_star=0.137; %Model 1A, alpha=1;
Natoms = 2.006e22; %
sigma0 = 2.37e-30; % model 1A, alpha=1;

p_muons=p_rate_calc2(f_star,Natoms,sigma0,p,mindepth,maxdepth);

fs.P14_att_l_fm=p_muons.L(2)*10;
fs.P14_att_l_nm=p_muons.L(1)*10;
fs.P14_top_fm=p_muons.P(2)*1000;
fs.P14_top_nm=p_muons.P(1)*1000;

%Muons 26Al
f_star=0.0133; %Model 1A, alpha=1;
Natoms = 1.003e22; %
sigma0 = 3.89e-30; % model 1A, alpha=1;

p_muons=p_rate_calc2(f_star,Natoms,sigma0,p,mindepth,maxdepth);

fs.P26_att_l_fm=p_muons.L(1)*10;
fs.P26_att_l_nm=p_muons.L(2)*10;
fs.P26_top_fm=p_muons.P(1)*1000;
fs.P26_top_nm=p_muons.P(2)*1000;


return
% %When using 26Al insted of 14C
% fs.P14_top_spal=fs.P26_top_spal;
% fs.P14_top_nm=fs.P26_top_nm;
% fs.P14_top_fm=fs.P26_top_fm;
% fs.L14=fs.L26;

%When using 26Al insted of 10Be
% fs.P10_top_spal=fs.P26_top_spal;
% fs.P10_top_nm=fs.P26_top_nm;
% fs.P10_top_fm=fs.P26_top_fm;
% fs.L10=fs.L26;