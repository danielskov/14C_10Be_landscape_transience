function [c10Bes,c14Cs,lump] = ...
    gCosmoLongsteps(ErateInt,ErateGla,tDegla,fixed_stuff)
%also saves concentration histories in "lump".
%gCosmoLongsteps2: Allows computation with fixed intervals
%Now, it is possible to compute the forward solution as a linear
%approximation:

erate_gla = ErateGla;
erate_int = ErateInt;

rho=fixed_stuff.density;
L10=fixed_stuff.L10;
L14=fixed_stuff.L14;

Tau_spal=fixed_stuff.att_l_spal;

Tau_10fm=fixed_stuff.P10_att_l_fm;
Tau_10nm=fixed_stuff.P10_att_l_nm;

Tau_14fm=fixed_stuff.P14_att_l_fm;
Tau_14nm=fixed_stuff.P14_att_l_nm;

P10_top_spal=fixed_stuff.production.P10_top_spal; %atoms/kg/yr
P10_top_nm=fixed_stuff.production.P10_top_nm; %atoms/kg/yr
P10_top_fm=fixed_stuff.production.P10_top_fm; %atoms/kg/yr

% %14C production Advanced geomodelling (samople GU 111)
P14_top_spal=fixed_stuff.production.P14_top_spal; %atoms/kg/yr
P14_top_nm=fixed_stuff.production.P14_top_nm; %atoms/kg/yr
P14_top_fm=fixed_stuff.production.P14_top_fm; %atoms/kg/yr


tStarts  = fixed_stuff.tStarts;
relExpos = fixed_stuff.relExpos;
relExpos(end+1) = 1;      % !!!!!! We hereby impose full exposure during the Holocene
%d18Oth should not determine the "start of Holocene": update by tDegla
tStarts(end+1) = -tDegla;
tsLong = [-inf,tStarts(:)',0];
is_ints2 = relExpos; %If relExpos==1 then is_ints2 is 1, i.e. it is an interglacial
is_ints2 = [1,is_ints2(:)']; %Exposed before the Quaternary

ers2 = erate_int*is_ints2 + erate_gla*(1-is_ints2);
is_ints2=ones(length(is_ints2),1)';

tau_10Be = 1/L10;
tau_14C = 1/L14;

L_spal = Tau_spal/rho; %Decay depth, exp(-z/L)

L_10nm = Tau_10nm/rho;
L_10fm = Tau_10fm/rho;

L_14nm = Tau_14nm/rho;
L_14fm = Tau_14fm/rho;



K_10Be_spal = is_ints2*P10_top_spal;
K_10Be_nm   = is_ints2*P10_top_nm;
K_10Be_fm   = is_ints2*P10_top_fm;

K_14C_spal  = is_ints2*P14_top_spal;
K_14C_nm    = is_ints2*P14_top_nm;
K_14C_fm    = is_ints2*P14_top_fm;

Css_10Be_spal=0;
Css_10Be_nm=0;
Css_10Be_fm=0;
c10Bes=0;

Css_14C_spal=0;
Css_14C_nm=0;
Css_14C_fm=0;
c14Cs=0;

zobs = fixed_stuff.zobs;


[zss,tss,ExposureTimeSinceNow,dts,z0s,zs2]=zss_tss_ExposureTime_calculator(tsLong,zobs,ers2,is_ints2);

for iNucl = 1:length(fixed_stuff.Nucleides)
    switch fixed_stuff.Nucleides{iNucl}
        case '10Be',
            
            [Css_10Be_spal] = C_all_intervals3(tsLong,zobs,K_10Be_spal,tau_10Be,ers2,L_spal,dts,z0s,zs2);
            [Css_10Be_nm  ] = C_all_intervals3(tsLong,zobs,K_10Be_nm  ,tau_10Be,ers2,L_10nm,dts,z0s,zs2);
            [Css_10Be_fm  ] = C_all_intervals3(tsLong,zobs,K_10Be_fm  ,tau_10Be,ers2,L_10fm,dts,z0s,zs2);
            c10Bes = Css_10Be_spal(:,end)+Css_10Be_nm(:,end)+Css_10Be_fm(:,end);
        case '14C'
            [Css_14C_spal] = C_all_intervals3(tsLong,zobs,K_14C_spal,tau_14C,ers2,L_spal,dts,z0s,zs2);
            [Css_14C_nm  ] = C_all_intervals3(tsLong,zobs,K_14C_nm  ,tau_14C,ers2,L_14nm,dts,z0s,zs2);
            [Css_14C_fm  ] = C_all_intervals3(tsLong,zobs,K_14C_fm  ,tau_14C,ers2,L_14fm,dts,z0s,zs2);
            c14Cs = Css_14C_spal(:,end)+Css_14C_nm(:,end)+Css_14C_fm(:,end);
    end
end

lump.zss = zss;
lump.ts = tss(1,:);
lump.c14Css = Css_14C_spal+Css_14C_nm+Css_14C_fm;
lump.c10Bess = Css_10Be_spal+Css_10Be_nm+Css_10Be_fm;

return