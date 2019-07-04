function [c10Bes,c26Als,c21Nes,c14Cs,lump] = ...
    gCosmoLongsteps(ErateInt,ErateGla,tDegla,fixed_stuff)
%also saves concentration histories in "lump".
%gCosmoLongsteps2: Allows computation with fixed intervals
%Now, it is possible to compute the forward solution as a linear
%approximation:

erate_gla = ErateGla;
erate_int = ErateInt;

% Rock density in kg/m3
% rho = 2650;
%
% %Half lives
% H10=1.39e6;
% L10=log(2)/H10;
% H26=0.705e6;
% L26=log(2)/H26;
% H14=5730;
% L14=log(2)/H14;
%
% %Absorption depth scale in kg/m2
%  Tau_spal=1600;
%  Tau_nm = 9900;
%  Tau_fm = 35000;

rho=fixed_stuff.density;
L10=fixed_stuff.L10;
L26=fixed_stuff.L26;
L14=fixed_stuff.L14;

Tau_spal=fixed_stuff.att_l_spal;

Tau_10fm=fixed_stuff.P10_att_l_fm;
Tau_10nm=fixed_stuff.P10_att_l_nm;

Tau_14fm=fixed_stuff.P14_att_l_fm;
Tau_14nm=fixed_stuff.P14_att_l_nm;

Tau_26fm=fixed_stuff.P26_att_l_fm;
Tau_26nm=fixed_stuff.P26_att_l_nm;

P10_top_spal=fixed_stuff.production.P10_top_spal; %atoms/kg/yr
P10_top_nm=fixed_stuff.production.P10_top_nm; %atoms/kg/yr
P10_top_fm=fixed_stuff.production.P10_top_fm; %atoms/kg/yr


P26_top_spal=fixed_stuff.production.P26_top_spal; %atoms/kg/yr
P26_top_nm=fixed_stuff.production.P26_top_nm; %atoms/kg/yr
P26_top_fm=fixed_stuff.production.P26_top_fm; %atoms/kg/yr


% %21Ne production Advanced geomodelling (samople GU 111)
%P21_top_spal=fixed_stuff.production.P21_top_spal; %atoms/kg/yr
%P21_top_nm=fixed_stuff.production.P21_top_nm; %atoms/kg/yr
%P21_top_fm=fixed_stuff.production.P21_top_fm; %atoms/kg/yr
%

% %14C production Advanced geomodelling (samople GU 111)
P14_top_spal=fixed_stuff.production.P14_top_spal; %atoms/kg/yr
P14_top_nm=fixed_stuff.production.P14_top_nm; %atoms/kg/yr
P14_top_fm=fixed_stuff.production.P14_top_fm; %atoms/kg/yr



%>>>>> BHJ: And now the analytical solution
%INPUTS:
% tau: Decay time of nucleide
% ts: End times of intervals. Normally ts(1)=0=now, and ts(2:end)<0, so
% that ts is a decreasing vector, going more and more negative.
% zs: Present lamina depths below present surface z=0.
% ers: ers(i) is the erosion rate in interval between ts(i) and ts(i+1),
% Ks: Ks(i) is the production rate in interval between ts(i) and ts(i+1)
% L: Penetration depth of radiation in play. dC/dt(z) = K*exp(-z/L)
%OUTPUTS:
%Css: Concentrations of nucleide. Css(it,iz) applies to ts(it) and zs(iz)
%zss,tss: depths and times so that mesh(zss,tss,Css) works.%The call:
%[Css,zss,tss] = C_all_intervals(ts,zs,Ks,tau,ers,L)

% ts = -fliplr(time); %First element in ts is the start of the model run
% zs = fixed_stuff.z;


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
tau_26Al = 1/L26;
tau_21Ne = inf;
tau_14C = 1/L14;

L_spal = Tau_spal/rho; %Decay depth, exp(-z/L)


L_10nm = Tau_10nm/rho;
L_10fm = Tau_10fm/rho;

L_14nm = Tau_14nm/rho;
L_14fm = Tau_14fm/rho;

L_26nm = Tau_26nm/rho;
L_26fm = Tau_26fm/rho;

% is_ints(1) = 1; %Ice free before simulation period
% ers(1)= erate_int; %same as we used yo initiate the advective solution

K_10Be_spal = is_ints2*P10_top_spal;
K_10Be_nm   = is_ints2*P10_top_nm;
K_10Be_fm   = is_ints2*P10_top_fm;

K_26Al_spal = is_ints2*P26_top_spal;
K_26Al_nm   = is_ints2*P26_top_nm;
K_26Al_fm   = is_ints2*P26_top_fm;

%K_21Ne_spal = is_ints2*P21_top_spal;
%K_21Ne_nm   = is_ints2*P21_top_nm;
%K_21Ne_fm   = is_ints2*P21_top_fm;

K_14C_spal  = is_ints2*P14_top_spal;
K_14C_nm    = is_ints2*P14_top_nm;
K_14C_fm    = is_ints2*P14_top_fm;

Css_10Be_spal=0;
Css_10Be_nm=0;
Css_10Be_fm=0;
c10Bes=0;

Css_26Al_spal=0;
Css_26Al_nm=0;
Css_26Al_fm=0;
c26Als=0;

Css_21Ne_spal=0;
Css_21Ne_nm=0;
Css_21Ne_fm=0;
c21Nes=0;

Css_14C_spal=0;
Css_14C_nm=0;
Css_14C_fm=0;
c14Cs=0;

zobs = fixed_stuff.zobs;


[zss,tss,ExposureTimeSinceNow,dts,z0s,zs2]=zss_tss_ExposureTime_calculator(tsLong,zobs,ers2,is_ints2);
%zss
%zss./abs(tss)

for iNucl = 1:length(fixed_stuff.Nucleides)
    switch fixed_stuff.Nucleides{iNucl}
        case '10Be',
            
            
            %              [Css_10Be_spal,zss,tss,ExposureTimeSinceNow] ...
            %                  = C_all_intervals2(tsLong,zobs,K_10Be_spal,tau_10Be,ers2,L_spal);
            %              [Css_10Be_nm  ] = C_all_intervals2(tsLong,zobs,K_10Be_nm  ,tau_10Be,ers2,L_nm  );
            %              [Css_10Be_fm  ] = C_all_intervals2(tsLong,zobs,K_10Be_fm  ,tau_10Be,ers2,L_fm  );
            %              c10Bes = Css_10Be_spal(:,end)+Css_10Be_nm(:,end)+Css_10Be_fm(:,end);
            
            [Css_10Be_spal] = C_all_intervals3(tsLong,zobs,K_10Be_spal,tau_10Be,ers2,L_spal,dts,z0s,zs2);
            [Css_10Be_nm  ] = C_all_intervals3(tsLong,zobs,K_10Be_nm  ,tau_10Be,ers2,L_10nm,dts,z0s,zs2);
            [Css_10Be_fm  ] = C_all_intervals3(tsLong,zobs,K_10Be_fm  ,tau_10Be,ers2,L_10fm,dts,z0s,zs2);
            c10Bes = Css_10Be_spal(:,end)+Css_10Be_nm(:,end)+Css_10Be_fm(:,end);
            %
            % cBe_analyt = c10Be_prof(1);
        case '26Al'
            
            %Model 26Al:
            %             [Css_26Al_spal] = C_all_intervals2(tsLong,zobs,K_26Al_spal,tau_26Al,ers2,L_spal);
            %             [Css_26Al_nm  ] = C_all_intervals2(tsLong,zobs,K_26Al_nm  ,tau_26Al,ers2,L_nm  );
            %             [Css_26Al_fm  ] = C_all_intervals2(tsLong,zobs,K_26Al_fm  ,tau_26Al,ers2,L_fm  );
            %             c26Als = Css_26Al_spal(:,end)+Css_26Al_nm(:,end)+Css_26Al_fm(:,end);
            % cAl_analyt = c26Al_prof(1);
            
            [Css_26Al_spal] = C_all_intervals3(tsLong,zobs,K_26Al_spal,tau_26Al,ers2,L_spal,dts,z0s,zs2);
            [Css_26Al_nm  ] = C_all_intervals3(tsLong,zobs,K_26Al_nm  ,tau_26Al,ers2,L_26nm,dts,z0s,zs2);
            [Css_26Al_fm  ] = C_all_intervals3(tsLong,zobs,K_26Al_fm  ,tau_26Al,ers2,L_26fm,dts,z0s,zs2);
            c26Als = Css_26Al_spal(:,end)+Css_26Al_nm(:,end)+Css_26Al_fm(:,end);
            
            
            
        case '21Ne'
            
            %Model 21Ne:
            %             [Css_21Ne_spal] = C_all_intervals2(tsLong,zobs,K_21Ne_spal,tau_21Ne,ers2,L_spal);
            %             [Css_21Ne_nm  ] = C_all_intervals2(tsLong,zobs,K_21Ne_nm  ,tau_21Ne,ers2,L_nm  );
            %             [Css_21Ne_fm  ] = C_all_intervals2(tsLong,zobs,K_21Ne_fm  ,tau_21Ne,ers2,L_fm  );
            %             c21Nes = Css_21Ne_spal(:,end)+Css_21Ne_nm(:,end)+Css_21Ne_fm(:,end);
            %cNe_analyt = c21Ne_prof(1);
            
            [Css_21Ne_spal] = C_all_intervals3(tsLong,zobs,K_21Ne_spal,tau_21Ne,ers2,L_spal,dts,z0s,zs2);
            [Css_21Ne_nm  ] = C_all_intervals3(tsLong,zobs,K_21Ne_nm  ,tau_21Ne,ers2,L_nm,dts,z0s,zs2);
            [Css_21Ne_fm  ] = C_all_intervals3(tsLong,zobs,K_21Ne_fm  ,tau_21Ne,ers2,L_fm,dts,z0s,zs2);
            c21Nes = Css_21Ne_spal(:,end)+Css_21Ne_nm(:,end)+Css_21Ne_fm(:,end);
            
        case '14C'
            
            %Model 14C:
            %             [Css_14C_spal] = C_all_intervals2(tsLong,zobs,K_14C_spal,tau_14C,ers2,L_spal);
            %             [Css_14C_nm  ] = C_all_intervals2(tsLong,zobs,K_14C_nm  ,tau_14C,ers2,L_nm  );
            %             [Css_14C_fm  ] = C_all_intervals2(tsLong,zobs,K_14C_fm  ,tau_14C,ers2,L_fm  );
            %             c14Cs = Css_14C_spal(:,end)+Css_14C_nm(:,end)+Css_14C_fm(:,end);
            % cC_analyt = c14C_prof(1);
            
            [Css_14C_spal] = C_all_intervals3(tsLong,zobs,K_14C_spal,tau_14C,ers2,L_spal,dts,z0s,zs2);
            [Css_14C_nm  ] = C_all_intervals3(tsLong,zobs,K_14C_nm  ,tau_14C,ers2,L_14nm,dts,z0s,zs2);
            [Css_14C_fm  ] = C_all_intervals3(tsLong,zobs,K_14C_fm  ,tau_14C,ers2,L_14fm,dts,z0s,zs2);
            c14Cs = Css_14C_spal(:,end)+Css_14C_nm(:,end)+Css_14C_fm(:,end);
    end
end

if nargout==5
    lump.zss = zss;
    lump.ts = tss(1,:);
    lump.ExposureTimeSinceNow =ExposureTimeSinceNow;
    lump.c14Css = Css_14C_spal+Css_14C_nm+Css_14C_fm;
    lump.c10Bess = Css_10Be_spal+Css_10Be_nm+Css_10Be_fm;
    lump.c26Alss = Css_26Al_spal+Css_26Al_nm+Css_26Al_fm;
    lump.c21Ness = Css_21Ne_spal+Css_21Ne_nm+Css_21Ne_fm;
end
% disp(['analyt time:']); toc

return