%C_all_interval2.m
%Note: C_all_intervals2 has starting concentration zero.
%Thus, if you want to start with steady state you must specify
% ts(1)=-inf (easy!)
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
%zss,tss: depths and times so that mesh(zss,tss,Css) works.



%function [Css,zss,tss,ExposureTimeSinceNow] = C_all_intervals2(ts,zs,Ks,tau,ers,L)
function [Css] = C_all_intervals3(ts,zs,Ks,tau,ers,L,dts,z0s,zs2)

if isinf(tau)
    tau=1/eps; %to avoid inf/inf below,
end

%ts = ts(:)'; %force row
%zs = zs(:); %force column
%[zss,tss]=ndgrid(zs(:),ts(:)');

%Compute surface levels at times ts
%dts = diff(ts); %Time intervals, typically positive
%dzs = ers.*dts; %Erosion. Note that diff(ts) are positive
%z0s = zs(1)-fliplr(cumsum([0,fliplr(dzs)])); %surface levels at times ts, typically negative

%zs2=zs*ones(1,length(z0s)-1);
%zs2=zs2-ones(length(zs),1)*z0s(2:end);

Css=zeros(length(zs),length(ts)); %Pre allocation of Css matrix, which contains the nuclide concentrations through time

decay_rate=exp(-dts/tau);
%Css(:,1)=Cinf2(zs-z0s(1),Ks(1),tau,ers(1),L);
t=0;
for it=1:length(ts)-1
    %C_produced = Cinterval2(zs-z0s(it+1),-dts(it),Ks(it),tau,ers(it),L); %concentration at end of interval produced in interval
    %C_from_before = exp(-dts(it)/tau)*Css(:,it);%concentration at end of interval present at beginning of interval
    %Css(:,it+1)= C_produced + C_from_before %Concentration at end of this interval
    %Css(:,it+1)=Cinterval2(zs-z0s(it+1),-dts(it),Ks(it),tau,ers(it),L)+exp(-dts(it)/tau)*Css(:,it);
    %Css(:,it+1)=Cinterval2(zs-z0s(it+1),-dts(it),Ks(it),tau,ers(it),L)+decay_rate(it)*Css(:,it);
    Css(:,it+1)= Ks(it)/(1/tau+ers(it)/L)*exp(-((zs-z0s(it+1))-ers(it)*t)/L)*(1-exp(-(1/tau+ers(it)/L)*(t+dts(it))))+decay_rate(it)*Css(:,it);
 end

%zss(:,length(ts)) = zs; %depths at end of last interval

% dExposure = dts.*Ks;
%dExposureTime = dts.*((Ks>0)+1e-6); %1e-6 added. ExposureTime mus be monotonic
%ExposureTimeSinceNow = fliplr([0,cumsum(fliplr(dExposureTime))]);

%disp('ny');

function C = Cinf2(zs,K,tau,er,L)
%Steady state concentration with constant erosion since t=-inf
a=er;
t=0; 
C = K/(1/tau+a/L)*exp(-(zs-a*t)/L);

function C = Cinterval2(zs,t0,K,tau,er,L)
%Steady state concentration with constant erosion since t=t0
%If both decay times and erosion rates are very small, then this expression
%is unstable. Make Taylor expansion to first order, and the result will probably be
%C=K*t0*exp(-zss/L);
a=er;
t = 0;
C = K/(1/tau+a/L)*exp(-(zs-a*t)/L)*(1-exp(-(1/tau+a/L)*(t-t0)));
%C = K/(1/tau+er/L)*exp(-(zs-er*0)/L)*(1-exp(-(1/tau+er/L)*(0-t0)));
%C=K/(1/tau+er/L)*exp(-zs/L)*(1-exp(t0/(tau+er/L)));

%zs is two columns, that is why we get to column concentrations

