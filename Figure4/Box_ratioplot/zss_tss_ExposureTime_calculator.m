function [zss,tss,ExposureTimeSinceNow,dts,z0s,zs2] = zss_tss_ExposureTime_calculator(ts,zs,ers,is_ints)
ts = ts(:)'; %force row
zs = zs(:); %force column

[zss,tss]=ndgrid(zs(:),ts(:)');


%Compute surface levels at times ts
dts = diff(ts); %Time intervals, typically positive

dzs = ers.*dts; %Erosion. Note that diff(ts) are positive
z0s = zs(1)-fliplr(cumsum([0,fliplr(dzs)])); %surface levels at times ts, typically negative

zs2=zs*ones(1,length(z0s)-1);
zs2=zs2-ones(length(zs),1)*z0s(2:end);


for it=1:length(ts)-1
%     it, -dts(it), zs-z0s(it)  
%     C_produced = Cinterval(zss(it:end,it),ts(it),Ks(it),tau,ers(it),L); %concentration at end of interval produced in interval
    %C_produced = Cinterval(zs-z0s(it+1),-dts(it),Ks(it),tau,ers(it),L); %concentration at end of interval produced in interval
    %C_from_before = exp(-dts(it)/tau)*Css(:,it);%concentration at end of interval present at beginning of interval
    %Css(:,it+1)= C_produced + C_from_before; %Concentration at end of this interval
    
    %Css(:,it+1)=Cinterval(zs-z0s(it+1),-dts(it),Ks(it),tau,ers(it),L)+exp(-dts(it)/tau)*Css(:,it);
    
    zss(:,it)=zs - z0s(it); %depth at beginning if this interval
%     [C_produced(1),C_from_before(1)]
end

zss(:,length(ts)) = zs; %depths at end of last interval

% dExposure = dts.*Ks;
dExposureTime = dts.*((is_ints>0)+1e-6); %1e-6 added. ExposureTime mus be monotonic

%This is has been changed to have continous exposure
%dExposureTime=dts;

ExposureTimeSinceNow = fliplr([0,cumsum(fliplr(dExposureTime))]);
