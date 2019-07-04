function [erosion_rates,depths,exposed,interval_time] = erosion_rates_depths_exposure(ti,d18O_triang,d18Oth,tDegla,k1_var,k2_mean,d18O_mean,d18O_min,d18O_max,d18O_std)

ti=ti;
%erosion_rates=k1*(d18O_triang-d18Oth)+k0;
%erosion_rates=k1_var*(d18O_triang-d18O_mean)/d18O_std+k2_mean;
%pause(200)
erosion_rates=k2_mean*(k1_var*(d18O_triang-d18O_mean)+1);
dt=diff(ti(1:2));
interval_time=dt*ones(1,size(ti,1));
interval_time(1)=0.5*dt;
interval_time(end)=0.5*dt;

erosion_rates=erosion_rates(:);

I=find(erosion_rates<0);
if (length(I)>0)
    disp('negative erosion rates in erosion_rates_depths_exposure.m')
    k1_var
    k2_mean
    return
end

interval_time=interval_time(:);

erosion_per_time=erosion_rates.*interval_time;
depths=cumsum(erosion_per_time);

%I=find(d18O_triang>d18Oth & (ti')>=tDegla);
%J=find(d18O_triang<=d18Oth | (ti')<tDegla);

%exposed=zeros(size(depths));
%exposed(I)=0;
%exposed(J)=1;

exposed=ones(size(depths));


%plot(ti*1e6,exposed)
%hold on
%plot(ti*1e6,d18O_triang)
%plot([0 2.6e6],[d18Oth d18Oth])
