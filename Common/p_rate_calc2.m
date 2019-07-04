function [out]=p_rate_calc2(f_star,Natoms,sigma0,p,mindepth,maxdepth)

mc.k_neg = f_star .* 0.704 .* 0.1828; % From BCO fit,   f_star*f_C*f_D
mc.Natoms=Natoms;
mc.sigma0=sigma0;

n=2;
plotFlag=0;

% Define z vector
fitx = linspace(mindepth,maxdepth,500); % Not sure how many to use...
y = zeros(size(fitx));

% Compute "actual" production rates using model 1A
for a = 1:length(fitx);
    y(a) = P_mu_total_alpha1(fitx(a),p,mc);
end;

% Redundant assignments to output arg
out.z = fitx;
out.actP = y;

% Do fitting

% Always do 1-order fit to get starting point.
fity = log(y);
pf1 = polyfit(fitx,fity,1);
P1 = exp(pf1(2));
L1 = -1./pf1(1);
L1=L1*10;

%figure
%f=polyval(pf1,fitx)
%plot(fitx,fity,'o',fitx,f,'-')
%legend('data','linear fit')


% Starting guess
options = optimset('MaxFunEvals',50000);
x0 = [P1/1.5 P1/3 L1/2 L1.*1.5]; %[neg_P fast_P neg_L fast_L]
xopt = fminsearch(@(x) sum(((x(1).*exp(-fitx./x(3)) + x(2).*exp(-fitx./x(4)))-y).^2),x0,options);
out.P = xopt([1 2]);
out.L = xopt([3 4]);
out.predP = out.P(1).*exp(-fitx./out.L(1)) + out.P(2).*exp(-fitx./out.L(2));
ts = 'Two exponentials';

% Do plotting

if plotFlag == 1;
    figure;
    plot(out.actP,out.z,'go','markerfacecolor','g');
    hold on;
    plot(out.predP,out.z,'k');
    title(ts);
    xlabel('Pmu (atoms/g/yr)');
    set(gca,'ydir','reverse','xscale','log');
    ylabel('Depth (g/cm2)');
    figure
    plot(out.actP./out.predP,out.z)
    set(gca,'ydir','reverse');
end;


return
