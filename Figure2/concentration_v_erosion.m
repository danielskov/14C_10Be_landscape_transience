%%Plots concentration v. apparent erosion rates.
%Daniel S. Skov, December 2018

clear
close all
path(pathdef)

addpath('../common')
addpath('../common/export_fig')

fs=constants();
n=1e4;

%%
%Calculate 14C
C=linspace(1e5,3e8,n);
for i=1:length(C)
    e_app(i)=bisection_e_app('14C',C(i),fs);
end

e_app_C=e_app;
C_C=C;

%%
%Calculate 10Be
C=logspace(5,12,n);

clear e_app
for i=1:length(C)
    e_app(i)=bisection_e_app('10Be',C(i),fs);
    e_app_Be_plus5(i)=bisection_e_app('10Be',C(i)*1.05,fs);
    e_app_Be_plus10(i)=bisection_e_app('10Be',C(i)*1.1,fs);
    e_app_Be_minus5(i)=bisection_e_app('10Be',C(i)*0.95,fs);
    e_app_Be_minus10(i)=bisection_e_app('10Be',C(i)*0.9,fs);
end
e_app_Be=e_app;
C_Be=C;

%%
%Calculate 26Al
clear e_app
for i=1:length(C)
    e_app(i)=bisection_e_app('26Al',C(i),fs);
    e_app_Al_plus5(i)=bisection_e_app('26Al',C(i)*1.05,fs);
    e_app_Al_plus10(i)=bisection_e_app('26Al',C(i)*1.1,fs);
    e_app_Al_minus5(i)=bisection_e_app('26Al',C(i)*0.95,fs);
    e_app_Al_minus10(i)=bisection_e_app('26Al',C(i)*0.9,fs);
end
e_app_Al=e_app;
C_Al=C;

%%
hh=figure(4)
set(hh,'units','centimeters','position',[0,0,9.0,9.00])
h1=loglog(C_C./1e3,1E6*e_app_C,'b','linewidth',1)
hold on
h2=loglog(C_Be./1e3,1E6*e_app_Be,'r','linewidth',1)
h3=loglog(C_Al./1e3,1E6*e_app_Al,'k','linewidth',1)

set(gca,'fontsize',10)
ylim([1e-1 1e6])
xlim([1e2 1e8])
h=legend([h1,h2,h3],'^{14}C','^{10}Be','^{26}Al','Location','northeast');
set(h,'fontsize',10)
title('\epsilon_{app} from CN concentrations','fontsize',10)
xlabel('CN concentration [atoms/g]','fontsize',10)
ylabel('\epsilon_{app} [mm/kyr]','fontsize',10)
set(gca,'xtick',[1E1 1E2 1E3 1E4 1E5 1e6 1e7 1e8 1e9 1e10 1e11 1e12],'ytick',[0.1 1 10 100 1000 10000 100000 1e6 1e7])

print_string=['../Figures/Figure2_conc_v_erosion'];

save='y'
if save=='y'
    export_fig(print_string,'-pdf','-r1000','-transparent')
end

return