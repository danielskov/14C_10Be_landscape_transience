clear; close all; clc
path(pathdef)
%This file was written by Daniel S. Skov
%Aarhus University, 2018
%Ver. 1.1

addpath('../common')
addpath('../common/export_fig')
addpath('./GC')

%Constants are defined
[fs]=constants(); %Defines halflives, and calculates production rate parameters. Spallation is calculated using parameters from Phillips2016. Muon poduction rate parameters are calculated the formulation of Heisinger 2002, with constants from Balco 2017, model 1A alpha=1 for 10Be and 26Al, and Heisinger 2002 for 14C

min_ratio1=1E10;
min_ratio2=1E10;
min_ratio3=1E10;
min_ratio4=1E10;
max_ratio1=-1E10;
max_ratio2=-1E10;
max_ratio3=-1E10;
max_ratio4=-1E10;

%%
ero_min=1e-9;
ero_max=1e-0;

min_ratio=0.4;
max_ratio=2.5;
save='y'

n=3e2;
[ratio_list_slower,ratio_list_slower_11,ratio_list_slower_105,ratio_list_slower_1025,ratio_list_faster_11,ratio_list_faster_105,ratio_list_faster_1025,e_now,ka20_over_now_slower2]=GC_tester(fs,ero_min,ero_max,n);


hhh=figure(2)
    set(hhh,'units','centimeters','position',[0,0,14.00,9.00]);
    print_string=['../Figures/Figure8_GC'];
h=contourf(1E6*e_now,ka20_over_now_slower2,log10(ratio_list_slower),500,'LineStyle','none');    
ylim([0.04 15])
xlim([0.1 10000])
set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'fontsize',10)
    axis off
    shading interp
    caxis(log10([min_ratio max_ratio]));
    if save=='y'
        export_fig(print_string,'-transparent','-jpeg','-r1000')
    end



hh=figure(1)
set(hh,'units','centimeters','position',[0,0,14.00,9.00])
hold on
[C,h1]=contour(1E6*e_now,ka20_over_now_slower2,ratio_list_slower_11,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
[C,h2]=contour(1E6*e_now,ka20_over_now_slower2,ratio_list_slower_105,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
[C,h3]=contour(1E6*e_now,ka20_over_now_slower2,ratio_list_slower_1025,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');

[C,h]=contour(1E6*e_now,ka20_over_now_slower2,ratio_list_faster_11,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
[C,h]=contour(1E6*e_now,ka20_over_now_slower2,ratio_list_faster_105,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
[C,h]=contour(1E6*e_now,ka20_over_now_slower2,ratio_list_faster_1025,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');
plot([1e-1 1e4],[1 1],'k','Linewidth',1)

h=legend([h1 h2 h3],'10%','5%','2.5%')
set(h,'location','southeast')

view(2)
set(gca,'Xscale','log');
set(gca,'Yscale','log');
set(gca,'fontsize',10);
ylim([0.04 15])
xlim([0.1 10000])
set(gca,'ytick',[0.04 0.1 0.3 1 3 10 15],'xtick',[0.1 0.3 1 3 10 100 1000 10000]);

shading interp
ylabel('\epsilon_{now}/\epsilon_{20kyr}','FontSize',10)
xlabel('\epsilon_{now} [mm/kyr]','FontSize',10)
title('Gradual Change model','FontSize',10)
set(gca,'fontsize',10)
set(gca,'TickDir','out');
ax = gca;
properties(ax)
ax.TickLength = [0.02, 0.02]; % Make tick marks longer.
ax.LineWidth = 100*0.012; % Make tick marks thicker.


colormap('default')
min_ratio=0.4;
max_ratio=2.5;

caxis([log10(min_ratio) log10(max_ratio)]);
colorbar
a=linspace(log10(min_ratio),log10(max_ratio),10);
c=colorbar('FontSize',10,'YTick',a,'YTickLabel',round(10.^a,1));
c.Label.String='\epsilon_{C}/\epsilon_{Be}';
c.FontSize=10;

print_string=['../Figures/Figure8_GC']

save='y'
if save=='y'
    export_fig(print_string,'-pdf','-r1000','-transparent')
end

rmpath('./GC')


return