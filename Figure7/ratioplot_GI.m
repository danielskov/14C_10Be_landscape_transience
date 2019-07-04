clear; close all; clc
path(pathdef)
%This file was written by Daniel S. Skov
%Aarhus University, 2018
%Ver. 1.1

%The script maps the ratio between apparent erosion rates for different
%nuclides given ranges of input parameters for 3 different models for erosion histories.
%The outputs outlines the parameter combinations for which landscape
%transience is detectable at different uncertainties.

addpath('../common')
addpath('../common/export_fig')
addpath('./GI')

%Constants are defined
[fs]=constants(); %Defines halflives, and calculates production rate parameters. Spallation is calculated using parameters from Phillips2016. Muon poduction rate parameters are calculated the formulation of Heisinger 2002, with constants from Balco 2017, model 1A alpha=1 for 10Be and 26Al, and Heisinger 2002 for 14C

times=[12e3;15e3;16e3;18e3];

%Automatic scaling of the subfigures
columns=ceil(sqrt(length(times)));
rows=columns;
if(rows*columns<length(times))
    columns=columns+1;
end
if((rows-1)*columns>=length(times))
    rows=rows-1;
end


min_ratio1=1E10;
min_ratio2=1E10;
min_ratio3=1E10;
min_ratio4=1E10;
max_ratio1=-1E10;
max_ratio2=-1E10;
max_ratio3=-1E10;
max_ratio4=-1E10;

ero_min=1e-7;
ero_max=1e-2;
n=3e2;
save='y';
        min_ratio2=0.3;
        max_ratio2=3;

if(times(end)>20000)
    return
end

h=figure(1);
set(h,'units','centimeters','position',[0,0,19.0,24.00])


for jj=1:length(times)
    
    time=times(jj);
    
    [glacial_erosionrates,fractions_faster,ratio_list_faster_BeC,ratio_list_faster_BeC_102,ratio_list_faster_BeC_105,ratio_list_faster_BeC_11,ratio_list_slower_BeC_102,ratio_list_slower_BeC_105,ratio_list_slower_BeC_11]=GI_tester(fs,time,ero_min,ero_max,n);
    ratio_list_faster_BeC(isnan(ratio_list_faster_BeC))=1;
    
    
    
    hhh=figure(2*jj)
    set(hhh,'units','centimeters','position',[0,0,19.0,24.00]);
    print_string=['../Figures/Figure7_GI_time_',num2str(time)];
    h=contourf(1E6*glacial_erosionrates,fractions_faster,log10(ratio_list_faster_BeC),500,'LineStyle','none');
    ylim([1e-2,1e2])
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'fontsize',10)
    axis off
    shading interp
    caxis(log10([min_ratio2 max_ratio2]));
    if save=='y'
        export_fig(print_string,'-transparent','-jpeg','-r1000')
    end
    
    
    figure(1)
    subplot(rows,columns,jj)
    hold on
    [C,h]=contour(1E6*glacial_erosionrates,fractions_faster,ratio_list_faster_BeC_102,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');
    [C,h]=contour(1E6*glacial_erosionrates,fractions_faster,ratio_list_faster_BeC_105,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
    [C,h]=contour(1E6*glacial_erosionrates,fractions_faster,ratio_list_faster_BeC_11,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
    
    [C,h1]=contour(1E6*glacial_erosionrates,fractions_faster,ratio_list_slower_BeC_102,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');
    [C,h2]=contour(1E6*glacial_erosionrates,fractions_faster,ratio_list_slower_BeC_105,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
    [C,h3]=contour(1E6*glacial_erosionrates,fractions_faster,ratio_list_slower_BeC_11,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
    if(jj==4)
        hh=legend([h1 h2 h3],'2.5%','5%','10%');
        set(hh,'fontsize',10,'location','southeast')
    end
    
    view(2)
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    set(gca,'fontsize',10);
    shading interp
    
    if(j==1)
        title('Glacial Interglacial Model')
        set(gca,'FontSize',10)
    end
    ylim([1e-2 1e2])
    
    if(jj==length(times))

        
        for ii=1:jj
            subplot(rows,columns,ii)
            caxis(log10([min_ratio2 max_ratio2]));
        end
        
        hp4 = get(subplot(rows,columns,rows*columns),'Position');
        colormap('default')
        colorbar
        a=linspace(log10(min_ratio2),log10(max_ratio2),10);
        a=[min_ratio2 0.5 1 2 3];
        a=log10(a);
        c=colorbar('FontSize',10,'YTick',a,'YTickLabel',round(10.^a,2),'Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(4)*(rows+0.1)]);
        c.Label.String='\epsilon_{C}/\epsilon_{Be}';
        c.FontSize=10;
    end
    title([num2str(time/1e3) 'kyr'],'fontsize',10);
    
    ylabel('\epsilon_{warm}/\epsilon_{cold}','fontsize',10);
    xlabel('\epsilon_{warm} [mm/kyr]','fontsize',10);
    set(gca,'ytick',[0.01 0.1 1 10 100],'xtick',[0.1 1 10 100 1000 10000]);
    set(gca,'TickDir','out');
    ax = gca;
    ax.TickLength = [0.02, 0.02]; % Make tick marks longer.
    ax.LineWidth = 100*0.012; % Make tick marks thicker.
end

print_string=['../Figures/Figure7_GI'];

if save=='y'
    export_fig(print_string,'-pdf','-r1000','-transparent')
end

rmpath('./GI')