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

%Constants are defined
[fs]=constants(); %Defines halflives, and calculates production rate parameters. Spallation is calculated using parameters from Phillips2016. Muon poduction rate parameters are calculated the formulation of Heisinger 2002, with constants from Balco 2017, model 1A alpha=1 for 10Be and 26Al, and Heisinger 2002 for 14C

times=[1e3;12e3;20e3;100e3]; %The times for which subplots are calculated

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
save='y';    
n=3e2;
min_ratio1=0.1;
        max_ratio1=10;

for jj=1:length(times) %Looping the times
    time=times(jj);
    
    
    %%
    %The SC_tester calculates the ratio e_app_14C/e_app_10Be for a
    %step function, where an erosion rate has persisted for infinite times, and
    %is then changed at the time of interest, defined by the varible time set
    %above.
    addpath('./SC')
    
    [ratio_list_faster,e_new_fractions_faster,e_old_list,ratio_list_faster_10pct,ratio_list_faster_5pct,ratio_list_faster_25pct,ratio_list_slower_10pct,ratio_list_slower_5pct,ratio_list_slower_25pct]=SC_test(fs,time,ero_min,ero_max,n);
    ratio_list_faster(isnan(ratio_list_faster))=1;
    
        %plot the 14C/10Be ratio of apparent erosion rates
    hhh=figure(2*jj)
    set(hhh,'units','centimeters','position',[0,0,19.0,24.00]);
    print_string=['../Figures/Figure6_SC_time_',num2str(time)];
    h=contourf(1E6*e_old_list/fs.density,e_new_fractions_faster,log10(ratio_list_faster),100,'LineStyle','none');
    ylim([1e-2,1e2])
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'fontsize',10)
    axis off
    shading interp
    caxis(log10([min_ratio1 max_ratio1]));
    if save=='y'
    export_fig(print_string,'-transparent','-jpeg','-r1000')
end
    
    hh=figure(1);
    set(hh,'units','centimeters','position',[0,0,19.0,24.00]);
    subplot(rows,columns,jj)
    hold on
    
    %Plot the 10% detection limit contour
    [C,h3]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_faster_10pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
    [C,h]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_slower_10pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
    %Plot the 5% detection limit contour
    [C,h2]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_faster_5pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
    [C,h]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_slower_5pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
    %Plot the 2.5% detection limit contour
    [C,h1]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_faster_25pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');
    [C,h]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_slower_25pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');
    
    if(jj==4)
        hh=legend([h1 h2 h3],'2.5%','5%','10%')
        set(hh,'fontsize',10,'location','southeast')
    end
    
    ylim([1e-2,1e2])
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'fontsize',10)
    shading interp
    
    
    if(jj==length(times))
    
        %Scale the colorbar
        
        
        for ii=1:jj
            subplot(rows,columns,ii)
            caxis(log10([min_ratio1 max_ratio1]));
        end
        
        hp4 = get(subplot(rows,columns,rows*columns),'Position');
        
        colormap('default')
        
        %Make logaritmic colorbar
        a=[0.1 0.3 1 3 10];
        a=log10(a);
        c=colorbar('FontSize',10,'YTick',a,'YTickLabel',round(10.^a,2),'Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(4)*(rows+0.1)]);
        c.Label.String='\epsilon_{C}/\epsilon_{Be}';
        c.FontSize=10;
    end
    
    view(2)
    hold on
    
    ylabel('\epsilon_{present}/\epsilon_{past}','fontsize',10)
    xlabel('\epsilon_{present} [mm/kyr]','fontsize',10)
    title([num2str(time/1e3) 'kyr'],'fontsize',10)
    
    set(gca,'ytick',[0.01 0.1 1 10 100],'xtick',[0.1 1.0 10 100 1000 10000]);
    set(gca,'TickDir','out');
    ax = gca;
    ax.TickLength = [0.02, 0.02]; % Make tick marks longer.
    ax.LineWidth = 100*0.012; % Make tick marks thicker.
end

print_string=['../Figures/Figure6_SC'];

if save=='y'
    export_fig(print_string,'-transparent','-jpeg','-r1000')
    export_fig(print_string,'-transparent','-pdf','-r1000')
end

rmpath('./SC')
return