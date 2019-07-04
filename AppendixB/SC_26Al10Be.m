clear; close all; clc
path(pathdef)

addpath('../common')
addpath('../common/export_fig')
addpath('../Figure6/SC')

%Constants are defined
[fs]=constants();

times=[1e3;12e3;20e3;100e3];

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

fs.P14_att_l_fm=fs.P26_att_l_fm
fs.P14_att_l_nm=fs.P26_att_l_nm;
fs.P14_top_fm=fs.P26_top_fm;
fs.P14_top_nm=fs.P26_top_nm;
fs.P14_top_spal=fs.P26_top_spal;
fs.L14=fs.L26;

ero_min=1e-7;
ero_max=1e-2;
save='y';
n=3e2;
min_ratio1=0.5;
max_ratio1=2;


for jj=1:length(times)
    time=times(jj)
    
    [ratio_list_faster,e_new_fractions_faster,e_old_list,ratio_list_faster_10pct,ratio_list_faster_5pct,ratio_list_faster_25pct,ratio_list_slower_10pct,ratio_list_slower_5pct,ratio_list_slower_25pct]=SC_test(fs,time,ero_min,ero_max,n);
    ratio_list_faster(isnan(ratio_list_faster))=1;
    
    hhh=figure(2*jj)
    set(hhh,'units','centimeters','position',[0,0,19.0,24.00]);
    print_string=['../Figures/AppC_time_',num2str(time)];
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
    set(hh,'units','centimeters','position',[0,0,19.0,24.00])
    subplot(rows,columns,jj)
    hold on
    [C,h3]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_faster_10pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
    [C,h1]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_faster_25pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');
    [C,h2]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_faster_5pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
    
    [C,h]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_slower_10pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle',':');
    [C,h]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_slower_25pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','-');
    [C,h]=contour(1E6*e_old_list/fs.density,e_new_fractions_faster,ratio_list_slower_5pct,[1.0 1.0],'LineWidth',1,'LineColor','r','linestyle','--');
    
    if(jj==4)
        hh=legend([h1 h2 h3],'2.5%','5%','10%');
        set(hh,'fontsize',10,'location','southeast')
    end
    
    ylim([1e-2,1e2])
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'fontsize',10)
    shading interp
    
    if(jj==length(times))
        min_ratio1=0.5;
        max_ratio1=2;
        
        for ii=1:jj
            subplot(rows,columns,ii)
            caxis(log10([min_ratio1 max_ratio1]));
        end
        
        hp4 = get(subplot(rows,columns,rows*columns),'Position');
        colormap('default')
        
        a=[0.5 0.8 1 1.25 2];
        a=log10(a);
        c=colorbar('FontSize',10,'YTick',a,'YTickLabel',round(10.^a,2),'Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(4)*(rows+0.1)]);
        c.Label.String='\epsilon_{Al}/\epsilon_{Be}';
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
    properties(ax)
    ax.TickLength = [0.02, 0.02]; % Make tick marks longer.
    ax.LineWidth = 100*0.012; % Make tick marks thicker.
end

print_string=['../Figures/AppC_26Al10Be']

if save=='y'
    export_fig(print_string,'-pdf','-r1000','-transparent')
end