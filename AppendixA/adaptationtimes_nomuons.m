clear; close all; clc
path(pathdef)
addpath('../common')
addpath('../common/export_fig')

fs=constants();

density=fs.density;
decay_factor_C=fs.L14;
decay_factor_Be=fs.L10;

fs.P10_top_nm=0;
fs.P10_top_fm=0;

fs.P14_top_nm=0;
fs.P14_top_fm=0;


cm3_to_m=1/density;

times=logspace(1,5,1e3);

e_old=1.0e-5*density;
e_new_faster=e_old*10;

times_=times;

att_l_spal=fs.att_l_spal;

hh=figure(3)
set(hh,'units','centimeters','position',[0,0,14.0,14.00])

subplot(2,2,1)
hold on

%patch([times_(1) times_(end) times_(end) times_(1)],[0.01 0.01 1 1],[0.7 0.7 0.7])
set(gca,'xscale','log','yscale','log')
subplot(2,2,2)
hold on
%patch([times_(1) times_(end) times_(end) times_(1)],[0.01 0.01 1 1],[0.7 0.7 0.7])
set(gca,'xscale','log','yscale','log')

rows=2;
for j=1:rows
    %10Be
    clear e_app_Be e_app_Be_09 e_app_Be_095 e_app_Be_11 e_app_Be_105
    clear e_app_C e_app_C09 e_app_C095 e_app_C11 e_app_C105
    clear ratio_list_faster ratio_list_10minus ratio_list_10plus
    
    for i=1:length(times_)
        times=times_(i);
        
        C_ss=C_ss_calculator('10Be',e_new_faster*times,e_old,fs);
        
        C_spal=fs.P10_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_Be+e_new_faster/att_l_spal)))./(e_new_faster+att_l_spal*decay_factor_Be);
        C_nm=fs.P10_top_nm*fs.P10_att_l_nm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_nm)))./(e_new_faster+fs.P10_att_l_nm*decay_factor_Be);
        C_fm=fs.P10_top_fm*fs.P10_att_l_fm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_fm)))./(e_new_faster+fs.P10_att_l_fm*decay_factor_Be);
        
        C_content=C_ss.*exp(-times*decay_factor_Be)+C_spal+C_nm+C_fm;

        e_app_Be(i)=bisection_e_app('10Be',C_content,fs);
        e_app_Be_09(i)=bisection_e_app('10Be',C_content*0.9,fs);
        e_app_Be_11(i)=bisection_e_app('10Be',C_content*1.1,fs);
        e_app_Be_095(i)=bisection_e_app('10Be',C_content*0.95,fs);
        e_app_Be_105(i)=bisection_e_app('10Be',C_content*1.05,fs);
        e_app_Be_0975(i)=bisection_e_app('10Be',C_content*0.975,fs);
        e_app_Be_1025(i)=bisection_e_app('10Be',C_content*1.025,fs);
        
        %14C
        clear C_ss C_spal C_nm C_fm C_content
        
        C_ss=C_ss_calculator('14C',e_new_faster*times,e_old,fs);
        
        C_spal=fs.P14_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_C+e_new_faster/att_l_spal)))./(e_new_faster+att_l_spal*decay_factor_C);
        C_nm=fs.P14_top_nm*fs.P14_att_l_nm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_nm)))./(e_new_faster+fs.P14_att_l_nm*decay_factor_C);
        C_fm=fs.P14_top_fm*fs.P14_att_l_fm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_fm)))./(e_new_faster+fs.P14_att_l_fm*decay_factor_C);
        
        C_content=C_ss.*exp(-times*decay_factor_C)+C_spal+C_nm+C_fm;
        
        e_app_C(i)=bisection_e_app('14C',C_content,fs);
        e_app_C09(i)=bisection_e_app('14C',C_content*0.9,fs);
        e_app_C11(i)=bisection_e_app('14C',C_content*1.1,fs);
        e_app_C095(i)=bisection_e_app('14C',C_content*0.95,fs);
        e_app_C105(i)=bisection_e_app('14C',C_content*1.05,fs);
        e_app_C1025(i)=bisection_e_app('14C',C_content*1.025,fs);
        e_app_C0975(i)=bisection_e_app('14C',C_content*0.975,fs);
        
        ratio_list_faster(i)=e_app_C(i)/e_app_Be(i);
        ratio_list_10minus(i)=e_app_C11(i)/e_app_Be_09(i);
        ratio_list_10plus(i)=e_app_C09(i)/e_app_Be_11(i);
        ratio_list_5minus(i)=e_app_C105(i)/e_app_Be_095(i);
        ratio_list_5plus(i)=e_app_C095(i)/e_app_Be_105(i);
        ratio_list_25minus(i)=e_app_C1025(i)/e_app_Be_0975(i);
        ratio_list_25plus(i)=e_app_C0975(i)/e_app_Be_1025(i);
    end
    
    figure(3)
    
    subplot(2,rows,j)
    h1=loglog(times_,ratio_list_faster,'b','linewidth',1)
    hold on
    h2=loglog(times_,ratio_list_10minus,'r:','linewidth',1)
    h3=loglog(times_,ratio_list_5minus,'r--','linewidth',1)
    h4=loglog(times_,ratio_list_25minus,'r','linewidth',1)
    loglog(times_,ratio_list_10plus,'r:','linewidth',1)
    loglog(times_,ratio_list_5plus,'r--','linewidth',1)
    loglog(times_,ratio_list_25plus,'r','linewidth',1)
    yl=ylim;
    ylim([yl(1) yl(2)])
    ylim([0.2 5.01])
    
    set(gca,'XTick',[10 100 1000 10000 100000],'YTick',[0.2 0.5 1.00000001 2 5])
    set(gca,'TickDir','out');
    ax.TickLength = [0.04, 0.04]; % Make tick marks longer.
    ax.LineWidth = 200*0.012; % Make tick marks thicker.
    grid minor
    set(gca,'fontsize',10')
    xlim([times_(1) times_(end)])
    
    e_old=e_old*10;
    e_new_faster=e_old*10;
    
end

figure(3)


subplot(2,2,1)
title('10 mm/kyr --> 100 mm/kyr','fontsize',10)
ylabel('\epsilon_C/\epsilon_{Be}','fontsize',10)
xlabel('Time since erosion rate change [yr]','fontsize',10)

subplot(2,2,2)
title('100 mm/kyr --> 1000 mm/kyr','fontsize',10)
ylabel('\epsilon_C/\epsilon_{Be}','fontsize',10)
xlabel('Time since erosion rate change [yr]','fontsize',10)

%%
e_old=1.0e-4*density;
e_new_faster=e_old/10;


%for i=1:length(e_old_list) %looping all the erosion rates prior to the change
%    for j=1:length(e_new_fractions_faster) %looping all the erosion rates prior to the change
%e_old=e_old_list(i);
%e_new_faster=e_old*e_new_fractions_faster(j); %The vector defining the new erosion rates is calculated
figure(3)
subplot(2,2,3)
hold on
%patch([times_(1) times_(end) times_(end) times_(1)],[20 20 1 1],[0.7 0.7 0.7])
set(gca,'xscale','log','yscale','log')
subplot(2,2,4)
hold on
%patch([times_(1) times_(end) times_(end) times_(1)],[20.01 20.01 1 1],[0.7 0.7 0.7])
set(gca,'xscale','log','yscale','log')

for j=1:rows
    %10Be
    clear e_app_Be e_app_Be_09 e_app_Be_095 e_app_Be_11 e_app_Be_105
    clear e_app_C e_app_C09 e_app_C095 e_app_C11 e_app_C105
    clear ratio_list_faster ratio_list_10minus ratio_list_10plus
    
    for(i=1:length(times_))
        times=times_(i);       
        C_ss=C_ss_calculator('10Be',e_new_faster*times,e_old,fs);
        
        C_spal=fs.P10_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_Be+e_new_faster/att_l_spal)))./(e_new_faster+att_l_spal*decay_factor_Be);
        C_nm=fs.P10_top_nm*fs.P10_att_l_nm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_nm)))./(e_new_faster+fs.P10_att_l_nm*decay_factor_Be);
        C_fm=fs.P10_top_fm*fs.P10_att_l_fm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_fm)))./(e_new_faster+fs.P10_att_l_fm*decay_factor_Be);
        
        C_content=C_ss.*exp(-times*decay_factor_Be)+C_spal+C_nm+C_fm;
        
        
        e_app_Be(i)=bisection_e_app('10Be',C_content,fs);
        e_app_Be_09(i)=bisection_e_app('10Be',C_content*0.9,fs);
        e_app_Be_11(i)=bisection_e_app('10Be',C_content*1.1,fs);
        e_app_Be_095(i)=bisection_e_app('10Be',C_content*0.95,fs);
        e_app_Be_105(i)=bisection_e_app('10Be',C_content*1.05,fs);
        e_app_Be_0975(i)=bisection_e_app('10Be',C_content*0.975,fs);
        e_app_Be_1025(i)=bisection_e_app('10Be',C_content*1.025,fs);
        
        
        %14C 
        clear C_ss C_spal C_nm C_fm C_content
        C_ss=C_ss_calculator('14C',e_new_faster*times,e_old,fs);
        
        C_spal=fs.P14_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_C+e_new_faster/att_l_spal)))./(e_new_faster+att_l_spal*decay_factor_C);
        C_nm=fs.P14_top_nm*fs.P14_att_l_nm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_nm)))./(e_new_faster+fs.P14_att_l_nm*decay_factor_C);
        C_fm=fs.P14_top_fm*fs.P14_att_l_fm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_fm)))./(e_new_faster+fs.P14_att_l_fm*decay_factor_C);
        
        C_content=C_ss.*exp(-times*decay_factor_C)+C_spal+C_nm+C_fm;
        
        e_app_C(i)=bisection_e_app('14C',C_content,fs);
        e_app_C09(i)=bisection_e_app('14C',C_content*0.9,fs);
        e_app_C11(i)=bisection_e_app('14C',C_content*1.1,fs);
        e_app_C095(i)=bisection_e_app('14C',C_content*0.95,fs);
        e_app_C105(i)=bisection_e_app('14C',C_content*1.05,fs);
        e_app_C0975(i)=bisection_e_app('14C',C_content*0.975,fs);
        e_app_C1025(i)=bisection_e_app('14C',C_content*1.025,fs);
        
        ratio_list_faster(i)=e_app_C(i)/e_app_Be(i);
        ratio_list_10minus(i)=e_app_C11(i)/e_app_Be_09(i);
        ratio_list_10plus(i)=e_app_C09(i)/e_app_Be_11(i);
        ratio_list_5minus(i)=e_app_C105(i)/e_app_Be_095(i);
        ratio_list_5plus(i)=e_app_C095(i)/e_app_Be_105(i);
        ratio_list_25minus(i)=e_app_C1025(i)/e_app_Be_0975(i);
        ratio_list_25plus(i)=e_app_C0975(i)/e_app_Be_1025(i);
    end
    
    figure(3)
    
    subplot(2,2,j+2)
    h1=loglog(times_,ratio_list_faster,'linewidth',1)
    hold on
    h2=loglog(times_,ratio_list_10minus,'r:','linewidth',1)
    h3=loglog(times_,ratio_list_5minus,'r--','linewidth',1)
    loglog(times_,ratio_list_10plus,'r:','linewidth',1)
    loglog(times_,ratio_list_5plus,'r--','linewidth',1)
    h4=loglog(times_,ratio_list_25minus,'r-','linewidth',1)
    loglog(times_,ratio_list_25plus,'r-','linewidth',1)
    grid minor
    set(gca,'fontsize',10)
    set(gca,'XTick',[10 100 1000 10000 100000],'YTick',[0.2 0.5 1.00000000001 2 5])
    set(gca,'TickDir','out');
    ax.TickLength = [0.04, 0.04]; % Make tick marks longer.
    ax.LineWidth = 200*0.012; % Make tick marks thicker.
    xlim([times_(1) times_(end)])
    ylim([0.2 5.01])
    e_old=e_old*10;
    e_new_faster=e_old/10;
    
end


subplot(2,2,3)
title('100 mm/kyr --> 10 mm/kyr','fontsize',10)
ylabel('\epsilon_C/\epsilon_{Be}','fontsize',10)
xlabel('Time since erosion rate change [yr]','fontsize',10)

subplot(2,2,4)
title('1000 mm/kyr --> 100 mm/kyr','fontsize',10)
ylabel('\epsilon_C/\epsilon_{Be}','fontsize',10)
xlabel('Time since erosion rate change [yr]','fontsize',10)

subplot(2,2,4)
h=legend([h1,h2,h3,h4],'\epsilon_{C}/\epsilon_{Be}','10%','5%','2.5%')
set(h,'fontsize',10,'location','northwest')


print_string=['../Figures/AppA_adaption_times_nomuons']

save='y'
if save=='y'
    export_fig(print_string,'-pdf','-r1000','-transparent')
end
return