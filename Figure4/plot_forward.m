%Plots the forward models, and corresponding apparent erosion rates at
%given times
%Daniel S. Skov, December 2018
close all; clear

path(pathdef)
addpath('../common')
addpath('../common/export_fig')

load lisiecki_triinterp_2p6Ma_5ky.mat %Load the d18O curve (5kyr triangle smoothed)
tdegla=12; %The timing of the last change in erosion rate (in the SC and GI models)

%Define constants
fs=constants();
density=fs.density; %Density
decay_factor_C=fs.L14; %Decay constant 14C
decay_factor_Be=fs.L10; %Decay constant 10Be

%Production rates for 10Be and 14C
P10_top_spal=fs.P10_top_spal;
P14_top_spal=fs.P14_top_spal;
P10_top_nm=fs.P10_top_nm;
P10_top_fm=fs.P10_top_fm;
P14_top_nm=fs.P14_top_nm;
P14_top_fm=fs.P14_top_fm;

h=figure(2)
set(h,'units','centimeters','position',[0,0,19.00,11.25])


%%
%Plot forward SC model
plot(ti*1000,d18O_triang,'linewidth',1,'color','r')
set(gca,'fontsize',10)
xlim([0 200])
set(gca,'Xtick',[0 50 100 150 200])

%%
%Calculate apparent erosion rates in SC model
old_erosion_rate=1e-5;
new_erosion_rate=old_erosion_rate*10;

subplot(5,1,1)
plot([0 tdegla tdegla 300],[1E6*new_erosion_rate 1E6*new_erosion_rate 1E6*old_erosion_rate 1E6*old_erosion_rate],'linewidth',1)
set(gca,'fontsize',10,'yscale','log')
ylabel('\epsilon [mm/kyr]','fontsize',10)
xlim([0 200])
ylim([1E6*5e-6 1E6*3e-4])
set(gca,'Xtick',[0 50 100 150 200],'Ytick',[10.00000000001 100])

times=linspace(0,tdegla*1e3,1e3);

e_old=old_erosion_rate*fs.density;
e_new_faster=new_erosion_rate*fs.density;

times_=times;


att_l_spal=fs.att_l_spal;

for i=1:length(times_) %Calculate apparent erosion rates for times after tchange
    times=times_(i);
    
    %10Be
    C_ss=C_ss_calculator('10Be',e_new_faster*times,e_old,fs); %Calculate steady state erosion rate at the depth of the sample prior to change in erosion rate
    
    C_spal=fs.P10_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_Be+e_new_faster/att_l_spal)))./(e_new_faster+att_l_spal*decay_factor_Be);
    C_nm=fs.P10_top_nm*fs.P10_att_l_nm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_nm)))./(e_new_faster+fs.P10_att_l_nm*decay_factor_Be);
    C_fm=fs.P10_top_fm*fs.P10_att_l_fm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_fm)))./(e_new_faster+fs.P10_att_l_fm*decay_factor_Be);
    
    C_content=C_ss.*exp(-times*decay_factor_Be)+C_spal+C_nm+C_fm;
    
    e_app_Be(i)=bisection_e_app('10Be',C_content,fs);
    
    %14C
    clear C_ss C_spal C_nm C_fm C_content
    C_ss=C_ss_calculator('14C',e_new_faster*times,e_old,fs);
    
    C_spal=fs.P14_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_C+e_new_faster/att_l_spal)))./(e_new_faster+att_l_spal*decay_factor_C);
    C_nm=fs.P14_top_nm*fs.P14_att_l_nm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_nm)))./(e_new_faster+fs.P14_att_l_nm*decay_factor_C);
    C_fm=fs.P14_top_fm*fs.P14_att_l_fm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_fm)))./(e_new_faster+fs.P14_att_l_fm*decay_factor_C);
    
    C_content=C_ss.*exp(-times*decay_factor_C)+C_spal+C_nm+C_fm;
    e_app_C(i)=bisection_e_app('14C',C_content,fs);
    
    ratio_list_faster(i)=e_app_C(i)/e_app_Be(i);
end

%Plot the evolution of apparent erosion rates for 10Be and 14C
hold on
plot([200 tdegla-times_./1000],[1E6*old_erosion_rate 1E6*e_app_Be],'r','linewidth',1)
plot([200 tdegla-times_./1000],[1E6*old_erosion_rate 1E6*e_app_C],'g','linewidth',1)
ylim([1E6*old_erosion_rate/3 1E6*new_erosion_rate*3])

%%
%Plot the ratio of apparent erosion rates for the SC model
subplot(5,1,4)

semilogy([200 tdegla-times_./1000],[1 ratio_list_faster],'m','linewidth',1)
hold on
ylabel('\epsilon_{C}/\epsilon_{Be}','fontsize',10)
xlim([0 30])
set(gca,'fontsize',10)
ylim([0.4 1.1])


%%
%GI model
addpath('./Box_ratioplot')

%Plot actual erosion rates model
subplot(5,1,2)
thresh=d18O_triang(tdegla+1);
x=[0]
y=[new_erosion_rate]
for i=1:length(d18O_triang)-1
    if(d18O_triang(i)>=thresh && d18O_triang(i+1)<thresh)
        x=[x ti(i)*1e3 ti(i)*1e3];
        y=[y old_erosion_rate new_erosion_rate];
    elseif(d18O_triang(i)<=thresh && d18O_triang(i+1)>thresh)
        x=[x ti(i)*1e3 ti(i)*1e3];
        y=[y new_erosion_rate old_erosion_rate];
    end
end
plot(x,1E6*y,'linewidth',1)
set(gca,'fontsize',10,'yscale','log')
ylabel('\epsilon [mm/kyr]','fontsize',10)
xlim([0 200])
ylim([1E6*old_erosion_rate/3 1E6*new_erosion_rate*3])
set(gca,'Xtick',[0 50 100 150 200],'Ytick',[10.000000001 100])


%%
%Plot calculated apparent erosion rates for the GI model
fs.d18O_filename = 'lisiecki_triinterp_2p6Ma_5ky.mat';
d18Ofn = fs.d18O_filename;
load(d18Ofn); %must contain a variable d18O_triang, sampled in steps of 0.001 Ma
if ~exist('d18O_triang','var')
    error(['the filename ',d18Ofn,' did not contain the variable d18O_triang'])
end
fs.ti=ti;
fs.d18O_triang=d18O_triang;

fs.zobs=0;

fs.production.P14_top_spal=fs.P14_top_spal;
fs.production.P10_top_spal=fs.P10_top_spal;
fs.production.P26_top_spal=fs.P26_top_spal;

fs.production.P14_top_nm=fs.P14_top_nm;
fs.production.P10_top_nm=fs.P10_top_nm;
fs.production.P26_top_nm=fs.P26_top_nm;

fs.production.P14_top_fm=fs.P14_top_fm;
fs.production.P10_top_fm=fs.P10_top_fm;
fs.production.P26_top_fm=fs.P26_top_fm;

original=fs.d18O_triang;

fs.m_true(2)=old_erosion_rate;
fs.m_true(1)=new_erosion_rate;

fs.m_true(3)=tdegla*1000;

fs.m_true(4)=fs.d18O_triang(floor(tdegla)+1)+1e-5; %The 1e-5 is deducted to make

start=2000;
for i=start:length(original)-1
    i
    clear fs.d18O_triang
    fs.d18O_triang=ones(size(original))*original(end);
    fs.d18O_triang(1:i)=original(length(original)-i:end-1);
    
    [d_true,fs.lump_m_true]= g(fs.m_true,fs);
    e_C10(i-start+1)=bisection_e_app('10Be',d_true(1),fs);
    e_C14(i-start+1)=bisection_e_app('14C',d_true(2),fs);
    
    ratio_list_faster_BeC(i-start+1)=e_C14(i-start+1)/e_C10(i-start+1);
end
hold on

%Plot the apparent erosion rates
plot(0:1:length(e_C10)-1,1E6*flip(e_C10),'r','linewidth',1)
plot(0:1:length(e_C10)-1,1E6*flip(e_C14),'g','linewidth',1)
j=legend('\epsilon','\epsilon_{Be}','\epsilon_{C}');
set(j,'fontsize',10,'location','northeast')

%%
%Plot the ratio of apparent erosion rates for the GI model
subplot(5,1,4)
semilogy(0:1:length(e_C10)-1,flip(ratio_list_faster_BeC),'k','linewidth',1)
ylabel('\epsilon_{C} / \epsilon_{Be}','fontsize',10)

%%
%Plot forward GC model
subplot(5,1,3)

kvar=-0.86;
mean_ero=5*old_erosion_rate;
min(mean_ero*(kvar*(d18O_triang-mean(d18O_triang(1:135)))+1));
max(mean_ero*(kvar*(d18O_triang-mean(d18O_triang(1:135)))+1));
fs.d18O_triang=(fs.d18O_triang-(max(fs.d18O_triang)+min(fs.d18O_triang))/2)/((max(fs.d18O_triang)-min(fs.d18O_triang))/2);
original=fs.d18O_triang;
plot(ti*1000,1E6*mean_ero*(kvar*(fs.d18O_triang)+1),'linewidth',1)
set(gca,'fontsize',10,'yscale','log')
ylabel('\epsilon [mm/kyr]','fontsize',10)
xlim([0 200])
ylim([1E6*old_erosion_rate/3 1E6*new_erosion_rate*3])
set(gca,'Xtick',[0 50 100 150 200],'Ytick',[10.000000001 100])

%%
%Calculate apparent erosion rates for the GC model
rmpath('./Box_ratioplot')

addpath('../Figure8/GC')

fs.d18O_min=min(fs.d18O_triang);
fs.d18O_max=max(fs.d18O_triang);
fs.d18O_std=std(fs.d18O_triang);
fs.d18O_mean=0*mean(fs.d18O_triang(1:135));
fs.Nucleides = {'10Be','14C'}; %We may switch nucleides on and off

fs.m_true(2)=mean_ero;
fs.m_true(1)=kvar;

clear e_C14
clear e_C10

start=2400;
for i=start:length(original)-1
    i
    clear fs.d18O_triang
    fs.d18O_triang=ones(size(original))*original(end);
    fs.d18O_triang(1:i)=original(length(original)-i:end-1);
    
    [d_true,fs.lump_m_true]= g_forward(fs.m_true,fs);
    e_C14(i-start+1)=bisection_e_app('14C',d_true(2),fs);
    e_C10(i-start+1)=bisection_e_app('10Be',d_true(1),fs);
    
    ratio_list_slower(i-start+1)=e_C14(i-start+1)/e_C10(i-start+1);
end

hold on
plot(0:1:length(e_C10)-1,1E6*flip(e_C10),'r','linewidth',1)
plot(0:1:length(e_C10)-1,1E6*flip(e_C14),'g','linewidth',1)

%%
%Plot apparent erosion ratio for the GC model
subplot(5,1,4)
semilogy(0:1:length(e_C10)-1,flip(ratio_list_slower),'c','linewidth',1)


ylabel('\epsilon_{C} / \epsilon_{Be}','fontsize',10)
j=legend('SC','GI','GC');
set(j,'fontsize',10,'Location','northeast')
ylim([0.4 1.5])
ylim([0.2 5])
xlim([0 200])
set(gca,'fontsize',10)
set(gca,'Xtick',[0 50 100 150 200],'Ytick',[0.2 0.5 1 2 5])

%%
%Plot the normalized d18O curve
subplot(5,1,5)
plot(0:1:length(e_C10)-1,fs.d18O_triang(1:201),'b','linewidth',1)
xlabel('Time before present [kyr]','fontsize',10)
ylabel(['\delta^{18}O [',char(8240),']'],'fontsize',10)


%%
%Print figure
print_string=['../Figures/Figure4_plot_forward']
save='y'
if save=='y'
    export_fig(print_string,'-tiff','-r1000','-transparent')
    export_fig(print_string,'-pdf','-r3000','-transparent')
end

rmpath('../Figure9/Variable_ratioplot')
rmpath('../Figure9/Variable_ratioplot/Functions')