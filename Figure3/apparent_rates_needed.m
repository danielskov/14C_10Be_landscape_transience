%Calculates the apparent erosion rate ratio of 14C/10Be needed to detect landscape transience for a
%range of 10Be apparent erosion rates
%Daniel S. Skov December 2018

clear; close all
path(pathdef)

addpath('../common')
addpath('../common/export_fig')


fs=constants();
n=2e4;

%%
%Calculate 14C
C=logspace(5,8.4,n*2); %14C Concentrations

for i=1:length(C) %Calculate corresponding apparent erosion rates
    e_app_C(i)=bisection_e_app('14C',C(i),fs);
    e_app_C_plus25(i)=bisection_e_app('14C',C(i)*1.025,fs);
    e_app_C_plus5(i)=bisection_e_app('14C',C(i)*1.05,fs);
    e_app_C_plus10(i)=bisection_e_app('14C',C(i)*1.1,fs);
    e_app_C_minus25(i)=bisection_e_app('14C',C(i)*0.975,fs);
    e_app_C_minus5(i)=bisection_e_app('14C',C(i)*0.95,fs);
    e_app_C_minus10(i)=bisection_e_app('14C',C(i)*0.9,fs);
end

%%
%Calculate 10Be
C=logspace(5,12,n); %Calculate corresponding apparent erosion rates
for i=1:length(C)
    e_app_Be(i)=bisection_e_app('10Be',C(i),fs);
    e_app_Be_plus25(i)=bisection_e_app('10Be',C(i)*1.025,fs);
    e_app_Be_plus5(i)=bisection_e_app('10Be',C(i)*1.05,fs);
    e_app_Be_plus10(i)=bisection_e_app('10Be',C(i)*1.1,fs);
    e_app_Be_minus25(i)=bisection_e_app('10Be',C(i)*0.975,fs);
    e_app_Be_minus5(i)=bisection_e_app('10Be',C(i)*0.95,fs);
    e_app_Be_minus10(i)=bisection_e_app('10Be',C(i)*0.9,fs);
end


%%
%Map when offset if significant for increases in erosion rate
increase_Be_25=[];
increase_C_25=[];
increase_Be_5=[];
increase_C_5=[];
increase_Be_10=[];
increase_C_10=[];

for i=1:length(e_app_Be)-1
    for j=1:length(e_app_C)
        if(e_app_Be_minus25(i)<e_app_C_plus25(j) && e_app_Be_minus25(i)>=e_app_C_plus25(j+1)) %Dececting increase in erosion rate
            increase_Be_25(i)=e_app_Be(i);
            increase_C_25(i)=e_app_C(j);
        end
        if(e_app_Be_minus5(i)<e_app_C_plus5(j) && e_app_Be_minus5(i)>=e_app_C_plus5(j+1)) %Dececting increase in erosion rate
            increase_Be_5(i)=e_app_Be(i);
            increase_C_5(i)=e_app_C(j);
        end
        if(e_app_Be_minus10(i)<e_app_C_plus10(j) && e_app_Be_minus10(i)>=e_app_C_plus10(j+1)) %Dececting increase in erosion rate
            increase_Be_10(i)=e_app_Be(i);
            increase_C_10(i)=e_app_C(j);
        end
    end
end
I=find(isnan(increase_C_10)==0);
increase_C_10=increase_C_10(I);
increase_Be_10=increase_Be_10(I);
I=find(isnan(increase_C_5)==0);
increase_C_5=increase_C_5(I);
increase_Be_5=increase_Be_5(I);
I=find(isnan(increase_C_25)==0);
increase_C_25=increase_C_25(I);
increase_Be_25=increase_Be_25(I);

I=find(isnan(increase_Be_10)==0);
increase_C_10=increase_C_10(I);
increase_Be_10=increase_Be_10(I);
I=find(isnan(increase_Be_5)==0);
increase_C_5=increase_C_5(I);
increase_Be_5=increase_Be_5(I);
I=find(isnan(increase_Be_25)==0);
increase_C_25=increase_C_25(I);
increase_Be_25=increase_Be_25(I);


hh=figure(5)
set(hh,'units','centimeters','position',[0,0,9.00,9.00])

hold on
fill([increase_Be_25*1E6 increase_Be_25(1)*1E6 increase_Be_25(1)*1E6],[increase_C_25./increase_Be_25 200  increase_C_25(1)./increase_Be_25(1)],[0.8 0.8 0.8],'LineStyle','none')
fill([increase_Be_5*1E6 increase_Be_5(1)*1E6 increase_Be_5(1)*1E6],[increase_C_5./increase_Be_5 200  increase_C_5(1)./increase_Be_5(1)],[0.6 0.6 0.6],'LineStyle','none')
fill([increase_Be_10*1E6 increase_Be_10(1)*1E6 increase_Be_10(1)*1E6],[increase_C_10./increase_Be_10 200  increase_C_10(1)./increase_Be_10(1)],[0.4 0.4 0.4],'LineStyle','none')

%%
%Map when offset is significant for decreases in erosion rate
decrease_Be_25=[];
decrease_C_25=[];
decrease_Be_5=[];
decrease_C_5=[];
decrease_Be_10=[];
decrease_C_10=[];
for i=1:length(e_app_Be)-1
    for j=1:length(e_app_C)
        if(e_app_Be_plus25(i)<e_app_C_minus25(j) && e_app_Be_plus25(i)>=e_app_C_minus25(j+1)) %Dececting decrease in erosion rate
            decrease_Be_25(i)=e_app_Be(i);
            decrease_C_25(i)=e_app_C(j+1);
        end
        if(e_app_Be_plus5(i)<e_app_C_minus5(j) && e_app_Be_plus5(i)>=e_app_C_minus5(j+1)) %Dececting decrease in erosion rate
            decrease_Be_5(i)=e_app_Be(i);
            decrease_C_5(i)=e_app_C(j+1);
        end
        if(e_app_Be_plus10(i)<e_app_C_minus10(j) && e_app_Be_plus10(i)>=e_app_C_minus10(j+1)) %Dececting decrease in erosion rate
            decrease_Be_10(i)=e_app_Be(i);
            decrease_C_10(i)=e_app_C(j+1);
        end
    end
end
I=find(isnan(decrease_C_10)==0);
decrease_C_10=decrease_C_10(I);
decrease_Be_10=decrease_Be_10(I);
I=find(isnan(decrease_C_5)==0);
decrease_C_5=decrease_C_5(I);
decrease_Be_5=decrease_Be_5(I);
I=find(isnan(decrease_C_25)==0);
decrease_C_25=decrease_C_25(I);
decrease_Be_25=decrease_Be_25(I);

figure(5)
h1=fill([decrease_Be_25*1E6 decrease_Be_25(1)*1E6 decrease_Be_25(1)*1E6],[decrease_C_25./decrease_Be_25 1/200  decrease_C_25(1)./decrease_Be_25(1)],[0.8 0.8 0.8],'LineStyle','none')
h2=fill([decrease_Be_5*1E6 decrease_Be_5(1)*1E6 decrease_Be_5(1)*1E6],[decrease_C_5./decrease_Be_5 1/200  decrease_C_5(1)./decrease_Be_5(1)],[0.6 0.6 0.6],'LineStyle','none')
h3=fill([decrease_Be_10*1E6 decrease_Be_10(1)*1E6 decrease_Be_10(1)*1E6],[decrease_C_10./decrease_Be_10 1/200  decrease_C_10(1)./decrease_Be_10(1)],[0.4 0.4 0.4],'LineStyle','none')
hold on
set(gca,'xscale','log','yscale','log','fontsize',10)
xlim([1e-6 1e-3])
ylim([0.01 100])
xlabel('\epsilon_{Be} [mm/kyr]','fontsize',10)
ylabel('\epsilon_{C}/\epsilon_{Be}','fontsize',10)
set(gca,'fontsize',10,'ytick',[0.010000001 0.1 1 10 100],'xtick',[0.10000000001 1.000000000001 10.000000000001 100.000000000001 1000.000000000001])
h=legend([h1 h2 h3],'2.5%','5%','10%')
set(h,'fontsize',10)
ylim([0.01 100])
xlim([0.1 1000.000000000001])
grid on


print_string=['../Figures/Figure3_adaption_needed']
save='y'
if save=='y'
    export_fig(print_string,'-jpg','-r1000','-transparent')
end

return
figure(6)
h1=fill([1e-7 1e-7],[1e-7 1e-7],[0.8 0.8 0.8],'LineStyle','none')
h2=fill([1e-7 1e-7],[1e-7 1e-7],[0.6 0.6 0.6],'LineStyle','none')
h3=fill([1e-7 1e-7],[1e-7 1e-7],[0.4 0.4 0.4],'LineStyle','none')
plot([1e-7 1e-7],[1e-7 1e-7])
set(gca,'xscale','log','yscale','log','fontsize',10)
xlim([1e-6 1e-3])
ylim([0.01 100])
xlabel('\epsilon_{Be} [mm/kyr]','fontsize',10)
ylabel('\epsilon_{C}/\epsilon_{Be}','fontsize',10)
set(gca,'fontsize',10,'ytick',[0.010000001 0.1 1 10 100],'xtick',[0.10000000001 1.000000000001 10.000000000001 100.000000000001 1000.000000000001])
h=legend([h1 h2 h3],'2.5%','5%','10%')
set(h,'fontsize',10)
ylim([0.01 100])
xlim([0.1 1000.000000000001])
grid on


    export_fig(print_string,'-pdf','-r1000','-transparent')

