function [ratio_list_faster,e_new_fractions_faster,e_old_list,ratio_list_faster_10pct,ratio_list_faster_5pct,ratio_list_faster_25pct,ratio_list_slower_10pct,ratio_list_slower_5pct,ratio_list_slower_25pct]=SC_test(fs,times,ero_min,ero_max,n)
%This file plots a surface showing the ratio between apparent erosion rates
%for different nuclides at a given time.

density=fs.density;
decay_factor_C=fs.L14;
decay_factor_Be=fs.L10;

e_old_list=logspace(log10(ero_min),log10(ero_max),n)*density; %The vector defining the erosion rates prior to a change
e_new_fractions_faster=logspace(2.4,-2.4,n);

for i=1:length(e_old_list) %looping all the erosion rates prior to the change
    for j=1:length(e_new_fractions_faster) %looping all the erosion rates prior to the change
        e_old=e_old_list(i);    
        e_new_faster=e_old*e_new_fractions_faster(j); %The vector defining the new erosion rates is calculated
        
        e_new_faster=e_old_list(i);
        e_old=e_new_faster/e_new_fractions_faster(j);
        
        
        %10Be
        C_ss=C_ss_calculator('10Be',e_new_faster*times,e_old,fs);
        
        C_spal=fs.P10_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.att_l_spal)))./(e_new_faster+fs.att_l_spal*decay_factor_Be);
        C_nm=fs.P10_top_nm*fs.P10_att_l_nm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_nm)))./(e_new_faster+fs.P10_att_l_nm*decay_factor_Be);
        C_fm=fs.P10_top_fm*fs.P10_att_l_fm.*(1-exp(-times.*(decay_factor_Be+e_new_faster/fs.P10_att_l_fm)))./(e_new_faster+fs.P10_att_l_fm*decay_factor_Be);
        
        C_content=C_ss.*exp(-times*decay_factor_Be)+C_spal+C_nm+C_fm;
        
        e_app_Be=bisection_e_app('10Be',C_content,fs);
        e_app_Be09=bisection_e_app('10Be',C_content*0.9,fs);
        e_app_Be095=bisection_e_app('10Be',C_content*0.95,fs);
        e_app_Be098=bisection_e_app('10Be',C_content*0.975,fs);
        
        e_app_Be11=bisection_e_app('10Be',C_content*1.1,fs);
        e_app_Be105=bisection_e_app('10Be',C_content*1.05,fs);
        e_app_Be102=bisection_e_app('10Be',C_content*1.025,fs);
        
        %14C        
        clear C_ss C_spal C_nm C_fm C_content
                
        C_ss=C_ss_calculator('14C',e_new_faster*times,e_old,fs);      
        
        C_spal=fs.P14_top_spal*fs.att_l_spal.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.att_l_spal)))./(e_new_faster+fs.att_l_spal*decay_factor_C);
        C_nm=fs.P14_top_nm*fs.P14_att_l_nm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_nm)))./(e_new_faster+fs.P14_att_l_nm*decay_factor_C);
        C_fm=fs.P14_top_fm*fs.P14_att_l_fm.*(1-exp(-times.*(decay_factor_C+e_new_faster/fs.P14_att_l_fm)))./(e_new_faster+fs.P14_att_l_fm*decay_factor_C);
                
        C_content=C_ss.*exp(-times*decay_factor_C)+C_spal+C_nm+C_fm;
        
        e_app_C=bisection_e_app('14C',C_content,fs);
        e_app_C11=bisection_e_app('14C',C_content*1.10,fs);
        e_app_C105=bisection_e_app('14C',C_content*1.05,fs);
        e_app_C102=bisection_e_app('14C',C_content*1.025,fs);
        
        e_app_C09=bisection_e_app('14C',C_content*0.9,fs);
        e_app_C095=bisection_e_app('14C',C_content*0.95,fs);
        e_app_C098=bisection_e_app('14C',C_content*0.975,fs);
        
        
        %Calculating ratio
        ratio_list_faster(j,i)=e_app_C/e_app_Be;
        ratio_list_faster_10pct(j,i)=e_app_C11/e_app_Be09;
        ratio_list_faster_5pct(j,i)=e_app_C105/e_app_Be095;
        ratio_list_faster_25pct(j,i)=e_app_C102/e_app_Be098;
        
        ratio_list_slower_10pct(j,i)=e_app_C09/e_app_Be11;
        ratio_list_slower_5pct(j,i)=e_app_C095/e_app_Be105;
        ratio_list_slower_25pct(j,i)=e_app_C098/e_app_Be102;
    end
end
   
return