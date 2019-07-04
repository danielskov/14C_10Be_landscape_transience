function [glacial_erosionrates,fractions_faster,ratio_list_faster_BeC,ratio_list_faster_BeC_102,ratio_list_faster_BeC_105,ratio_list_faster_BeC_11,ratio_list_slower_BeC_102,ratio_list_slower_BeC_105,ratio_list_slower_BeC_11]=GI_tester(fs,time,ero_min,ero_max,n);

fs.d18O_filename = 'lisiecki_triinterp_2p6Ma_5ky.mat';
d18Ofn = fs.d18O_filename;
load(d18Ofn); %must contain a variable d18O_triang, sampled in steps of 0.001 Ma
if ~exist('d18O_triang','var')
    error(['the filename ',d18Ofn,' did not contain the variable d18O_triang'])
end
fs.ti=ti;
fs.d18O_triang=d18O_triang;

fs.Nucleides = {'10Be','14C'}; %We may switch nucleides on and off
fs.zobs=0;

% Rock density in kg/m3
rho = fs.density;

%Half lives
L10=fs.L10;
L26=fs.L26;
L14=fs.L14;

fs.production.P14_top_spal=fs.P14_top_spal;
fs.production.P10_top_spal=fs.P10_top_spal;

fs.production.P14_top_nm=fs.P14_top_nm;
fs.production.P10_top_nm=fs.P10_top_nm;

fs.production.P14_top_fm=fs.P14_top_fm;
fs.production.P10_top_fm=fs.P10_top_fm;

glacial_erosionrates=logspace(log10(ero_min),log10(ero_max),n);
fractions_faster=linspace(1,11,n);
fractions_faster=logspace(2.1,-2.1,n);

ratio_list_faster_BeC=zeros(length(glacial_erosionrates),length(fractions_faster));

fs.m_true(3)=time;

fs.m_true(4)=fs.d18O_triang(floor(fs.m_true(3)./1e3)+1)+1e-5; %The 1e-5 is added to force a shift at the time of deglaciation

for i=1:length(glacial_erosionrates)
    
    fs.m_true(1)=glacial_erosionrates(i);
    if(mod(i,10)==0)
        i
    end
    for ii=1:length(fractions_faster)
        fs.m_true(2)=glacial_erosionrates(i)/fractions_faster(ii);
        [d_true,fs.lump_m_true]= g(fs.m_true,fs);
        
        %10Be
        e_C10=bisection_e_app('10Be',d_true(1),fs);
        
        e_app_Be0975=bisection_e_app('10Be',d_true(1)*0.975,fs);
        e_app_Be095=bisection_e_app('10Be',d_true(1)*0.95,fs);
        e_app_Be09=bisection_e_app('10Be',d_true(1)*0.9,fs);
        
        e_app_Be102=bisection_e_app('10Be',d_true(1)*1.025,fs);
        e_app_Be105=bisection_e_app('10Be',d_true(1)*1.05,fs);
        e_app_Be11=bisection_e_app('10Be',d_true(1)*1.1,fs);
        
        
        %14C
        e_C14=bisection_e_app('14C',d_true(2),fs);

        e_app_C1025=bisection_e_app('14C',d_true(2)*1.025,fs);
        e_app_C105=bisection_e_app('14C',d_true(2)*1.05,fs);
        e_app_C11=bisection_e_app('14C',d_true(2)*1.1,fs);
        
        e_app_C098=bisection_e_app('14C',d_true(2)*0.975,fs);
        e_app_C095=bisection_e_app('14C',d_true(2)*0.95,fs);
        e_app_C09=bisection_e_app('14C',d_true(2)*0.9,fs);
        
        ratio_list_faster_BeC(ii,i)=e_C14/e_C10;
        
        ratio_list_faster_BeC_102(ii,i)=e_app_C1025/e_app_Be0975;
        ratio_list_faster_BeC_105(ii,i)=e_app_C105/e_app_Be095;
        ratio_list_faster_BeC_11(ii,i)=e_app_C11/e_app_Be09;
        
        ratio_list_slower_BeC_102(ii,i)=e_app_C098/e_app_Be102;
        ratio_list_slower_BeC_105(ii,i)=e_app_C095/e_app_Be105;
        ratio_list_slower_BeC_11(ii,i)=e_app_C09/e_app_Be11;
    end
end
return