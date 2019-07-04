function [ratio_list_slower,ratio_list_slower_11,ratio_list_slower_105,ratio_list_slower_1025,ratio_list_faster_11,ratio_list_faster_105,ratio_list_faster_1025,e_now,ka20_over_now_slower2]=GC_tester(fs,ero_min,ero_max,n)

    fs.Evarminmax = [-0.999,0.999]; %Erosion domain
    fs.Emeanminmax = [1e-7,5e-3]; %Erosion domain
    fs.EvarDistr = 'uniform';
    fs.EmeanDistr = 'logunif';

fs.d18O_filename = 'lisiecki_triinterp_2p6Ma_5ky.mat';
fs.tStarts = NaN; %load or compute fixed times of more or less glaciated periods
fs.relExpos = NaN; %load or compute degree of exposure in periods
d18Ofn = fs.d18O_filename;
load(d18Ofn); %must contain a variable d18O_triang, sampled in steps of 0.001 Ma
if ~exist('d18O_triang','var')
    error(['the filename ',d18Ofn,' did not contain the variable d18O_triang'])
end
fs.ti=ti;
fs.d18O_triang=d18O_triang;
fs.d18O_triang_not_scaled=d18O_triang;

    fs.d18O_runav=d18O_runav;
    fs.d18O_min=min(fs.d18O_triang);
    fs.d18O_max=max(fs.d18O_triang);
    
    d18O_triang2=2*(fs.d18O_triang)/(fs.d18O_max-fs.d18O_min);
    fs.d18O_triang=(fs.d18O_triang-(max(fs.d18O_triang)+min(fs.d18O_triang))/2)/((max(fs.d18O_triang)-min(fs.d18O_triang))/2);
    
    fs.d18O_min=min(fs.d18O_triang);
    fs.d18O_max=max(fs.d18O_triang);
    fs.d18O_std=std(fs.d18O_triang);
    fs.d18O_mean=0; 

fs.Nucleides = {'10Be','14C'}; %We may switch nucleides on and off
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

fs.mname{1} ='Evar'; %The variation in erosion rate
    fs.mname{2} ='Emean'; %The mean erosion rate

interglacial_erosionrates=logspace(log10(ero_min),log10(ero_max),n);
fractions_slower=linspace(-0.9999,.9999,n);

ratio_list_slower=zeros(length(interglacial_erosionrates),length(fractions_slower));

    for i=1:length(interglacial_erosionrates)
        fs.m_true(2)=interglacial_erosionrates(i);
        if(mod(i,10)==0)
            i
        end
        for ii=1:length(fractions_slower)
            fs.m_true(1)=fractions_slower(ii);
            [d_true,fs.lump_m_true]= g_forward(fs.m_true,fs);
            e_C14=bisection_e_app('14C',d_true(2),fs);
            
            e_app_C0975=bisection_e_app('14C',d_true(2)*0.975,fs);
        e_app_C095=bisection_e_app('14C',d_true(2)*0.95,fs);
        e_app_C09=bisection_e_app('14C',d_true(2)*0.9,fs);
            
        e_app_C1025=bisection_e_app('14C',d_true(2)*1.025,fs);
        e_app_C105=bisection_e_app('14C',d_true(2)*1.05,fs);
        e_app_C11=bisection_e_app('14C',d_true(2)*1.1,fs);
            
            e_C10=bisection_e_app('10Be',d_true(1),fs);
            
            e_app_Be1025=bisection_e_app('10Be',d_true(1)*1.025,fs);
        e_app_Be105=bisection_e_app('10Be',d_true(1)*1.05,fs);
        e_app_Be11=bisection_e_app('10Be',d_true(1)*1.1,fs);
            
             e_app_Be0975=bisection_e_app('10Be',d_true(1)*0.975,fs);
        e_app_Be095=bisection_e_app('10Be',d_true(1)*0.95,fs);
        e_app_Be09=bisection_e_app('10Be',d_true(1)*0.9,fs);
        
        e_now(ii,i)=fs.m_true(2)*((fs.m_true(1)*(fs.d18O_triang(1))+1));
        
        
        ka20_over_now_slower2(ii,i)=((fs.m_true(1)*(fs.d18O_triang(1))+1)/(fs.m_true(1)*(fs.d18O_triang(21))+1));

            
            ratio_list_slower(ii,i)=e_C14/e_C10;
            ratio_list_slower_1025(ii,i)=e_app_C0975/e_app_Be1025;
            ratio_list_slower_105(ii,i)=e_app_C095/e_app_Be105;
            ratio_list_slower_11(ii,i)=e_app_C09/e_app_Be11;
            
            ratio_list_faster_1025(ii,i)=e_app_C1025/e_app_Be0975;
            ratio_list_faster_105(ii,i)=e_app_C105/e_app_Be095;
            ratio_list_faster_11(ii,i)=e_app_C11/e_app_Be09;
        end
    end
return