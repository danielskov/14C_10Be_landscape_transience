function erosion_rate=bisection_e_app(nuclide,C_content,fs)
%Daniel S. Skov
%December 2018
%Returns the apparent erosion rate in m/yr
%Calculated with bisection method

att_l_spal=fs.att_l_spal;
density=fs.density;

switch nuclide
    case '10Be'
        half_life=fs.L10;
        P_spal=fs.P10_top_spal;
        P_nm=fs.P10_top_nm;
        P_fm=fs.P10_top_fm;
        att_l_nm=fs.P10_att_l_nm;
        att_l_fm=fs.P10_att_l_fm;
    case '26Al'
        half_life=fs.L26;
        P_spal=fs.P26_top_spal;
        P_nm=fs.P26_top_nm;
        P_fm=fs.P26_top_fm;
        att_l_nm=fs.P26_att_l_nm;
        att_l_fm=fs.P26_att_l_fm;
    case '14C'
        half_life=fs.L14;
        P_spal=fs.P14_top_spal;
        P_nm=fs.P14_top_nm;
        P_fm=fs.P14_top_fm;
        att_l_nm=fs.P14_att_l_nm;
        att_l_fm=fs.P14_att_l_fm;
    case '21Ne'
        half_life=fs.L21;
        P_spal=fs.P21_top_spal;
        P_nm=fs.P21_top_nm;
        P_fm=fs.P21_top_fm;
        att_l_nm=fs.P21_att_l_nm;
        att_l_fm=fs.P21_att_l_fm;
    otherwise
        disp('Wrong nuclide chosen for bisection method')
        return
end


%for ai=1:100*100
lower_limit=1e-8;
upper_limit=1e1;


low_lim=P_spal*att_l_spal./(lower_limit+att_l_spal*half_life)+P_nm*att_l_nm./(lower_limit+att_l_nm*half_life)+P_fm*att_l_fm./(lower_limit+att_l_fm*half_life)-C_content;
up_lim=P_spal*att_l_spal./(upper_limit+att_l_spal*half_life)+P_nm*att_l_nm./(upper_limit+att_l_nm*half_life)+P_fm*att_l_fm./(upper_limit+att_l_fm*half_life)-C_content;

n_limit=10;
n1=1;
while(up_lim>0 && n1<n_limit)
    upper_limit=upper_limit*10;
    up_lim=P_spal*att_l_spal./(upper_limit+att_l_spal*half_life)+P_nm*att_l_nm./(upper_limit+att_l_nm*half_life)+P_fm*att_l_fm./(upper_limit+att_l_fm*half_life)-C_content;
    n1=n1+1;
end
n2=1;
while(low_lim<0 && n2<n_limit)
    lower_limit=lower_limit/10;
    low_lim=P_spal*att_l_spal./(lower_limit+att_l_spal*half_life)+P_nm*att_l_nm./(lower_limit+att_l_nm*half_life)+P_fm*att_l_fm./(lower_limit+att_l_fm*half_life)-C_content;
    n2=n2+1;
end
if(n2==n_limit)
    erosion_rate=1E-9;
    if(strcmp(nuclide,'10Be'))
        erosion_rate=1E-10;
    end
    if(strcmp(nuclide,'26Al'))
        erosion_rate=1E-11;
    end
    erosion_rate=nan;
    return
elseif(n1==n_limit)
    erosion_rate=1E18;
    disp('upper reached')
    disp('nuclide')
    nuclide
    erosion_rate=nan;
    return
end

n=50;
i=0;
a=1e8;

if(sign(low_lim)~=sign(up_lim))
    while(abs(a)>1E-3 && i<n)
        i=i+1;
        %for i=1:n
        
        mid=10^((log10(lower_limit)+log10(upper_limit))/2);

        if(i>1)
            if((mid-upper_limit)==0)
                %disp('break')
                break
            elseif((mid-lower_limit)==0)
                %disp('break')
                break
            end
        end
        
        %midd(i)=mid;
        
        a=P_spal*att_l_spal./(mid+att_l_spal*half_life)+P_nm*att_l_nm./(mid+att_l_nm*half_life)+P_fm*att_l_fm./(mid+att_l_fm*half_life)-C_content;
        %aa(i)=a;
        if(sign(low_lim)==sign(a))
            lower_limit=mid;
        else
            upper_limit=mid;
        end
        
        if(a==0)
            break
        end
    

        
    end
end

if(exist('mid','var')==0)
    sign(low_lim)
    sign(up_lim)
    figure
    
    erosion_rates=logspace(-8,1,1e5);
    
    a=P_spal*att_l_spal./(erosion_rates+att_l_spal*half_life)+P_nm*att_l_nm./(erosion_rates+att_l_nm*half_life)+P_fm*att_l_fm./(erosion_rates+att_l_fm*half_life)-C_content;
    semilogx(erosion_rates,a)
    low_lim=P_spal*att_l_spal./(lower_limit+att_l_spal*half_life)+P_nm*att_l_nm./(lower_limit+att_l_nm*half_life)+P_fm*att_l_fm./(lower_limit+att_l_fm*half_life)-C_content
    up_lim=P_spal*att_l_spal./(upper_limit+att_l_spal*half_life)+P_nm*att_l_nm./(upper_limit+att_l_nm*half_life)+P_fm*att_l_fm./(upper_limit+att_l_fm*half_life)-C_content
    disp('n1')
    n1
    disp('n2')
    n2
    nuclide
end
erosion_rate=mid/density;
return