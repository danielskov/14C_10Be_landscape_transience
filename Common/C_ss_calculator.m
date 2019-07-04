function C_ss=C_ss_calculator(nuclide,depth,e_old,fs)
%Calculates the steady state CN concentration

att_l_spal=fs.att_l_spal;

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

C_spal=exp(-depth/att_l_spal)*P_spal*att_l_spal/(e_old+att_l_spal*half_life);
C_nm=exp(-depth/att_l_nm)*P_nm*att_l_nm/(e_old+att_l_nm*half_life);
C_fm=exp(-depth/att_l_fm)*P_fm*att_l_fm/(e_old+att_l_fm*half_life);

C_ss=C_spal+C_nm+C_fm;
return


