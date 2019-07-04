%This function is a forward model for calculating the cosmogenic nuclide
%concentration given erosion rates that scale with the d18O value as
%e=k1*(d18O-d18O_thres)+k0, where k0 is a mean erosion rate, and k1 is scaled
%such that erosion rates are always positive.


function [d,lump] = g_forward(m,fixed_stuff)

k2_meanerosion=m(2);
k1_varerosion=m(1);
d18Ofn = fixed_stuff.d18O_filename;
ti=fixed_stuff.ti*1e6;
d18O_std=fixed_stuff.d18O_std;
d18O_min=fixed_stuff.d18O_min;
d18O_max=fixed_stuff.d18O_max;
d18O_mean=fixed_stuff.d18O_mean;

[erosion_rates,depths,relExpos,interval_time] = erosion_rates_depths_exposure(ti,fixed_stuff.d18O_triang,[],[],k1_varerosion,k2_meanerosion,d18O_mean,d18O_min,d18O_max,d18O_std);
erosion_rates=flipud(erosion_rates);
relExpos=flipud(relExpos);
relExpos=relExpos(:)';
erosion_rates=erosion_rates(:)'; %force row
fixed_stuff.tStarts = -flipud(ti);
fixed_stuff.relExpos = relExpos;

[cBes,cAls,cNes,cCs,lump] = ...
    gCosmoLongsteps_forward([],[],[],[],[],fixed_stuff,erosion_rates);
lump.erosion_rates=erosion_rates;

d = [];
for iNucl = 1:length(fixed_stuff.Nucleides)
    switch fixed_stuff.Nucleides{iNucl}
        case '10Be', d=[d;cBes(:)];
        case '26Al', d=[d;cAls(:)];
        case '21Ne', d=[d;cNes(:)];
        case '14C', d=[d;cCs(:)];
        case '26Al/10Be',d=[d;cAls(:)./cBes(:)];
        case '21Ne/10Be',d=[d;cNes(:)./cBes(:)];
        otherwise
            error([fixed_stuff.Nucleides{iNucl},' not implemented'])
    end
end