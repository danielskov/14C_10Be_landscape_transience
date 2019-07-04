function [d,lump] = g(m,fixed_stuff,iwalk)
ErateInt = m(1);
ErateGla = m(2);
tDegla   = m(3);
d18Oth   = m(4);
d18Ofn = fixed_stuff.d18O_filename;

i_t_trunc = find(fixed_stuff.ti*1e6>20e3);
tExtract = fixed_stuff.ti(i_t_trunc)*1e6;
[tStarts,relExpos] = ExtractHistory2(tExtract,fixed_stuff.d18O_triang(i_t_trunc),d18Oth,ErateInt,ErateGla); %should fix the ys=midvalue problem

fixed_stuff.tStarts = tStarts;
fixed_stuff.relExpos = relExpos;

[cBes,cCs,lump] = ...
    gCosmoLongsteps(ErateInt,ErateGla,tDegla,fixed_stuff);

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



