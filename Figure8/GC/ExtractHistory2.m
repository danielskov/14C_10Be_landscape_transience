function [tStarts,relExpos] = ExtractHistory2(ti,d18Oi,d18Oth,ErateInt,ErateGla);
%checks for ys==midvalue and fixes it with a nudge.
%Built from ExtractHistory.m
%First, use the code from fillredblue.m
midvalue = d18Oth; ys = d18Oi;

%Dealing with ys == midvalue
nudge_factor = 1e-10;
izero = find(ys==midvalue);
if midvalue == 0 % take scale from ys
  nudge_step = nudge_factor*std(ys)+nudge_factor*mean(abs(ys));
else %take scale from midvalue
  nudge_step=nudge_factor*midvalue;
end
ys(izero)=ys(izero)+nudge_step;
%Now, we shold have avoided the problem that diff(sign(ys-midvalue)) may be
%1 instead of 0 or 2.
i_cross = find(abs(diff(sign(ys-midvalue)))==2); 
%zero crossing between ys(i_cross) and ys(i_cross+1)

if isempty(i_cross)
  if d18Oi(1)>d18Oth %all glaciated
    tStarts = max(ti);
    relExpos = 0;
  else %no glaciation
    tStarts = max(ti);
    relExpos = 1;
  end
else %we have crossings to interpolate
  N_cross = length(i_cross);
  di_cross = (midvalue-ys(i_cross))./(ys(i_cross+1)-ys(i_cross));
  dt = diff(ti(1:2));
  tStarts = [ti(i_cross(:))+di_cross(:)*dt;ti(end)];
  relExpos(1:length(i_cross)) = (d18Oi(i_cross)<d18Oth); %if true, interglacial, and relExpo=1; 
  relExpos = [relExpos(:);(d18Oi(end)<d18Oth)];
end

tStarts = -flipud(tStarts); %we start at negative times at start of Quaternary
relExpos = flipud(relExpos);

%save relExpos.dat relExpos -ascii