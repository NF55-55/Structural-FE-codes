function [d,F,mapVect] = rearrangeUnknowns(d,F)

% re-organaizing the rows so that we end up with Fi(known) on the upper part:

mapVect = 1:length(F); % to be able to re-order everything afterwards

mapVect = [mapVect(isnan(d)), mapVect(~isnan(d))]';

F = F(mapVect);
d = d(mapVect);

end

