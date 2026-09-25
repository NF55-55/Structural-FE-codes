function [M,Fint,acck1,ordRow,ordCol,ind] = rearrangeUnknowns(M,F,d,Fint,acck1)

% re-organaizing the rows so that we end up with Fi(known) on the upper part:

ordRow = 1:length(F); % to be able to re-order everything afterwards

ordRow = [ordRow(~isnan(F)), ordRow(isnan(F))];
M = [M(~isnan(F),:); M(isnan(F),:)];
Fint = [Fint(~isnan(F)); Fint(isnan(F))];
ind = sum(~isnan(F));

% re-organaizing the col. so that we end up with ui(unknown) on the upper part of ui:

ordCol = 1:length(d); % to be able to re-order everything afterwards

ordCol = [ordCol(isnan(d)), ordCol(~isnan(d))];
acck1 = [acck1(isnan(d)); acck1(~isnan(d))];
M = [M(:,isnan(d)), M(:,~isnan(d))];

end

