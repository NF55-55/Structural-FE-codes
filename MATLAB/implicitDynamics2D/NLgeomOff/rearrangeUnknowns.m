function [M,K,d,ordRow,ordCol] = rearrangeUnknowns(M,K,d,F)

% re-organaizing the rows so that we end up with Fi(known) on the upper part:

ordRow = 1:length(F); % to be able to re-order everything afterwards

ordRow = [ordRow(~isnan(F)), ordRow(isnan(F))];
M = [M(~isnan(F),:); M(isnan(F),:)];
K = [K(~isnan(F),:); K(isnan(F),:)];
% F = [F(~isnan(F)); F(isnan(F))];


% re-organaizing the col. so that we end up with ui(unknown) on the upper part of ui:

ordCol = 1:length(d); % to be able to re-order everything afterwards

ordCol = [ordCol(isnan(d)), ordCol(~isnan(d))];
M = [M(:,isnan(d)), M(:,~isnan(d))];
K = [K(:,isnan(d)), K(:,~isnan(d))];
d = [d(isnan(d)); d(~isnan(d))];
% v = [v(isnan(d)); v(~isnan(d))];
% a = [a(isnan(d)); a(~isnan(d))];


end

