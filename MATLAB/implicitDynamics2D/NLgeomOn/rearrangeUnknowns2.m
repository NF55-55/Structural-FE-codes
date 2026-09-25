function [KtEff,dk1,ordCol] = rearrangeUnknowns2(F,d,KtEff,dk1)

% re-organaizing the rows so that we end up with Fi(known) on the upper part:

KtEff = [KtEff(~isnan(F),:); KtEff(isnan(F),:)];

% re-organaizing the col. so that we end up with ui(unknown) on the upper part of ui:

ordCol = 1:length(d); % to be able to re-order everything afterwards

ordCol = [ordCol(isnan(d)), ordCol(~isnan(d))];
dk1 = [dk1(isnan(d)); dk1(~isnan(d))];
KtEff = [KtEff(:,isnan(d)), KtEff(:,~isnan(d))];
ind2 = sum(isnan(d));

end

