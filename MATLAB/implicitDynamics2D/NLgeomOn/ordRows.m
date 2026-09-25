function [ordRow,Fi] = ordRows(ordRow,Fi)

if sum(int16([size(Fi)]==1))>=1
    Fi(ordRow) = Fi;
else % if Fi is a matrix
    Fi(ordRow,:) = Fi;
end

ordRow(ordRow) = ordRow;

return
