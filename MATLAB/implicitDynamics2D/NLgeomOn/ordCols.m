function [ordCol,ui] = ordCols(ordCol,ui)

if sum(int16([size(ui)]==1))>=1
    ui(ordCol) = ui;
else % is ui is a matrix
    ui(:,ordCol) = ui;
end

ordCol(ordCol) = ordCol;

return

