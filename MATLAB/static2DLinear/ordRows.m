function [ordRow,Fi] = ordRows(ordRow,Fi)

j = 1;

for i=1:length(ordRow)
    if ordRow(i) == j
        j = j+1;
    else
        for k=1:(length(ordRow)-i)
            vF = Fi(i+1:end);
            vord = ordRow(i+1:end);

            Fi(i:end) = [vF; Fi(i)];
            ordRow(i:end) = [vord, ordRow(i)];
            if ordRow(i) == j
                j = j+1;
                break;
            end
        end
    end
end

return
