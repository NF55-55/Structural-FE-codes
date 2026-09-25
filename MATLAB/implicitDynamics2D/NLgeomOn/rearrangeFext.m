function Fext = rearrangeFext(F,Fext)

% re-organaizing the rows so that we end up with Fi(known) on the upper part:
Fext = [Fext(~isnan(F)); Fext(isnan(F))];

end

