function [dist, accuracy] = local_ancestry_evaluate(trueZ, predZ, trueLoc, predLoc, Param)

%local_ancestry_evaluate    evaluate local ancestry prediction
%  trueZ is the true ancestry assignment
%  predZ is the predicted ancestry assignment
%  trueLoc is nAncestry * 2 matrix
%  predLoc is nAncestry * 2 matrix
%  the trueLoc and predLoc are assumed to be aligned

if Param.nMaternal == 0
  accuracy = sum(trueZ == predZ) / Param.nSnp;
  dist = zeros(1, Param.nSnp);
  for i = 1:Param.nPaternal
    indi = (trueZ == i);
    for j = 1:Param.nPaternal
      indj = (predZ == j);
      ind = indi & indj;
      dist(ind) = pos2dist(trueLoc(i, 1), trueLoc(i, 2), ...
                           predLoc(j, 1), predLoc(j, 2));
    end
  end
else
  accuracy = sum(trueZ(:) == predZ(:)) / (Param.nSnp * 2);
  dist = zeros(2, Param.nSnp);
  for i = 1:Param.nPaternal
    indi = (trueZ(1, :) == i);
    for j = 1:Param.nPaternal
      indj = (predZ(1, :) == j);
      ind = indi & indj;
      dist(1, ind) = ...
        pos2dist(trueLoc(i, 1), ...
                 trueLoc(i, 2), ...
                 predLoc(j, 1), ...
                 predLoc(j, 2));
    end
  end
  
  for i = 1:Param.nMaternal
    indi = (trueZ(2, :) == i);
    for j = 1:Param.nMaternal
      indj = (predZ(2, :) == j);
      ind = indi & indj;
      dist(2, ind) = ...
        pos2dist(trueLoc(i + Param.nPaternal, 1), ...
                 trueLoc(i + Param.nPaternal, 2), ...
                 predLoc(j + Param.nPaternal, 1), ...
                 predLoc(j + Param.nPaternal, 2));
    end
  end
end

