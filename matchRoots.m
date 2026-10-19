function ordered = matchRoots(previous, current)

    P = perms(1:numel(current));

    cost = zeros(size(P,1),1);

    for m = 1:size(P,1)
        cost(m) = sum(abs(previous - current(P(m,:))));
    end

    [~, ind] = min(cost);

    ordered = current(P(ind,:));
end