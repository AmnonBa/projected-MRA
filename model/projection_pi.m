function xProj = projection_pi(v)
    p = numel(v);
    q = (p - 1)/2;
    xProj = zeros(q,1);

    for j = 1:q
        idxPlus = j + 1;
        idxMinus = p - j + 1;
        xProj(j) = v(idxPlus) + v(idxMinus);
    end
end
