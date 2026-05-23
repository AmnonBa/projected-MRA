function Pmat = projection_shift_matrix(p, ell)
    q = (p - 1)/2;
    Pmat = zeros(q,p);

    for j = 1:q
        mPlus = j;
        mMinus = mod(-j,p);

        sourcePlus = mod(mPlus - ell, p);
        sourceMinus = mod(mMinus - ell, p);

        Pmat(j, sourcePlus + 1) = Pmat(j, sourcePlus + 1) + 1;
        Pmat(j, sourceMinus + 1) = Pmat(j, sourceMinus + 1) + 1;
    end
end

