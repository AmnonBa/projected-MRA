function A = cosine_matrix_A(p)
    q = (p - 1)/2;
    A = zeros(q,q);
    for j = 1:q
        for k = 1:q
            A(j,k) = 2*cos(2*pi*j*k/p);
        end
    end
end
