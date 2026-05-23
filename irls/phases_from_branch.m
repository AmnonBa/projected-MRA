function phi = phases_from_branch(beta, eps, p)
    q = numel(beta);
    rhs = 2*sum(eps(1:q-1).*beta(1:q-1)) + eps(q)*beta(q);
    phi1 = rhs / p;

    phi = zeros(q,1);
    phi(1) = phi1;

    for k = 2:q
        phi(k) = k*phi1 - sum(eps(1:k-1).*beta(1:k-1));
    end
end
