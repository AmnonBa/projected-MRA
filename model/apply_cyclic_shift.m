function shifted = apply_cyclic_shift(theta, ell)
    shifted = circshift(theta(:), ell);
end
