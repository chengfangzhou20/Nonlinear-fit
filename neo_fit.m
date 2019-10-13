%neo-hookean
function sigma = neo_fit(k_elastin,lam)
sigma  = lam.^2 .* k_elastin .* (1 - (1./lam.^6));
end