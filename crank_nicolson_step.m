function psi_new = crank_nicolson_step(psi, A, B)
rhs = B * psi;
psi_new = A \ rhs;   %psi_new = inv(A) * rhs;
end
