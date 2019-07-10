function dy = tov(r, y)
new_zeta = @(xi) barotropic(xi) - log10(y(2));
log_rho = fzero(new_zeta, y(2));
rho = 10^log_rho;
dy(1) = 4*pi*r^2*rho;
dy(2) = (-y(1)*rho/r^2*(1 + y(2) / rho)*(1 + 4*pi*r^3*y(2)/y(1))/(1-2*y(1)/r));
end