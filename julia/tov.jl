#=
tov:
- Julia version: 1.0.3
- Author: kavis
- Date: 2019-01-06
=#

using DifferentialEquations
using Plots

const K = 100;
const Γ = 2;
const ρ_c = 1.28E-3;

# TOV
function tov_rhs(r, m_old, Φ_old, p_old)
    r2 = r^2;

    ρ_old = rho_of_p(p_old);

    m_r = 0.;
    m_rr = 0.;
    m_rrr = 0.;

    if r < 1E-3
        m_rrr = 4/3*π*(ρ_c^2/(Γ - 1) + 4^(-1/Γ)*ρ_c^(2/Γ));
        m_rr = m_rrr * r;
        m_r = m_rr * r;
    else
        m_r = m_old / r;
        m_rr = m_r / r;
        m_rrr = m_rr / r;
    end

    m_new = 4*π*r2*ρ_old;
    p_new = -ρ_old*m_rr*(1 + p_old/ρ_old)*(1 + 4*π*p_old/m_rrr)*1/(1 - 2*m_r);
    Φ_new = -1/ρ_old * p_new / (1 + p_old/ρ_old);

    return m_new, Φ_new, p_new
end

function rho_of_p(p)
    if p < 0
        p = abs(p)
    end
    return (p/K)^(1/Γ) + p / (Γ - 1)
end

f = @ode_def TOV begin
    dm = 4*π*t^2*rho_of_p(p)
    dΦ = -1/rho_of_p(p) * p / (1 + p/rho_of_p(p))
    dp = -rho_of_p(p)*m/t^2 * (1 + p/rho_of_p(p)) * (1 + 4*π*p*t^3/m) / (1 - 2*m/t)
end

p_c = K * ρ_c^(Γ)
u0 = [0;0;p_c]
tspan = (0.0, 16.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob)
plot(sol)