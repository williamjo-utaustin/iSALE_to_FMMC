function N = particle2grain(m_sim_kg, d_um, rho_kg_m3)
%GRAINS_PER_SIM_PARTICLE  Number of 10-um grains represented by a sim particle.
%   N = GRAINS_PER_SIM_PARTICLE(m_sim_kg)               % d=10 um, rho=3100 kg/m^3
%   N = GRAINS_PER_SIM_PARTICLE(m_sim_kg, d_um)
%   N = GRAINS_PER_SIM_PARTICLE(m_sim_kg, d_um, rho_kg_m3)
%
%   m_sim_kg can be scalar or a vector/array of particle masses (kg).
%   d_um is grain DIAMETER in micrometers (default 10).
%   rho_kg_m3 is grain (0% porosity) density in kg/m^3 (default 3100).

    if nargin < 2 || isempty(d_um),       d_um = 10;       end
    if nargin < 3 || isempty(rho_kg_m3),  rho_kg_m3 = 3100; end

    r_m   = (d_um * 1e-6) / 2;                         % radius in meters
    m_grain = rho_kg_m3 * (4/3) * pi * r_m.^3;         % kg per real grain
    N = m_sim_kg ./ m_grain;
end
