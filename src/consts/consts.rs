use self::Dimension::*;

/// Fundamental Constants with various unit systems
#[derive(Debug, Copy, Clone)]
pub struct Units {
    pub c: f64,
    pub G: f64,
    pub e: f64,
    pub k_b: f64,
    pub N_A: f64,
    pub h: f64,
    pub m_u: f64,
}

pub enum Dimension {
    Time,
    Length,
    Mass,
    Velocity,
    AngularVelocity,
    Acceleration,
    Energy,
    EnergyDensity,
    AngularMomentum,
    Force,
    Power,
    Pressure,
    Density,
}

/// CGS Unit
pub const CGS: Units = Units {
    c: 2.99792458e+10,   // Speed of light (cm s^{-1})
    G: 6.6742e-8,        // Gravitational constants (cm^3 g^{-1} s^{-2})
    e: 1.60317653e-19,   // Elementary charge (C)
    k_b: 1.3806505e-16,  // Boltzmann constant (erg K^{-1})
    N_A: 6.0221415e+23,  // Avogadro constant (mol^{-1})
    h: 6.6260693e-27,    // Planck constant (erg s)
    m_u: 1.66053886e-24, // Atomic mass unit (g)
};

pub fn cgs_to_geom(value: f64, dim: Dimension) -> f64 {
    match dim {
        Time => value * CGS.c,
        Length => value,
        Mass => value * CGS.G / CGS.c.powi(2),
        Velocity => value / CGS.c,
        AngularVelocity => value / CGS.c,
        Acceleration => value / CGS.c.powi(2),
        Energy => value * CGS.G / CGS.c.powi(4),
        EnergyDensity => value * CGS.G / CGS.c.powi(4),
        AngularMomentum => value * CGS.G / CGS.c.powi(3),
        Force => value * CGS.G / CGS.c.powi(4),
        Power => value * CGS.G / CGS.c.powi(5),
        Pressure => value * CGS.G / CGS.c.powi(4),
        Density => value * CGS.G / CGS.c.powi(2),
    }
}