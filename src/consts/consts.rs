pub use self::Dimension::*;

pub enum UnitSystem {
    CGS,
    Geometrized,
    Natural
}

/// Fundamental Constants with various unit systems
#[allow(non_snake_case)]
#[derive(Debug, Copy, Clone)]
pub struct Units {
    pub c: f64,
    pub G: f64,
    pub e: f64,
    pub k_b: f64,
    pub N_A: f64,
    pub h: f64,
    pub hbar: f64,
    pub m_u: f64,
    pub eV: f64,
}

pub enum Dimension {
    Time,
    Length,
    Mass,
    Velocity,
    Momentum,
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
    hbar: 1.05457266e-27,// Planck Constant (erg s)
    m_u: 1.66053886e-24, // Atomic mass unit (g)
    eV: 1.60217733e-12,  // Electron Volt (erg)
};

pub fn cgs_to_geom(value: f64, dim: Dimension) -> f64 {
    match dim {
        Time => value * CGS.c,
        Length => value,
        Mass => value * CGS.G / CGS.c.powi(2),
        Velocity => value / CGS.c,
        Momentum => value * CGS.G / CGS.c.powi(3),
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

// R. L. Jaffe, Supplementary Notes for MIT's Quantum Theory Sequence
pub fn cgs_to_natural(value: f64, dim: Dimension) -> f64 {
    match dim {
        Time => value * CGS.eV / CGS.hbar,
        Length => value * CGS.eV / (CGS.hbar * CGS.c),
        Mass => value * CGS.c.powi(2) / CGS.eV,
        Velocity => value / CGS.c,
        Momentum => value * CGS.c / CGS.eV,
        AngularVelocity => unimplemented!(),
        AngularMomentum => value / CGS.hbar,
        Acceleration => value * CGS.hbar / (CGS.c * CGS.eV),
        Energy => value / CGS.eV,
        EnergyDensity => value * (CGS.hbar * CGS.c).powi(3) / CGS.eV.powi(4),
        Pressure => value * (CGS.hbar * CGS.c).powi(3) / CGS.eV.powi(4),
        Force => value * (CGS.hbar * CGS.c) / CGS.eV.powi(2),
        Power => value * CGS.hbar / CGS.eV.powi(2),
        Density => value * CGS.c.powi(5) * CGS.hbar.powi(3) / CGS.eV.powi(4),
    }
}

#[allow(non_snake_case)]
pub fn cgs_to_MeV(value: f64, dim: Dimension) -> f64 {
    match dim {
        Time => value * (CGS.eV * 1e+6) / CGS.hbar,
        Length => value * (CGS.eV * 1e+6) / (CGS.hbar * CGS.c),
        Mass => value * CGS.c.powi(2) / (CGS.eV * 1e+6),
        Velocity => value / CGS.c,
        Momentum => value * CGS.c / (CGS.eV * 1e+6),
        AngularVelocity => unimplemented!(),
        AngularMomentum => value / CGS.hbar,
        Acceleration => value * CGS.hbar / (CGS.c * (CGS.eV * 1e+6)),
        Energy => value / (CGS.eV * 1e+6),
        EnergyDensity => value * (CGS.hbar * CGS.c).powi(3) / (CGS.eV * 1e+6).powi(4),
        Pressure => value * (CGS.hbar * CGS.c).powi(3) / (CGS.eV * 1e+6).powi(4),
        Force => value * (CGS.hbar * CGS.c) / (CGS.eV * 1e+6).powi(2),
        Power => value * CGS.hbar / (CGS.eV * 1e+6).powi(2),
        Density => value * CGS.c.powi(5) * CGS.hbar.powi(3) / (CGS.eV * 1e+6).powi(4),
    }
}

#[allow(non_snake_case)]
pub fn cgs_to_GeV(value: f64, dim: Dimension) -> f64 {
    match dim {
        Time => value * (CGS.eV * 1e+9) / CGS.hbar,
        Length => value * (CGS.eV * 1e+9) / (CGS.hbar * CGS.c),
        Mass => value * CGS.c.powi(2) / (CGS.eV * 1e+9),
        Velocity => value / CGS.c,
        Momentum => value * CGS.c / (CGS.eV * 1e+9),
        AngularVelocity => unimplemented!(),
        AngularMomentum => value / CGS.hbar,
        Acceleration => value * CGS.hbar / (CGS.c * (CGS.eV * 1e+9)),
        Energy => value / (CGS.eV * 1e+9),
        EnergyDensity => value * (CGS.hbar * CGS.c).powi(3) / (CGS.eV * 1e+9).powi(4),
        Pressure => value * (CGS.hbar * CGS.c).powi(3) / (CGS.eV * 1e+9).powi(4),
        Force => value * (CGS.hbar * CGS.c) / (CGS.eV * 1e+9).powi(2),
        Power => value * CGS.hbar / (CGS.eV * 1e+9).powi(2),
        Density => value * CGS.c.powi(5) * CGS.hbar.powi(3) / (CGS.eV * 1e+9).powi(4),
    }
}
