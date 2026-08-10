pub mod particle;
pub mod universe;

pub use particle::{Particle, Reference};
pub use universe::{ConsiderEffects, IgnoreGravityTerms, Universe};

use serde::{Deserialize, Serialize};
pub fn calculate_spin(particles: &mut [Particle]) {
    for (i, particle) in particles.iter_mut().enumerate() {
        assert_ne!(
            particle.moment_of_inertia, 0.,
            "Moment of inertia for particle {i} is zero!"
        );
        particle.update_spin();
    }
}

#[derive(Debug, Default, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Axes {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Axes {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline(always)]
    pub fn from(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    #[inline(always)]
    pub fn add(&mut self, other: &Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }

    #[inline(always)]
    pub fn sub(&mut self, other: &Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }

    #[inline(always)]
    pub fn mul(&mut self, factor: f64) {
        self.x *= factor;
        self.y *= factor;
        self.z *= factor;
    }

    #[inline(always)]
    pub fn div(&mut self, factor: f64) {
        self.x /= factor;
        self.y /= factor;
        self.z /= factor;
    }

    #[inline(always)]
    pub fn norm(&self) -> f64 {
        sqrt!(self.dot(self))
    }

    #[inline(always)]
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    #[inline(always)]
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    #[inline(always)]
    pub fn negate(&mut self) {
        self.mul(-1.)
    }

    #[inline(always)]
    pub fn zero(&mut self) {
        self.x = 0.0;
        self.y = 0.0;
        self.z = 0.0;
    }

    #[inline(always)]
    pub fn unpack(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero() {
        let mut a = Axes::from(1., 2., 3.);
        let zero = Axes::from(0., 0., 0.);
        a.zero();
        assert_eq!(&a, &zero);
    }

    #[test]
    fn test_netage() {
        let mut a = Axes::from(1., 2., 3.);
        let b = Axes::from(-1., -2., -3.);
        a.negate();
        assert_eq!(&a, &b);
    }

    #[test]
    fn test_add() {
        let mut a = Axes::from(1., 2., 3.);
        let b = Axes::from(-1., -2., -3.);
        let zero = Axes::from(0., 0., 0.);
        a.add(&b);
        assert_eq!(&a, &zero);
    }

    #[test]
    fn test_sub() {
        let mut a = Axes::from(1., 2., 3.);
        let b = Axes::from(1., 2., 3.);
        let zero = Axes::from(0., 0., 0.);
        a.sub(&b);
        assert_eq!(&a, &zero);
    }

    #[test]
    fn test_mul() {
        let mut a = Axes::from(1., 2., 3.);
        let b = Axes::from(-1., -2., -3.);
        a.mul(-1.);
        assert_eq!(&a, &b);
    }

    #[test]
    fn test_div() {
        let mut a = Axes::from(1., 2., 3.);
        let b = Axes::from(0.5, 1., 1.5);
        a.div(2.);
        assert_eq!(&a, &b);
    }

    #[test]
    fn test_norm() {
        let a = Axes::from(3., 4., 12.);
        let res = a.norm();
        let expected = 13.;
        assert_eq!(expected, res);
    }

    #[test]
    fn test_dot() {
        let a = Axes::from(1., 2., 3.);
        let b = Axes::from(0.5, 1., 1.5);
        let res = a.dot(&b);
        let expected = 7.;
        assert_eq!(expected, res);
    }

    #[test]
    fn test_cross() {
        let a = Axes::from(1., 2., 3.);
        let b = Axes::from(0.5, 1., 1.5);
        let res = a.cross(&b);
        let zero = Axes::from(0., 0., 0.);
        assert_eq!(zero, res);
    }
}
