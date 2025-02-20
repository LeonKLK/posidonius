use serde::{Deserialize, Serialize};

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
    pub fn copy_from(&mut self, other: &Self) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }

    #[inline(always)]
    pub fn add(&mut self, other: &Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }

    #[inline(always)]
    pub fn mul(&mut self, factor: f64) {
        self.x *= factor;
        self.y *= factor;
        self.z *= factor;
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
