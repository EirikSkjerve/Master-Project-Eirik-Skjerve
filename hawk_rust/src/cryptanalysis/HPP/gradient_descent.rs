use nalgebra as na;
use na::*;
use crate::rngcontext::{RngContext, get_random_bytes};
use std::time::{Instant};

use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand::rngs::StdRng;

pub fn gradient_descent(u:Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>) {

}
