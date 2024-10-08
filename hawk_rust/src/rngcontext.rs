extern crate sha3;

use rand::Rng;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

// RngContext struct equivalent
pub struct RngContext {
    shake256: Shake256,
    i: usize,
}

impl RngContext {
    // Constructor for RngContext
    pub fn new(seed: &[u8]) -> Self {
        let mut shake256 = Shake256::default();
        shake256.update(seed);
        RngContext { shake256, i: 0 }
    }

    // Random function to generate deterministically randomness
    pub fn random(&mut self, size: usize) -> Vec<u8> {
        let mut reader = self.shake256.clone().finalize_xof();
        let mut result: Vec<u8> = vec![0; self.i + size];
        reader.read(&mut result);
        let r = result[self.i..self.i + size].to_vec();
        self.i += size;
        return r;
    }
}

// SHAKE256x4 equivalent in Rust
pub fn shake256x4(message: &[u8], num: usize) -> Vec<u64> {
    let mut digest = Vec::new();

    // Initialize 4 SHAKE256 instances
    for i in 0..4 {
        let mut shake = Shake256::default();
        let mut input = message.to_vec();
        input.push(i as u8);
        shake.update(&input);
        let mut reader = shake.finalize_xof();
        let mut out = vec![0u8; (num * 8) / 4];
        reader.read(&mut out);
        digest.push(out);
    }

    // Convert the four streams into interleaved form
    let mut y = vec![0u64; num];
    let mut j = 0;

    for i in 0..(num / 4) {
        y[j] = u64::from_le_bytes(digest[0][i * 8..(i + 1) * 8].try_into().unwrap());
        y[j + 1] = u64::from_le_bytes(digest[1][i * 8..(i + 1) * 8].try_into().unwrap());
        y[j + 2] = u64::from_le_bytes(digest[2][i * 8..(i + 1) * 8].try_into().unwrap());
        y[j + 3] = u64::from_le_bytes(digest[3][i * 8..(i + 1) * 8].try_into().unwrap());
        j += 4;
    }

    return y;
}

pub fn get_random_bytes(num_bytes: usize) -> Vec<u8> {
    let mut res = Vec::with_capacity(num_bytes);
    let mut rng = rand::thread_rng();

    for _ in 0..num_bytes {
        res.push(rng.gen_range(0..255));
    }

    return res;
}
