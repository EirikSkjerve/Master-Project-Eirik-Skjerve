use sha3::{
    digest::{ExtendableOutput, Update},
    Shake256,
};

use rand::{rngs::StdRng, Rng, SeedableRng};

pub struct RngContext {
    seed: u128,
    i: usize,
}

impl RngContext {
    pub fn new(init_seed: u128) -> Self {
        RngContext {
            seed: init_seed,
            i: 0,
        }
    }

    /*
    Returns size random bits,
     */
    pub fn rnd(&mut self, size: u8) -> u128 {
        self.i += size as usize;
        // println!("Recursion step {}", self.i);

        let temp = self.seed + self.i as u128;
        // set seed to hash of input seed + incremented variable
        let seed = shake256(&temp.to_ne_bytes());

        let mut rng = StdRng::from_seed(seed); // start rng from seed
        let rand_u128: u128 = rng.gen(); // generate random number from rng
        return rand_u128 >> (128 - size); // return number reduced to wanted bitsize
    }
}

pub fn shake256(input_data: &[u8]) -> [u8; 32] {
    //TODO shake256 should have variable output length. Should use a Vec instead of array

    // create a shake256 object
    let mut hasher = Shake256::default();

    // update the hasher with bytes
    hasher.update(input_data);

    // create byte array to store the digest in
    let mut res: [u8; 32] = [0; 32];
    //
    hasher.finalize_xof_into(&mut res);

    return res;
}

pub fn shake256x4(message: &[u8], num: usize) -> Vec<u64> {
    /*
    Inputs: message m
    Outputs: four shake256 digests interleaved, with num length
    */

    // This should ideally be inside a loop, but because of Rust's borrow system, it is
    // not possible to store the Shake256 instances inside an array.

    // four shake256 instances
    let mut s_1 = Shake256::default();
    let mut s_2 = Shake256::default();
    let mut s_3 = Shake256::default();
    let mut s_4 = Shake256::default();

    // initialize an array to store the digests
    let mut digest: [Vec<u8>; 4] = [vec![0; num], vec![0; num], vec![0; num], vec![0; num]];

    // insert message
    s_1.update(message);
    s_2.update(message);
    s_3.update(message);
    s_4.update(message);

    // insert extra data to make the digests different
    s_1.update(&[1 as u8]);
    s_2.update(&[2 as u8]);
    s_3.update(&[3 as u8]);
    s_4.update(&[4 as u8]);

    // store the digests
    digest[0] = s_1.finalize_boxed((num * 8) / 4).to_vec();
    digest[1] = s_2.finalize_boxed((num * 8) / 4).to_vec();
    digest[2] = s_3.finalize_boxed((num * 8) / 4).to_vec();
    digest[3] = s_4.finalize_boxed((num * 8) / 4).to_vec();

    // initialized a vector with the desired capacity
    let mut y: Vec<u64> = Vec::with_capacity(num);

    // interleaves the four digests
    for i in 0..(num / 4) {
        for (index, dig) in digest.iter().enumerate() {
            let start = i * 8;
            let end = start + 8;
            let bytes = &dig[start..end];
            y.push(bytes_to_u64(bytes));
        }
    }
    return y;
}

// simple function for converting bytes to an unsigned 64-bit integer
fn bytes_to_u64(bytes: &[u8]) -> u64 {
    let mut result: u64 = 0;
    for &byte in bytes.iter() {
        result = result << 8 | byte as u64;
    }
    result
}
