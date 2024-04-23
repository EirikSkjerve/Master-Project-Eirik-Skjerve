use sha3::{
    digest::{ExtendableOutput, Update},
    Shake256,
};

use rand::{rngs::StdRng, Rng, SeedableRng};

pub struct RngContext {
    seed: u128,
    i:u8
}

impl RngContext {
    
    pub fn new(init_seed:u128) -> Self{
        RngContext {
            seed: init_seed,
            i: 0 
        }
    }

    /*
    Returns size random bits,
     */
    pub fn rnd(&mut self, size: u8) -> u128{
        self.i += 1;
        let temp = self.seed + self.i as u128;
        let seed = shake256(&temp.to_ne_bytes());  // set seed to hash of input seed + incremented variable

        let mut rng = StdRng::from_seed(seed);  // start rng from seed
        let rand_u128: u128 = rng.gen();  // generate random number from rng
        return rand_u128 >> (128-size);  // return number reduced to wanted bitsize

    }
}


pub fn shake256(input_data: &[u8]) -> [u8;32] {
    
    //TODO shake256 should have variable output length. Don't know how to best to this

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
