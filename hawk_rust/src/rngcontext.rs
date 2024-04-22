use sha3::{
    digest::{ExtendableOutput, Update},
    Digest, Shake256,
};

use rand::{random, rngs::StdRng, Rng, SeedableRng};

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
        let seed_temp = shake256(&temp.to_ne_bytes());

        // println!("{:?}", seed_temp);

        let mut rng = StdRng::from_seed(seed_temp);
        let rand_u128: u128 = rng.gen();
        return rand_u128 >> (128-size);

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
