use hawklib::hawksign::hawksign_total_h;
use hawklib::hawkkeygen::hawkkeygen;

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use rand::Rng;

pub fn measure_uy_pipeline(n: usize, num_samples: usize) {

    let (privkey, _) = hawkkeygen(n);

    for i in 0..num_samples {
    }

}

fn get_random_bytes(num_bytes: usize) -> Vec<u8> {
    let mut res = Vec::with_capacity(num_bytes);
    let mut rng = rand::thread_rng();

    let mut mu = 0;
    for _ in 0..num_bytes {
        res.push(rng.gen_range(0..255));
    }

    res
}
