use crate::rngcontext::RngContext;



pub fn hawkkeygen(logn: u16, rng: Option<RngContext>) {
    
    // checks if rng-context is initialized or not. If not, initialize a new one and recusively call hawkkeygen
    let mut rng = match rng {
        Some(rng) => rng,
        None => {
            let new_rng = RngContext::new(1337); // this should be an actual random number
            return hawkkeygen(logn, Some(new_rng));
        },
    };

    // generate f and g

}

fn generate_f_g(seed: usize, logn: u16) {

}