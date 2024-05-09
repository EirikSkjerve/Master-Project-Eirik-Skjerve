use crate::rngcontext::{shake256x4, RngContext};



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
    generate_f_g(1337, logn)

}

fn generate_f_g(seed: usize, logn: u16) {
    let n = 1 << logn;
    let b = n/64;
    assert!(b==4 || b==8 || b==16);

    let y = shake256x4(&seed.to_ne_bytes(), 2*n*b/64);
    println!("Vector y: {:?}", y);
}