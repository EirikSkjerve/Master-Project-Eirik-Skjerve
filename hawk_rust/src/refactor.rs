use crate::{hawkkeygen, hawksign, hawkverify};

pub fn run_refactor() {
    // degree n
    let ns: [usize; 3] = [256, 512, 1024];
    let num_rounds = 100;
    let mut num_fails = 0;
    for n in ns.iter(){

        println!("Creating keypairs for degree {}",n);
        let (privkey, pubkey) = hawkkeygen::hawkkeygen(*n); 

        println!("Generating {num_rounds} signatures");
        for i in 0..num_rounds {
            let message: [u8; 5] = [i, i, i, i, i];
            let (signature, salt) = hawksign::hawksign(&privkey, &message, *n);

            let verification = hawkverify::hawkverify(&message, &pubkey, &signature, &salt, *n);

            if !verification {
                num_fails += 1;
            }
        }
    }
    println!("Number of failures: {num_fails}");
}
