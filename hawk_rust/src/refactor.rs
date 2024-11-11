use crate::{hawkkeygen, hawksign, hawkverify};

pub fn run_refactor() {
    // degree n
    let n = 256;
    // create public private keypair
    println!("Generating keypair");
    let (privkey, pubkey) = hawkkeygen::hawkkeygen(n);
    println!("Keypair generated");

    // create some message
    let message: [u8; 8] = [1,2,3,4,5,6,7,8];

    println!("Generating signature for message {:?}", message);
    let (signature, salt) = hawksign::hawksign(&privkey, &message, n);
    println!("Signature generated");

    let verification = hawkverify::hawkverify(&message, &pubkey, &signature, &salt, n);

    println!("Verification: {}", verification);

}
