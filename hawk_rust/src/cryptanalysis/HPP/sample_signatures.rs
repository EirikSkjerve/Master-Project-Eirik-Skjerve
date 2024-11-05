use crate::hawk256::{hawkkeygen_256::hawkkeygen_256, hawksign_256::hawksign_256};
use crate::rngcontext::get_random_bytes;
use crate::write_to_file::wvtf;

fn sample_sigs(n: usize) -> (Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<u8>) {
    let init_seed = get_random_bytes(15);

    // generate keypair
    let (privkey, pubkey) = hawkkeygen_256(&init_seed);

    // generate n random messages
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(n);
    for _ in 0..n {
        messages.push(get_random_bytes(100));
    }

    // vector to keep all signatures
    let mut sigs: Vec<Vec<u8>> = Vec::with_capacity(n);

    // create signatures for each of the messages
    for i in 0..n {
        let sig = hawksign_256(&privkey, &messages[i]);
        sigs.push(sig);
    }

    // return the signatures, messages, and the public key
    (sigs, messages, pubkey)
}

pub fn write_samples_to_file(n: usize, path: &str) {
    let (sigs, messages, pubkey) = sample_sigs(n);

    // write public key at top of file
    wvtf(path, &pubkey, "\npk");

    for i in 0..sigs.len() {
        wvtf(path, &messages[i], "\nmsg");
        wvtf(path, &sigs[i], "sig");
    }
}
