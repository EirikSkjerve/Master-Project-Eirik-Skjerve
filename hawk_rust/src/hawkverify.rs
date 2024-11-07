use sha3::{
    digest::{ExtendableOutputReset, Update},
    Shake256,
};

use crate::utils::{bytes_to_poly, modulo, poly_sub};
use crate::hawksign::symbreak;
use crate::verifyutils::{rebuildw0, poly_qnorm};
use crate::parameters::{hawk256_params, hawk512_params, hawk1024_params};

// poly times const
pub fn ptc(a: &Vec<i64>, b: i64) -> Vec<i64> {
    a.clone().iter().map(|&x| x*b).collect()
}

fn hawkverify_inner(
    msg: &[u8], 
    pub_key: (&Vec<i64>, &Vec<i64>), 
    signature: &Vec<u8>, 
    salt: &Vec<u8>,
    n: usize
    ) -> bool {

    // get the correct parameters
    let (highs0, highs1, high00, high01) = match n {
        256 => ( hawk256_params::HIGHS0,
                 hawk256_params::HIGHS1,
                 hawk256_params::HIGH00,
                 hawk256_params::HIGH01
            ),
        
        512 => ( hawk512_params::HIGHS0,
                 hawk512_params::HIGHS1,
                 hawk512_params::HIGH00,
                 hawk512_params::HIGH01
            ),

        _ => ( hawk1024_params::HIGHS0,
                 hawk1024_params::HIGHS1,
                 hawk1024_params::HIGH00,
                 hawk1024_params::HIGH01
            ),
    };

    // convert signature to Vec<i64>
    let s1: Vec<i64> = signature.iter().map(|&x| x as i64).collect();
    let (q00, q01) = pub_key;
    
    // compute hash digest of message m
    let mut shaker = Shake256::default();
    shaker.update(msg);
    let mut m: Vec<u8> = vec![0; 64];
    shaker.finalize_xof_reset_into(&mut m);

    // create buffer for vector h whose length depends on degree n
    let mut h: Vec<u8> = vec![0; n / 4];

    // digest h is digest of message+salt
    shaker.update(&m);
    shaker.update(&salt);

    // compute digest
    shaker.finalize_xof_reset_into(&mut h);

    // convert digest h to usable polynomials
    let (h0, h1) = (
        &bytes_to_poly(&h[0..n / 8], n),
        &bytes_to_poly(&h[n / 8..n / 4], n),
    );
    // reconstruct digest h

    let w1 = poly_sub(
        &h1,
        &ptc(&s1, 2)
        );

    if !symbreak(&w1) {
        println!("Symbreak failed");
        return false;
    }

    let w0 = rebuildw0(
        &q00, 
        &q01, 
        &w1, 
        &h0,
        highs0,
        highs1,
        high00,
        high01
        );


    false
}

pub fn hawkverify(
    msg: &[u8], 
    pub_key: (&Vec<i64>, &Vec<i64>), 
    signature: &Vec<u8>, 
    n: usize
) -> bool{

    assert!(n==256||n==512||n==1024);
    false
}
