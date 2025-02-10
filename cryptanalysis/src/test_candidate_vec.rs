use hawklib::ntru_solve::ntrusolve;
use hawklib::utils::rot_key;
use nalgebra::*;

use indicatif::{ProgressBar, ProgressStyle};

pub fn test_candidate_vec(candidate_vec: &DVector<i32>, pubkey: &DMatrix<i32>) -> bool {
    // given a test candidate from gradient descent, create all 2*n possible private keys B' by rot().
    //
    // For each created private key:
    // - check if B't B' = Q
    // - sign some message and check if signature is accepted by public key Q

    // there are two cases:
    // 1: vec is {G, -g}
    // 2: vec is {-F, f}

    let n = candidate_vec.len() / 2;

    let mut cvec = candidate_vec.clone();

    let pb = ProgressBar::new(2 * n as u64);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );

    // case 1
    for i in (0..2 * n) {
        // we need to use ntrusolve to find -F and f as fG - gF = 1
        // extract G and -g out from column of B inverse
        let (left, right) = cvec.as_slice().split_at(n);
        let (bigg, g): (Vec<i64>, Vec<i64>) = (
            left.to_vec().iter().map(|&x| x as i64).collect(),
            right.to_vec().iter().map(|&x| x as i64).collect(),
        );
        // find -F and f that satisfy ntru equation
        if let Some((f, bigf)) = ntrusolve(&bigg, &g) {
            println!("f: {:?}", f);
            println!("g: {:?}", g);
            println!("F: {:?}", bigf);
            println!("G: {:?}", bigg);
            // convert polynomials to 2n x 2n matrix
            let b = to_dmatrix(rot_key(&f, &g, &bigf, &bigg));
            println!("Det(B) = {}", b.map(|x| x as f64).determinant());

            // check if BtB = Q. If so we can forge signature
            if &b.transpose() * b == pubkey.map(|x| x as i64) {
                eprintln!("Found! \n{} is correct!", cvec);
                return true;
            }
        }
        break;
        pb.inc(1);

        // shift the vector
        // after 2n times this has gone through all rotations for +/- cvec
        shift(&mut cvec);
    }

    // reset our candidate vector again
    // cvec = candidate_vec.clone();
    // pb.reset();
    // // case 2
    // for i in (0..2*n) {
    //     // we need to use ntrusolve to find G and -g such that fG - gF = 1
    //
    //     // extract -F and f out from column of B inverse
    //     let (left, right) = cvec.as_slice().split_at(n);
    //     let (bigf, f): (Vec<i64>, Vec<i64>) = (
    //         left.to_vec().iter().map(|&x| -x as i64).collect(),
    //         (right).to_vec().iter().map(|&x| x as i64).collect()
    //         );
    //
    //     // find -F and f that satisfy ntru equation
    //     if let Some((bigg, g)) = ntrusolve(&f, &bigf){
    //
    //         // convert polynomials to 2n x 2n matrix
    //         let b = to_dmatrix(rot_key(&f, &g, &bigf, &bigg));
    //
    //         // check if BtB = Q. If so we can forge signature
    //         if &b.transpose() * b == pubkey.map(|x| x as i64) {
    //             eprintln!("Found! \n{} is correct!", cvec);
    //             return true;
    //         }
    //     }
    //
    //     // shift the vector
    //     // after 2n times this has gone through all rotations for +/- cvec
    //     pb.inc(1);
    //     shift(&mut cvec);
    // }

    false
}

fn to_dmatrix(b: Vec<Vec<i64>>) -> DMatrix<i64> {
    let n = b.len();
    let bflat: Vec<i64> = b.clone().into_iter().flatten().collect();
    DMatrix::from_column_slice(n, n, &bflat)
}
//
// fn to_dvector(f: &Vec<i64>) -> DVector<i64> {
//
// }

fn shift(f: &mut DVector<i32>) {
    let n = f.len() / 2;

    // Mutable slices for each half
    let slice1 = f.as_mut_slice();

    // Rotate first half
    let last1 = -slice1[n - 1]; // Last element of first half
    for i in (1..n).rev() {
        slice1[i] = slice1[i - 1]; // Shift right
    }
    slice1[0] = last1; // Move last to first

    // Rotate second half
    let last2 = -slice1[2 * n - 1]; // Last element of second half
    for i in ((n + 1)..(2 * n)).rev() {
        slice1[i] = slice1[i - 1]; // Shift right
    }
    slice1[n] = last2; // Move last to first of second half
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_shift() {
        let mut a = DVector::from_vec(vec![1, 2, 3, 4, 7, 8, 9, 10]);
        let b = DVector::from_vec(vec![-4, 1, 2, 3, -10, 7, 8, 9]);
        let c = DVector::from_vec(vec![-3, -4, 1, 2, -9, -10, 7, 8]);
        shift(&mut a);
        assert_eq!(b, a);
        shift(&mut a);
        assert_eq!(c, a);
    }
}
