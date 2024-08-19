






/*

pub fn ntrusolve(f: &Vec<i64>, g: &Vec<i64>) -> (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>){


    // return (f, g, f, g);
}

pub fn reduce(f: &Vec<i64>, g: &Vec<i64>, G: &Vec<i64>, G: &Vec<i64>) -> (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>){

    // return (f, g, F, G);
}
*/
use crate::utils::modulo;
pub fn xgcd(a_inp: i64, b_inp: i64) -> (i64, i64, i64){
    /*
     * Implements extended euclidean algorithm to find Bezout's coefficients
     * Inputs: integers a and b
     * Outputs: gcd(a,b) and s, t such that a*s + b*t = gcd(a,b)
    */

    // swap the order if a is less than b
    if a_inp < b_inp{
        let res = xgcd(b_inp, a_inp);
        return (res.0, res.2, res.1);
    }

    let mut cof: [i64; 4] = [1, 0, 0, 1];
    let mut a = a_inp;
    let mut b = b_inp;
    let mut q = 0;
    let mut temp = 0;

    while b != 0 {

        // rounded division
        q = a/b;

        // calculates the gcd
        temp = b;
        b = modulo(a, b);
        a = temp;

        // calculates the coefficients
        temp = cof[1];
        cof[1] = cof[0] - q*cof[1];
        cof[0] = temp;

        temp = cof[3];
        cof[3] = cof[2] - q*cof[3];
        cof[2] = temp;
    }

    return (a, cof[0], cof[2]);
}
