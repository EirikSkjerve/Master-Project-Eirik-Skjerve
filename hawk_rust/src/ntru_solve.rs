






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

    if a_inp < b_inp{
        return xgcd(b_inp, a_inp);
        // TODO need to reorder the coeficcients
    }
    let mut cof: [i64; 4] = [1, 0, 0, 1];
    let mut a = a_inp;
    let mut b = b_inp;
    let mut q = 0;
    let mut temp = 0;

    while b != 0 {

        // rounded division
        // q = (a as f64 / b as f64).round() as i64;
        q = a/b;

        temp = b;
        b = modulo(a, b);
        a = temp;

        temp = cof[1];
        cof[1] = cof[0] - q*cof[1];
        cof[0] = temp;
        temp = cof[3];
        cof[3] = cof[2] - q*cof[3];
        cof[2] = temp;
    }

    return (a, cof[0], cof[2]);
}
