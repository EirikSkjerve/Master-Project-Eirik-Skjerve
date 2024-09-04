use crate::grutils::*;

// using Golomb-Rice compression algorithm
pub fn compressgr(x: &Vec<i64>, low: usize, high: usize) -> Vec<bool> {

    /*
     Input: sequence x, low and high
     Output: compressed sequence of bits 
     If failure, [0] is returned
     */

    let k = x.len();
    assert_eq!(k%8, 0);

    for i in 0..k{
        assert!(x[i] < (1<<high) && x[i] > -(1<<high));
    }

    let mut y: Vec<bool> = Vec::new();
    let mut v: Vec<i16> = Vec::new();

    for i in 0..k {
        let s = x[i] < 0;
        let si = match s {
            false => 0,
            true => 1,
        };

        y.push(s);
        v.push((x[i]- si*(2*x[i] + 1)) as i16);
        if v[i] >= (1<<high){
            // empty vec indicates failed attempt
            return vec![false];
        }
    }

    for i in 0..k { 
        let res = bin(modulo(v[i] as i32, (1<<low) as i32), low);
        for r in res{
            if r==1{
                y.push(true);
            }
            if r==0{
                y.push(false)
            }
        }
    }

    for i in 0..k{
        for z in 0..(v[i]/(1<<low)) as usize {
            y.push(false);
        }
        y.push(true);
    }

    return y;
}

