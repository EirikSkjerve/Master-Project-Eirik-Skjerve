use crate::grutils::*;

pub fn decompressgr(y: &Vec<bool>, k: usize, low: usize, high: usize) -> (Vec<i64>, usize){
    assert!(k%8 == 0);

    if y.len() < k*(low+2) {
        // return some failure value here
        return (vec![0], 0);
    }

    let mut x: Vec<i64> = vec![0; k];

    // convert vector of bools to vector of u8's
    let mut y_u8: Vec<u8> = Vec::with_capacity(y.len());
    for i in 0..y.len() {
        if y[i]{
            y_u8.push(1);
        }
        if !y[i]{
            y_u8.push(0);
        }
    }

    for i in 0..k {
        x[i] = int(&y_u8[((i*low)+k)..((i+1)*low+k)].to_vec()) as i64;
    }

    let mut j = k*(low+1);

    for i in 0..k{

        let mut z: i64 = -1;
        let mut t = 0;

        while t != 1{
            z += 1;

            if j >= y.len() && z >= (i<<(high-low)) as i64 {
                // return some failure value here
                return (vec![0], 0);
            }
            t = y_u8[j];
            j += 1;
        }
        x[i] += z*(1<<low);
    }

    for i in 0..k{
        x[i] -= (y_u8[i] as i64) * ((2*x[i]) + 1);
    }

    return (x, j);

}
