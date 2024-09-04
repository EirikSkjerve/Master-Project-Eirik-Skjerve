use crate::grutils::*;

pub fn decompressgr(y: &Vec<u8>, k: usize, low: usize, high: usize) -> (Vec<i64>, usize){
    assert!(k%8 == 0);

    if y.len() < k*(low+2) {
        // return some failure value here
        return (vec![0], 0);
    }

    let mut x: Vec<i64> = vec![0; k];


    for i in 0..k {
        x[i] = int(&y[((i*low)+k)..((i+1)*low+k)].to_vec()) as i64;
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
            t = y[j];
            j += 1;
        }
        x[i] += z*(1<<low);
    }

    for i in 0..k{
        x[i] -= (y[i] as i64) * ((2*x[i]) + 1);
    }

    return (x, j);

}
