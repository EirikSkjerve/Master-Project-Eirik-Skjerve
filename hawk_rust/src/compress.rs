use crate::grutils;

// using Golomb-Rice compression algorithm
pub fn compressgr(x: &Vec<i64>, low: usize, high: usize) -> Vec<u8> {
    /*
    Input: sequence x, low and high
    Output: compressed sequence of bits
    If failure, [0] is returned
    */

    let k = x.len();
    assert_eq!(k % 8, 0);

    for i in 0..k {
        if !(x[i] < (1 << high) && x[i] > -(1 << high)) {
            println!("failure from compressgr 1");
            return vec![0];
        }
        // assert!(x[i] < (1<<high) && x[i] > -(1<<high));
    }

    let mut y: Vec<u8> = Vec::new();
    let mut v: Vec<i16> = Vec::new();

    for i in 0..k {
        let s = x[i] < 0;
        let si = match s {
            false => 0,
            true => 1,
        };

        y.push(si);
        v.push((x[i] - (si as i64) * (2 * x[i] + 1)) as i16);
        if v[i] >= (1 << high) {
            // empty vec indicates failed attempt
            println!("failure from compressgr 2");
            return vec![0];
        }
    }

    for i in 0..k {
        let res = grutils::bin(grutils::modulo(v[i] as i32, (1 << low) as i32), low);
        for r in res {
            y.push(r);
        }
    }

    for i in 0..k {
        let res = grutils::bin(0, (v[i] >> low) as usize);
        for r in res {
            y.push(r);
        }
        y.push(1);
    }

    return y;
}
