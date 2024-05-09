
pub fn bin(a: u128, x: usize) -> Vec<u8> {
    /* 
    Converts an integer to binary representation in a vector of arbitrary size
    */
    let mut res: Vec<u8> = vec![];
    let mut b = a;

    for i in 0..x{
        if b % 2 == 1 { res.push(1); }
        if b % 2 == 0 { res.push(0); }
        b /= 2;
    }

    res.reverse();
    return res;
}

pub fn int<T: AsRef<[u8]>>(input: T) -> u128 {

    let a = input.as_ref();
    /* 
    Converts an array representing a binary string to its decimal representation
    and returns an 128 bit integer
    */
    let mut res: u128 = 0;
    for i in 0..a.len() {
        if a[i] == 1 {
            res += (2_u128).pow((a.len() - (i + 1)) as u32);
        }
    }
    return res;
}