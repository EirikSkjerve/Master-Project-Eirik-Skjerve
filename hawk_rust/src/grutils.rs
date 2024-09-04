// utils file for the Golomb-Rice compression/decompression algorithms

use num::traits::{FromPrimitive, Num, PrimInt};

pub fn modulo<T: PrimInt>(a: T, b: T) -> T
where
    T: Num + FromPrimitive,
{
    // convert the inputs to u128
    let a_i32 = a.to_i32().unwrap();
    let b_i32 = b.to_i32().unwrap();

    // perform the calculations
    let result = ((a_i32 % b_i32) + b_i32) % b_i32;

    return T::from_i32(result).unwrap();
}

pub fn bin(a: i32, x: usize) -> Vec<u8> {
    /*
    Converts an integer to binary representation in a vector of arbitrary size
    Call bin(a, x=0) for default binary size
    */

    if x == 0 {
        // returns a as binary with *default* size
        // return bin(a, (a as f64).log2().ceil() as usize + 1);
        return vec![];
    }

    let mut res: Vec<u8> = vec![];
    let mut b = a;
    for i in 0..x {
        if b % 2 == 1 {
            res.push(1);
        }
        if b % 2 == 0 {
            res.push(0);
        }
        b /= 2;
    }

    // depending on bit-order, this should be reversed
    return res;
}

pub fn int(a: &Vec<u8>) -> i32 {
    let mut res: i32 = 0;
    for i in 0..a.len() {
        if a[i] == 1 {
            res += 1<<(a.len() - (i+1));
        }
    }
    return res;

}

pub fn packbits(bits: &Vec<u8>) -> Vec<u8> {
    let mut packed_bytes = Vec::new();

    for chunk in bits.chunks(8) {
        let mut byte: u8 = 0;

        for (i, &bit) in chunk.iter().enumerate() {

            if bit == 1 {
                byte |= 1 << i;
            }
        }
        packed_bytes.push(byte);
    }

    return packed_bytes;
}

pub fn unpackbits(arr: &Vec<u8>) -> Vec<u8> {
    let mut bits = Vec::new();

    for &byte in arr {
        for i in 0..8 {
            let bit = (byte>>i) & 1;
            bits.push(bit);
        }
    }
    return bits;
}
