/*
  This file contains parameters for HAWK 256, HAWK 512, and HAWK 1024.
  Parameter is retrieved by calling params_i or params_f for integer and
  rational parameter, respectively, with the size of log n (8,9,10).
  If parameter keyword is invalid, params_i will return 0.
  If parameter keyword is invalid, params_f will return NAN, and can be handled
  by checking if result.is_nan().
*/

pub fn params_i(logn: usize, param: &str) -> usize {
    match logn {
        8 => return hawk_256_params_i(param),
        // 9=> return hawk_512_params_i(param),
        // 10=> return hawk_1024_params_i(param),
        //default case
        _ => return 0,
    }
}

pub fn params_f(logn: usize, param: &str) -> f64 {
    match logn {
        8 => return hawk_256_params_f(param),
        // 9=> return hawk_512_params_f(param),
        // 10=> return hawk_1024_params_f(param),
        // default case
        _ => return f64::NAN,
    }
}

fn hawk_256_params_i(param: &str) -> usize {
    return match param {
        "lenpriv" => 96,
        "lenpub" => 450,
        "lensig" => 249,
        "n" => 256,
        "qs" => 1 << 32,
        "lensalt" => 112 / 8,
        "lenkgseed" => 128 / 8,
        "lenhpub" => 128 / 8,
        "low00" => 5,
        "high00" => 9,
        "low01" => 8,
        "high01" => 11,
        "high11" => 13,
        "highs0" => 12,
        "lows1" => 5,
        "highs1" => 9,
        _ => 0,
    };
}

fn hawk_256_params_f(param: &str) -> f64 {
    return match param {
        "sigmasign" => 1.010,
        "sigmaverify" => 1.042,
        "sigmakrsec" => 1.042,
        "beta0" => 1.0 / 250.0,
        _ => f64::NAN,
    };
}

fn hawk_512_params_f(param: &str) -> f64 {
    return 0.0;
}

fn hawk_1024_params_f(param: &str) -> f64 {
    return 0.0;
}
