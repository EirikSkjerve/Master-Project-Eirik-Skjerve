// use std::fmt::Debug;
// use std::fs::OpenOptions;
// use std::io::prelude::*;
//
// // write vector to file
// pub fn wvtf<T: Debug>(path: &str, v: &Vec<T>, prefix: &str) -> std::io::Result<()> {
//     // let mut file = File::create(path)?;
//     let mut file = OpenOptions::new().append(true).create(true).open(path)?;
//     file.write_all(format!("{}: {:?} \n", prefix, v).as_bytes())?;
//     Ok(())
// }
