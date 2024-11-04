use std::fmt::Debug;
use std::fs::File;
use std::io::prelude::*;

// write vector to file
pub fn wvtf<T: Debug>(path: &str, v: &Vec<T>, prefix: &str) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    println!("{path}");
    println!("{:?}", v);
    file.write_all(format!("{}, {:?}", prefix, v).as_bytes())?;
    Ok(())
}
