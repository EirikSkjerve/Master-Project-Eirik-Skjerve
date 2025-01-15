use bincode;
use std::fs::File;
use std::io::{Read, Write};

pub fn write_vectors_to_file(
    vs: Vec<Vec<i16>>,
    sk: (Vec<u8>, Vec<i64>, Vec<i64>),
    path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let encoded: Vec<u8> = bincode::serialize(&(vs, sk))?;
    let mut file = File::create(format!("{path}.bin"))?;
    file.write_all(&encoded)?;

    Ok(())
}

pub fn read_vectors_from_file(
    path: &str,
) -> Result<(Vec<Vec<i16>>, (Vec<u8>, Vec<i64>, Vec<i64>)), Box<dyn std::error::Error>> {
    let mut file = File::open(format!("{path}.bin"))?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    let vs: (Vec<Vec<i16>>, (Vec<u8>, Vec<i64>, Vec<i64>)) = bincode::deserialize(&buffer)?;

    return Ok((vs));
}
