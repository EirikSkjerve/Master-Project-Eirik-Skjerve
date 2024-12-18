// HAWK 1024 parameter set

// many are commented out because they are not used when encoding/decoding is not used
//
// pub const LENPRIV: usize = 360;
// pub const LENPUB: usize = 2440;
// pub const LENSIG: usize = 1221;
//
// pub const N: usize = 1024;
// pub const QS: u128 = 1 << 64;
//
pub const SIGMASIGN: f64 = 1.299;
pub const SIGMAVERIFY: f64 = 1.571;
pub const SIGMAKRSEC: f64 = 1.974;

pub const LENSALT: usize = 40;
pub const LENKGSEED: usize = 40;
// pub const LENHPUB: usize = 64;

// pub const LOW00: usize = 6;
pub const HIGH00: usize = 10;
// pub const LOW01: usize = 10;
pub const HIGH01: usize = 14;
pub const HIGH11: usize = 17;
pub const HIGHS0: usize = 14;
// pub const LOWS1: usize = 6;
pub const HIGHS1: usize = 10;

pub const BETA0: f64 = 1.0 / 3000.0;

pub const T0: [u128; 13] = [
    0x2C583AAA2EB76504E560,
    0x0F1D70E1C03E49BB683E,
    0x031955CDA662EF2D1C48,
    0x005E31E874B355421BB7,
    0x000657C0676C029895A7,
    0x00003D4D67696E51F820,
    0x0000014A1A8A93F20738,
    0x00000003DAF47E8DFB21,
    0x0000000006634617B3FF,
    0x000000000005DBEFB646,
    0x00000000000002F93038,
    0x0000000000000000D5A7,
    0x00000000000000000021,
];

pub const T1: [u128; 13] = [
    0x1B7F01AE2B17728DF2DE,
    0x07506A00B82C69624C93,
    0x01252685DB30348656A4,
    0x001A430192770E205503,
    0x00015353BD4091AA96DB,
    0x000009915A53D8667BEE,
    0x00000026670030160D5F,
    0x00000000557CD1C5F797,
    0x00000000006965E15B13,
    0x00000000000047E9AB38,
    0x000000000000001B2445,
    0x000000000000000005AA,
    0x00000000000000000000,
];
