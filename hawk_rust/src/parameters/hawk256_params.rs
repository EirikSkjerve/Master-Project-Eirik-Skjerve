// HAWK 256 parameter set

pub const LENPRIV: usize = 96;
pub const LENPUB: usize = 450;
pub const LENSIG: usize = 249;

pub const N: usize = 256;
pub const QS: usize = 1<<32;

pub const SIGMASIGN: f64 = 1.010;
pub const SIGMAVERIFY: f64 = 1.042;
pub const SIGMAKRSEC: f64 = 1.042;

pub const LENSALT: usize = 14;
pub const LENKGSEED: usize = 16;
pub const LENHPUB: usize = 16;

pub const LOW00: usize = 5;
pub const HIGH00: usize = 9;
pub const LOW01: usize = 8;
pub const HIGH01: usize = 11;
pub const HIGH11: usize = 13;
pub const HIGHS0: usize = 12;
pub const LOWS1: usize = 5;
pub const HIGHS1: usize = 9;

pub const BETA0: f64 = 1.0 / 250.0;

pub const T0: [u128; 10] = [
    0x26B871FBD58485D45050,
    0x07C054114F1DC2FA7AC9,
    0x00A242F74ADDA0B5AE61,
    0x0005252E2152AB5D758B,
    0x00000FDE62196C1718FC,
    0x000000127325DDF8CEBA,
    0x0000000008100822C548,
    0x00000000000152A6E9AE,
    0x0000000000000014DA4A,
    0x0000000000000000007B,
];

pub const T1: [u128; 10] = [
    0x13459408A4B181C718B1,
    0x027D614569CC54722DC9,
    0x0020951C5CDCBAFF49A3,
    0x0000A3460C30AC398322,
    0x000001355A8330C44097,
    0x00000000DC8DE401FD12,
    0x00000000003B0FFB28F0,
    0x00000000000005EFCD99,
    0x00000000000000003953,
    0x00000000000000000000,
];

