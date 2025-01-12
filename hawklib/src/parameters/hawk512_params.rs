// HAWK 512 parameter set

// many are commented out because they are not used when encoding/decoding is not used

// pub const LENPRIV: usize = 184;
// pub const LENPUB: usize = 1024;
// pub const LENSIG: usize = 555;
//
// pub const N: usize = 512;
// pub const QS: u128 = 1 << 64;
//
pub const SIGMASIGN: f64 = 1.278;
pub const SIGMAVERIFY: f64 = 1.425;
pub const SIGMAKRSEC: f64 = 1.425;

pub const LENSALT: usize = 24;
pub const LENKGSEED: usize = 24;
// pub const LENHPUB: usize = 32;

// pub const LOW00: usize = 5;
pub const HIGH00: usize = 9;
// pub const LOW01: usize = 9;
pub const HIGH01: usize = 12;
pub const HIGH11: usize = 15;
pub const HIGHS0: usize = 13;
// pub const LOWS1: usize = 5;
pub const HIGHS1: usize = 9;

pub const BETA0: f64 = 1.0 / 1000.0;

pub const T0: [u128; 13] = [
    0x2C058C27920A04F8F267,
    0x0E9A1C4FF17C204AA058,
    0x02DBDE63263BE0098FFD,
    0x005156AEDFB0876A3BD8,
    0x0005061E21D588CC61CC,
    0x00002BA568D92EEC18E7,
    0x000000CF0F8687D3B009,
    0x0000000216A0C344EB45,
    0x0000000002EDF0B98A84,
    0x0000000000023AF3B2E7,
    0x00000000000000EBCC6A,
    0x000000000000000034CF,
    0x00000000000000000006,
];

pub const T1: [u128; 13] = [
    0x1AFCBC689D9213449DC9,
    0x06EBFB908C81FCE3524F,
    0x01064EBEFD8FF4F07378,
    0x0015C628BC6B23887196,
    0x0000FF769211F07B326F,
    0x00000668F461693DFF8F,
    0x0000001670DB65964485,
    0x000000002AB6E11C2552,
    0x00000000002C253C7E81,
    0x00000000000018C14ABF,
    0x0000000000000007876E,
    0x0000000000000000013D,
    0x00000000000000000000,
];
