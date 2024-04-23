use std::collections::HashMap;

// Parameters used in Hawk

const Hawk256_T0: [u128; 13] = [
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

    /// these are redundant
    0x00000000000000000000,
    0x00000000000000000000,
    0x00000000000000000000,
];

const Hawk512_T0: [u128; 13] = [
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

const Hawk1024_T0: [u128; 13] = [
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

const Hawk256_T1: [u128; 13] = [
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
    
    // these are redundant
    0x00000000000000000000,
    0x00000000000000000000,
    0x00000000000000000000,
];

const Hawk512_T1: [u128; 13] = [
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

const Hawk1024_T1: [u128; 13] = [
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

#[derive(Debug)]
pub enum HawkParams<'a>{
    Lenpriv(u64),
    Lenpub(u64),
    Lensig(u64),
    n(u64),
    qs(u128),
    SigmaSign(f32),
    SigmaVerify(f32),
    SigmaKrSec(f32),
    LenSalt(u64),
    LenkgSeed(u64),
    LenhPub(u64),
    Low00(u64),
    High00(u64),
    Low01(u64),
    High01(u64),
    High11(u64),
    Highs0(u64),
    Lows1(u64),
    Highs1(u64),
    Beta0(f32),
    T0(&'a [u128; 13]),
    T1(&'a [u128; 13]),

}

pub fn initialize_params<'a> (logn:u8) -> HashMap<&'a str, HawkParams<'a>>{
    if logn != 8 && logn != 9 && logn != 10 {
        panic!("Invalid param");
    }

    // initialize the parameters hashmap
    let mut params = HashMap::new();

    if logn == 8 {
        
        // update the parameter set with correct values for n=8
        params.insert("lenpriv", HawkParams::Lenpriv(96));
        params.insert("lenpub", HawkParams::Lenpub(450));
        params.insert("lensig", HawkParams::Lensig(249));
        params.insert("n", HawkParams::n(256));
        params.insert("qs", HawkParams::qs(2_u128.pow(32)));
        params.insert("sigmasign", HawkParams::SigmaSign(1.010));
        params.insert("sigmaverify", HawkParams::SigmaVerify(1.042));
        params.insert("sigmakrsec", HawkParams::SigmaKrSec(1.042));
        params.insert("lensalt", HawkParams::LenSalt(112/8));
        params.insert("lenkgseed", HawkParams::LenkgSeed(128/8));
        params.insert("lenhpub", HawkParams::LenhPub(128/8));
        params.insert("low00", HawkParams::Low00(5));
        params.insert("high00", HawkParams::High00(9));
        params.insert("low01", HawkParams::Low01(8));
        params.insert("high01", HawkParams::High01(11));
        params.insert("high11", HawkParams::High11(13));
        params.insert("highs0", HawkParams::Highs0(12));
        params.insert("lows1", HawkParams::Lows1(8));
        params.insert("highs1", HawkParams::Highs1(9));
        params.insert("beta0", HawkParams::Beta0(1_f32/250_f32));
        params.insert("T0", HawkParams::T0(&Hawk256_T0));
        params.insert("T1", HawkParams::T1(&Hawk256_T1));

    }
    if logn == 9 {

    }
    if logn == 10 {

    }

    return params;


}

/*

HAWK_256_PARAMS = {
    "lenpriv": 96,
    "lenpub": 450,
    "lensig": 249,
    "n": 256,
    "qs": 2**32,
    "sigmasign": 1.010,
    "sigmaverify": 1.042,
    "sigmakrsec": 1.042,
    "lensalt": 112 // 8,
    "lenkgseed": 128 // 8,
    "lenhpub": 128 // 8,
    "low00": 5,
    "high00": 9,
    "low01": 8,
    "high01": 11,
    "high11": 13,
    "highs0": 12,
    "lows1": 5,
    "highs1": 9,
    "beta0": 1 / 250,
    "T0": Hawk256_T0,
    "T1": Hawk256_T1,
}
*/
