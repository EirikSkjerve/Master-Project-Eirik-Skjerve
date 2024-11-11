use crate::{hawkkeygen, hawksign, hawkverify};

pub fn run_refactor() {
    hawkkeygen::hawkkeygen_256();
    hawkkeygen::hawkkeygen_512();
    hawkkeygen::hawkkeygen_1024();
}
