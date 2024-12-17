mod collect_signatures;
mod dgd_estimation;
mod hpp;

use collect_signatures::generate_t_signatures;
use dgd_estimation::estimate_sigma;
use hpp::hpp::run_hpp_attack;

fn main() {
    // run_hpp_attack();
    // generate_t_signatures(100000,256);
    estimate_sigma(10000, 256);
}
