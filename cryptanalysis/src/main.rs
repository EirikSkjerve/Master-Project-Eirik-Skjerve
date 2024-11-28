mod collect_signatures;
mod hpp;

use hpp::hpp::run_hpp_attack;
use collect_signatures::generate_t_signatures;

fn main() {
    // run_hpp_attack();
    generate_t_signatures(900000,256);
}
