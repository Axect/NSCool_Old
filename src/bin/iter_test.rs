extern crate NSCool;
use NSCool::structure::eos::iter_candidate;

fn main() {
    let a = iter_candidate(6, 2);
    for x in a {
        println!("{:?}", x);
    }
}