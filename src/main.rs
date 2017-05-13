extern crate num;
extern crate core;
extern crate primal;
extern crate rand;
extern crate num_cpus;
extern crate ansi_term;

use num::bigint::{BigUint, RandBigInt, ToBigUint};
use num::traits::{One, Zero};
use num::Integer;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime};
use std::thread;
use ansi_term::Colour::{Red, Green, Blue};

fn main() {
    let n0 = &BigUint::parse_bytes(b"833810193564967701912362955539789451139872863794534923259743419423089229206473091408403560311191545764221310666338878019", 10).unwrap(); // time: [impossibru]s [120],
    let n1 = &BigUint::parse_bytes(b"100433627766186892221372630609062766858404681029709092356097",10).unwrap(); // time: [impossibru]s lenght[61], r: 618970019642690137449562111 * 162259276829213363391578010288127
    let n2 = &BigUint::parse_bytes(b"72191103626161875816648846887", 10).unwrap(); // time: []s [29], r: 1500450271 *  48112959837082048697
    let n3 = &BigUint::parse_bytes(b"4292017463532640823", 10).unwrap(); // time: [97]s [19], r: 1500450271 *  2860486313
    let n4 = &BigUint::parse_bytes(b"608364911153957", 10).unwrap(); // time: [1]s [15], r: 4093082899 * 102841
    let n5 = &BigUint::parse_bytes(b"999962000357", 10).unwrap(); // time: 0s [13], r: 999979 * 999983
    let n6 = &1429229.to_biguint().unwrap(); // time: 0s [7], r: 2477 * 577

    run_factorization(n6);
    run_factorization(n5);
    run_factorization(n4);
    run_factorization(n3);
    run_factorization(n2);
    run_factorization(n1);
    run_factorization(n0);
}

fn run_factorization(n: &BigUint) {
    let now = SystemTime::now();

    println!("\nN: {}\n", Red.bold().paint(n.to_str_radix(10)));
    println!("{:?}", factor(&n));

    match now.elapsed() {
        Ok(elapsed) => {
            println!("Time elapsed {}s", Red.paint(elapsed.as_secs().to_string()));
        }
        Err(e) => {
            println!("Error: {:?}", e);
        }
    };

}

// Find an approximation of a square root (seems to return the good result plus one)
fn approx_sqrt(number: &BigUint, iterations: usize) -> BigUint {
    let mut approx = number.clone();

    for _ in 0..iterations {
        approx = (&approx + (number / &approx)) / 2.to_biguint().unwrap();
    }

    approx
}

/* Miller Rabin Functions 
 * 
 * Probablistic approch: is_prime?
 * If p is an odd prime number, and p - 1 = 2^s d, with d odd,
 * then for every a prime to p, either a*d = 1 mod p, or
 * there exists t such that 0 ≤ t < s and a2td = 1 mod p
 */

// find r and d s.t. n - 1 = 2^r * d
fn find_r_and_d(i: &BigUint) -> (u64, BigUint) {
    let mut d = i.clone();
    let mut r = 0;
    let two = BigUint::new(vec![2]);

    while d.is_even() {
        d = d / two.clone();
        r += 1;
    }
    (r, d)
}

fn miller_rabin(candidate: &BigUint, limit: u64) -> bool {
    // precompute fixed numbers
    let one: BigUint = BigUint::one();
    let mut two: BigUint = BigUint::one() + BigUint::one();
    let mone: BigUint = candidate - one.clone();

    if *candidate == BigUint::one() {
        return false;
    }

    // find r and d s.t. n - 1 = 2^r * d
    let (r, mut d) = find_r_and_d(&mone.clone());

    // start random generator
    let mut rng = rand::thread_rng();
    let mut phi: BigUint;
    let mut app: BigUint;
    let mut index = 0;
    'wit: while index <= limit {
        index = index + 1;
        app = rng.gen_biguint_range(&two, &mone); // [start, end)
        phi = mod_exp(&mut app, &mut d, &candidate);
        if phi == one || phi == mone {
            continue 'wit;
        }
        for _ in 0..r {
            phi = mod_exp(&mut phi, &mut two, &candidate);
            if phi == one {
                return false;
            } else if phi == mone {
                continue 'wit;
            }
        }
        return false;
    }
    return true;
}
/* End Miller Rabin */


/* Little Fermat Function
 *
 * Probablistic approch: is_prime?
 * If p is a prime number, then for any integer a, the number a * p − a is an integer multiple of p.
 * This is expressed as a^p = a (mod p)
 */
fn little_fermat(num: &BigUint) -> bool {
    let mut rng = rand::thread_rng();
    // Pick up a random a
    let mut a: BigUint = rng.gen_biguint_below(num);
    let result = mod_exp(&mut a, &mut (num - BigUint::one()), num);

    // If the result is 1 the equality holds for a
    result == BigUint::one()
}

/* End Little Fermat Function */

// Compute factorization
fn factor(num: &BigUint) -> Vec<(BigUint, BigUint)> {
    // Output vector with lock system on write/read
    let factors: Arc<Mutex<Vec<(BigUint, BigUint)>>> = Arc::new(Mutex::new(Vec::new()));
    let cpus = num_cpus::get();
    // Split the solution set
    let step: BigUint = approx_sqrt(num, 1000) / (cpus.to_biguint().unwrap());
    //let step: BigUint = (num / 2.to_biguint().unwrap()) / (cpus.to_biguint().unwrap());
    // TODO: if p or q are lower than the square fail
    let mut start = 3.to_biguint().unwrap();
    let mut end = start.clone() + step.clone();
    let mut threads = Vec::new();

    for _ in 0..cpus {
        // needs to copy values, for lifetime
        let st = start.clone();
        let en = end.clone();
        let n = num.clone();
        let lock: Arc<Mutex<Vec<(BigUint, BigUint)>>> = factors.clone();
        let zero = BigUint::zero();
        threads.push(thread::spawn(move || {
            // run only for odd numbers
            let mut index = if st.is_even() { st + BigUint::one() } else { st };
            while index <= en {
                // n mod i and i is prime -> found p or q
                if n.clone() % index.clone() == zero && is_prime(&index) {
                    let q = n.clone() / index.clone();
                    let mut _lock = lock.lock().unwrap();
                    if !_lock.iter().any(|p| p.0 == index || p.1 == index) && is_prime(&q) {
                        // TODO: maybe block and exit?
                        println!("{} FOUND: {} and {}",
                                 Blue.bold().paint("[!]"),
                                 Green.bold().paint(index.to_str_radix(10)),
                                 Green.bold().paint(q.to_str_radix(10)));
                        _lock.push((index.clone(), q.clone()));
                    }

                }
                index = index + 2.to_biguint().unwrap();
            }
        }));
        start = end + BigUint::one();
        end = start.clone() + step.clone();
    }

    for t in threads {
        handle(t.join());
    }

    // Needs to first check if the lock is free, then unwrap the guard
    let guard = Arc::try_unwrap(factors).expect("Lock is owned");
    // Get into the real vector
    guard.into_inner().expect("Error on acquiring lock")
}

// Check if the number is a prime using Fermat and/or Miller Rabin
fn is_prime(num: &BigUint) -> bool {
    //little_fermat(num) && miller_rabin(num, 128);
    // TODO: little_fermat or miller_rabin have equal performances
    little_fermat(num)
}

/* Utils Functions */

// Compute base^exponent mod modulus
fn mod_exp(base: &mut BigUint, exponent: &mut BigUint, modulus: &BigUint) -> BigUint {
    let one = BigUint::one();
    let zero = BigUint::zero();
    let mut result = BigUint::one();

    while *exponent > zero {
        if &*exponent & one.clone() == one {
            result = (result * &*base) % &*modulus;
        }
        *base = (&*base * &*base) % &*modulus;
        *exponent = &*exponent >> 1usize;
    }
    result
}

// Handle the threads exit and joins
fn handle(r: thread::Result<()>) {
    match r {
        Ok(r) => r,
        Err(e) => {
            if let Some(e) = e.downcast_ref::<&'static str>() {
                println!("Got an error: {}", e);
            } else {
                println!("Got an unknown error: {:?}", e);
            }
        }
    }
}

/* End Utils Functions */
