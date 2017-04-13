extern crate num;
extern crate core;
extern crate primal;
extern crate rand;
extern crate num_cpus;
extern crate ansi_term;

use num::bigint::{BigUint, RandBigInt, ToBigUint};
use num::traits::{One, Zero, ToPrimitive};
use num::Integer;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, Instant};
use std::thread;
use ansi_term::Colour::{Red, Green, Blue};

fn main() {
    let now = SystemTime::now();

    let mut n = BigUint::parse_bytes(b"833810193564967701912362955539789451139872863794534923259743419423089229206473091408403560311191545764221310666338878019", 10).unwrap(); // time: impossibru

    //n = BigUint::parse_bytes(b"999962000357", 10).unwrap(); // time: 0s [13], r: 999979 * 999983

    //n = BigUint::parse_bytes(b"100433627766186892221372630609062766858404681029709092356097",10).unwrap();
    // time: []s lenght[61], r: 618970019642690137449562111 * 162259276829213363391578010288127

    n = (701 * 719).to_biguint().unwrap(); // time: 0s

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

fn bench(n: &BigUint) {
    let s = Instant::now();
    let mut i = BigUint::one();
    while i < *n {
        i = i + BigUint::one();
    }
    let elapsed = s.elapsed();

    // for 1billion -> 24143ms = 40348612.00 it/s
    println!("Elapsed: {} ms",
             (elapsed.as_secs() * 1_000) + (elapsed.subsec_nanos() / 1_000_000) as u64);
}

fn approx_sqrt(number: &BigUint, iterations: usize) -> BigUint {
    let start: BigUint = number.clone();
    let mut approx = start.clone();

    for _ in 0..iterations {
        approx = (&approx + (&start / &approx)) / 2.to_biguint().unwrap();
    }

    return approx;
}

// probabilistic approch: is_prime?
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

// base^exponent mod modulus
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
    return result;
}

fn find_r_and_d(i: &BigUint) -> (u64, BigUint) {
    let mut d = i.clone();
    let mut r = 0;
    let two = BigUint::new(vec![2]);

    while d.is_even() {
        d = d / two.clone();
        r += 1;
    }
    return (r, d);
}

fn little_fermat(num: &BigUint) -> bool {
    let mut rng = rand::thread_rng();
    let mut random: BigUint = rng.gen_biguint_below(num);
    let result = mod_exp(&mut random, &mut (num - BigUint::one()), num);
    return result == BigUint::one();
}

fn factor(num: &BigUint) -> Vec<(BigUint, BigUint)> {
    let factors: Arc<Mutex<Vec<(BigUint, BigUint)>>> = Arc::new(Mutex::new(Vec::new()));
    let cpus = num_cpus::get();
    let step: BigUint = approx_sqrt(num, 1000) / (cpus.to_biguint().unwrap());
    let mut start = step.clone();
    let mut end = start.clone() + step.clone();
    let mut threads = Vec::new();

    if *num > 5000000.to_biguint().unwrap() {
        for _ in 0..cpus {
            // needs to copy the value, for lifetime
            let st = start.clone();
            let en = start + step.clone();
            let n = num.clone();
            let lock: Arc<Mutex<Vec<(BigUint, BigUint)>>> = factors.clone();
            let zero = BigUint::zero();
            threads.push(thread::spawn(move || {
                let mut index = st;
                while index <= en && index < n {
                    if n.clone() % index.clone() == zero && is_prime(&index) {
                        let d = n.clone() / index.clone();
                        let mut _lock = lock.lock().unwrap();
                        if !_lock.iter().any(|p| p.0 == index || p.1 == index) && is_prime(&d) {
                            println!("{} FOUND: {} and {}",
                                     Blue.bold().paint("[!]"),
                                     Green.bold().paint(index.to_str_radix(10)),
                                     Green.bold().paint(d.to_str_radix(10)));
                            _lock.push((index.clone(), d.clone()));
                        }

                    }
                    index = index + BigUint::one();
                }
            }));
            start = end.clone() + BigUint::one();
            end = start.clone() + step.clone();
        }

        for t in threads {
            handle(t.join());
        }
    } else {
        let mut i = BigUint::new(vec![3]);
        let sqrt = approx_sqrt(num, 10000);
        while i < sqrt {
            if num % i.clone() == BigUint::zero() && is_prime(&i) {
                let d = &*num / i.clone();
                let mut _lock = factors.lock().unwrap();
                if !_lock.iter().any(|p| p.0 == i || p.1 == i) && is_prime(&d) {
                    _lock.push((i.clone(), d.clone()));
                    println!("{} FOUND: {} and {}",
                             Blue.bold().paint("[!]"),
                             Green.bold().paint(i.to_str_radix(10)),
                             Green.bold().paint(d.to_str_radix(10)));

                }
            }
            i = i + BigUint::one();
        }
    }

    let guard = Arc::try_unwrap(factors).expect("Lock is owned");
    return guard.into_inner().expect("Error on acquiring lock");
}

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

fn is_prime(num: &BigUint) -> bool {
    let max = 10000;
    let sieve = primal::Sieve::new(max);
    if *num <= max.to_biguint().unwrap() && sieve.is_prime(num.to_usize().unwrap()) {
        return true;
    }
    return little_fermat(num) && miller_rabin(num, 128);
}
