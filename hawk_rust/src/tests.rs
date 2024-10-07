

// fn sample_test() {
//     let num_samples = 10;
//     let mut result: Vec<i8> = vec![];
//
//     let n = 256;
//     for i in 0..num_samples{
//         let seed = get_random_bytes(10);
//         let t = get_random_bytes(2*n).iter().map(|&x| modulo(x, 2)).collect();
//         let x = sample(&seed, t, n); 
//         x.iter().for_each(|&el| {
//             result.push(el);
//         });
//     }
//
//     println!("{} samples: {:?}",num_samples, result);
// }
