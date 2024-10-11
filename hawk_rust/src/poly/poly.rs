use std::ops;

// normal polynomials in our ring Z[x]/(x^n + 1)

#[derive(Debug, PartialEq)]
enum PolyType {
    R, // our ring Z[x] / (x^n + 1)
    NTT, // ntt representation
    FFT, // fft representation
    Q, // polynomials over field Q[x]
    VEC, // only vector representation
}

struct Poly {
    degree: usize,
    coefficients: Vec<i64>, // this can be different based on poly_type
    poly_type: PolyType,
}

impl Poly{
    pub fn new(n: usize, coefficients: &Vec<i64>, poly_type: PolyType) -> Poly {
        Poly{
            degree: n,
            coefficients: coefficients.to_vec(),
            poly_type,
        }
    }
}

impl ops::Add<Poly> for Poly {
    type Output = Poly;

    /// implements addition for polynomials
    /// no modulation yet
    fn add(self, rhs: Poly) -> Poly {
        // check that these are compatible
        assert_eq!(self.degree, rhs.degree);
        assert_eq!(self.poly_type, rhs.poly_type);

        // sum together coefficient for coefficient
        let coefficient_sum: Vec<i64> = (0..self.degree)
            .into_iter()
            .map(|i| self.coefficients[i] + rhs.coefficients[i])
            .collect();

        // return a new polynomial struct
        Poly{
            degree: self.degree,
            coefficients: coefficient_sum,
            poly_type: self.poly_type,
        }
    }
}

// TODO Implement == operator for polynomials for completeness

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_poly_make() {
        let new_poly = Poly::new(4, &vec![1,2,3,4], PolyType::R);
        assert_eq!(vec![1,2,3,4], new_poly.coefficients);
    }

    #[test]
    fn test_poly_add() {
        let a = Poly::new(4, &vec![1,2,3,4], PolyType::R);
        let b = Poly::new(4, &vec![4,3,2,1], PolyType::R);

        let c = a+b;

        assert_eq!(c.coefficients, vec![5,5,5,5]);

    }
    //
    // #[test]
    // fn test_poly_sub() {
    //
    // }
    //
    // #[test]
    // fn test_poly_ring_mult() {
    //
    // }
}
