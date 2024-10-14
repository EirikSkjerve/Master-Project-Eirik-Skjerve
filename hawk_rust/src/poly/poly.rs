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

#[derive(Debug)]
struct Poly<T: Copy + Default> {
    degree: usize,
    coefficients: Vec<T>, // this can be different based on poly_type
    poly_type: PolyType,
}

impl <T: Copy + Default> Poly<T>{

    pub fn new(n: usize, coefficients: &Vec<T>, poly_type: PolyType) -> Self {
        Poly{
            degree: n,
            coefficients: coefficients.to_vec(),
            poly_type,
        }
    }


    // Helper function to perform both addition and subtraction
    fn operate<F>(self, rhs: &Poly<T>, op: F) -> Poly<T>
    where
        F: Fn(T, T) -> T,
    {

        // assert that degree and type of polynomials are equal
        // in theory, degree does not need to be equal, but for this implementation
        // I want to assume that they are predefined to the correct size

        assert_eq!(self.degree, rhs.degree);
        assert_eq!(self.poly_type, rhs.poly_type);

        // create empy coefficient vector to keep the result coefficients
        let mut res: Vec<T> = Vec::with_capacity(self.degree);

        // loop through i from 0 to n
        for i in 0..self.degree {
            // add result from operation 
            res.push(op(self.coefficients[i], rhs.coefficients[i]));
        }

        // return a newly created polynomial instance with the resulting 
        Poly::new(self.degree, &res, self.poly_type)
    }
}


impl<T> ops::Add<Poly<T>> for Poly<T> 
where
    T: ops::Add<Output = T> + Copy + Default,
{
    type Output = Poly<T>;

    fn add(self, rhs: Poly<T>) -> Poly<T> {
        self.operate(&rhs, |a, b| a + b)
    }
}

impl<T> ops::Sub<Poly<T>> for Poly<T> 
where
    T: ops::Sub<Output = T> + Copy + Default,
{
    type Output = Poly<T>;

    fn sub(self, rhs: Poly<T>) -> Poly<T> {
        self.operate(&rhs, |a, b| a - b)
    }
}

// TODO Implement == operator for polynomials for completeness

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_poly_make() {
        let a: Vec<u8> = vec![1,2,3,4];
        let b: Vec<i64> = vec![8,16,32,64];
        let new_poly_a = Poly::new(4, &a, PolyType::R);
        let new_poly_b = Poly::new(4, &b, PolyType::NTT);

        assert_eq!(new_poly_a.coefficients, vec![1,2,3,4]);
        assert_eq!(new_poly_a.degree, new_poly_a.coefficients.len());

        assert_eq!(new_poly_b.coefficients, vec![8,16,32,64]);
        assert_eq!(new_poly_b.degree, new_poly_b.coefficients.len());
    }

    #[test]
    fn test_poly_add() {
        let a = Poly::new(4, &vec![1,2,3,4], PolyType::R);
        let b = Poly::new(4, &vec![4,3,2,1], PolyType::R);

        let c = a+b;

        assert_eq!(c.coefficients, vec![5,5,5,5]);

    }

    #[test]
    fn test_poly_sub() {

        let a = Poly::new(4, &vec![1,2,3,4], PolyType::R);
        let b = Poly::new(4, &vec![4,3,2,1], PolyType::R);

        let c = a-b;

        assert_eq!(c.coefficients, vec![-3,-1,1,3]);

    }
    //
    // #[test]
    // fn test_poly_ring_mult() {
    //
    // }
}
