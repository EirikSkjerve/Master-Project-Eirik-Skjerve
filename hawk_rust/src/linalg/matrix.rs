
// struct of a square matrix
#[derive(Debug)]
struct Matrix<T: Copy + Default> {
    degree: usize,
    entries: Vec<Vec<T>>,
}

impl <T: Copy + Default> Matrix<T> {
    ///
    /// initialize a new matrix given entries
    /// entries need to be square
    pub fn new(entries: &Vec<Vec<T>>) -> Matrix<T> {
        for e in entries {
            assert_eq!(entries.len(), e.len());
        }

        Matrix {
            degree: entries.len(),
            entries: entries.to_vec()
        }
    }

    // pub fn matmul(self, rhs: Matrix<T>) -> Matrix<T>
    // {
    //     let mut res: Vec<Vec<T>> = vec![Vec::with_capacity(self.degree); self.degree];
    //
    //     let mut sum: T = 0 as T;
    //     for i in 0..self.degree {
    //         for j in 0..self.degree {
    //         sum = 0 as T;
    //         for k in 0..self.degree {
    //             sum += (self.entries[i][k] * rhs.entries[k][j]);
    //         }
    //         }
    //     }       
    //
    //
    //     self
    // }
}
