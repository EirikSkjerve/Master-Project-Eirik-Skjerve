use std::ops::{Add, Mul};

// Define a struct for the matrix
#[derive(Debug, Clone)]
struct Matrix<T> {
    size: usize,
    data: Vec<Vec<T>>,
}

// Implement methods for the matrix
impl<T> Matrix<T> {
    // Constructor for square matrix
    fn new(size: usize, initial_value: T) -> Matrix<T>
    where
        T: Clone,
    {
        let data = vec![vec![initial_value; size]; size];
        Matrix { size, data }
    }

    // Get the size of the matrix
    fn size(&self) -> usize {
        self.size
    }

    // Get an element of the matrix
    fn get(&self, row: usize, col: usize) -> &T {
        &self.data[row][col]
    }

    // Set an element of the matrix
    fn set(&mut self, row: usize, col: usize, value: T) {
        self.data[row][col] = value;
    }
}

// Implement matrix addition
impl<T> Add for Matrix<T>
where
    T: Add<Output = T> + Clone,
{
    type Output = Matrix<T>;

    fn add(self, other: Matrix<T>) -> Matrix<T> {
        let mut result = Matrix::new(self.size, self.data[0][0].clone());
        for i in 0..self.size {
            for j in 0..self.size {
                result.set(i, j, self.get(i, j).clone() + other.get(i, j).clone());
            }
        }
        result
    }
}

// Implement matrix multiplication
impl<T> Mul for Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Clone + Default,
{
    type Output = Matrix<T>;

    fn mul(self, other: Matrix<T>) -> Matrix<T> {
        let mut result = Matrix::new(self.size, T::default());
        for i in 0..self.size {
            for j in 0..self.size {
                let mut sum = T::default();
                for k in 0..self.size {
                    sum = sum + (self.get(i, k).clone() * other.get(k, j).clone());
                }
                result.set(i, j, sum);
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test() {
        let mut matrix_a = Matrix::new(2, 0);
        matrix_a.set(0, 0, 1);
        matrix_a.set(0, 1, -1);
        matrix_a.set(1, 0, 2);
        matrix_a.set(1, 1, 0);

        let mut matrix_b = Matrix::new(2, 0);
        matrix_b.set(0, 0, 2);
        matrix_b.set(0, 1, 1);
        matrix_b.set(1, 0, 1);
        matrix_b.set(1, 1, -1);

        let matrix_sum = matrix_a.clone() + matrix_b.clone();
        let matrix_product = matrix_a * matrix_b;

        println!("Matrix A + B = {:?}", matrix_sum);
        println!("Matrix A * B = {:?}", matrix_product);
    }
}
