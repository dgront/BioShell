use std::ops::{Index, IndexMut};
use crate::vec3::Vec3;

/// Represents a 3x3 matrix, e.g. for linear transformations.
///
/// Internally the matrix elements are stored as an `[f64; 9]` array.
///
/// # Example
/// ```rust
/// use bioshell_numerical::matrix::Matrix3x3;
/// use bioshell_numerical::Vec3;
/// let vx = Vec3::new(1.0, 0.0, 0.0);
/// let vy = Vec3::new(0.0, 1.0, 0.0);
/// let vz = Vec3::new(1.0, 0.0, 1.0);
/// let unit_mtx = Matrix3x3::from_column_vectors(&vx, &vy, &vz);
///
/// assert_eq!(unit_mtx[0], 1.0);
/// assert_eq!(unit_mtx[4], 1.0);
/// assert_eq!(unit_mtx[8], 1.0);
/// ```
pub struct Matrix3x3
{
    _array: [f64; 9],
}

impl Index<usize> for Matrix3x3
{
    type Output = f64;
    fn index<'a>(&'a self, i: usize) -> &'a f64
    {
        &self._array[i]
    }
}

impl IndexMut<usize> for Matrix3x3
{
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut f64
    {
        &mut self._array[i]
    }
}

impl Matrix3x3
{
    pub fn show(&self)
    {
        println!("[[{} {} {}]", self._array[0], self._array[1], self._array[2]);
        println!(" [{} {} {}]", self._array[3], self._array[4], self._array[5]);
        println!(" [{} {} {}]]", self._array[6], self._array[7], self._array[8]);
    }
    pub fn new(m: [f64; 9]) -> Self
    {
        Matrix3x3 { _array: m }
    }

    pub fn from_column_values(cx_x: f64, cx_y: f64, cx_z: f64,
                              cy_x: f64, cy_y: f64, cy_z: f64,
                              cz_x: f64, cz_y: f64, cz_z: f64) -> Self
    {
        Self::new([
            cx_x, cx_y, cx_z,
            cy_x, cy_y, cy_z,
            cz_x, cz_y, cz_z,
        ])
    }

    pub fn from_column_vectors(cx: &Vec3, cy: &Vec3, cz: &Vec3) -> Self
    {
        let cx1 = cx.normalized();
        let cy1 = cy.normalized();
        let cz1 = cz.normalized();

        Self::new([
            cx1.x, cx1.y, cx1.z,
            cy1.x, cy1.y, cy1.z,
            cz1.x, cz1.y, cz1.z,
        ])
    }

    pub fn from_row_vectors(rx: &Vec3, ry: &Vec3, rz: &Vec3) -> Self
    {
        let rx1 = rx.normalized();
        let ry1 = ry.normalized();
        let rz1 = rz.normalized();

        Self::new([
            rx1.x, ry1.x, rz1.x,
            rx1.y, ry1.y, rz1.y,
            rx1.z, ry1.z, rz1.z,
        ])
    }
}