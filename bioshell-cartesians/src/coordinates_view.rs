use crate::Coordinates;

/// Stateless immutable view of coordinates
pub struct CoordinatesView<'a> {
    pub points: &'a Coordinates,
}