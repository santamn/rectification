use nalgebra::{RealField, Unit, Vector2, convert};
use std::ops::Mul;

pub(super) trait EdgeDisplacements<T> {
    fn constrain_along(&self, e_r: &Unit<Vector2<T>>) -> ConstrainedEdgeDisplacements<T>;
}

impl<T: RealField + Copy> EdgeDisplacements<T> for (Vector2<T>, Vector2<T>) {
    /// dir方向についての長さを変化させないように両端の変位を制約する
    #[inline]
    fn constrain_along(&self, dir: &Unit<Vector2<T>>) -> ConstrainedEdgeDisplacements<T>
    where
        T: RealField + Copy,
    {
        let (dr_1, dr_2) = self;
        let dv = dir.as_ref() * dir.as_ref().dot(&(dr_1 - dr_2)) * convert::<_, T>(0.5);
        ConstrainedEdgeDisplacements(dr_1 - dv, dr_2 + dv)
    }
}

/// 棒の長さを変化させないような両端の変位を表す型
#[derive(Clone, Copy, Debug)]
pub(super) struct ConstrainedEdgeDisplacements<T>(Vector2<T>, Vector2<T>);

impl<T> ConstrainedEdgeDisplacements<T>
where
    T: RealField + Copy,
{
    #[inline]
    pub(super) fn into_inner(self) -> (Vector2<T>, Vector2<T>) {
        (self.0, self.1)
    }
}

impl<T> Mul<T> for &ConstrainedEdgeDisplacements<T>
where
    T: RealField + Copy,
{
    type Output = ConstrainedEdgeDisplacements<T>;

    fn mul(self, rhs: T) -> Self::Output {
        ConstrainedEdgeDisplacements(self.0 * rhs, self.1 * rhs)
    }
}

impl<T> Mul<T> for ConstrainedEdgeDisplacements<T>
where
    T: RealField + Copy,
{
    type Output = ConstrainedEdgeDisplacements<T>;

    fn mul(self, rhs: T) -> Self::Output {
        ConstrainedEdgeDisplacements(self.0 * rhs, self.1 * rhs)
    }
}
