use nalgebra::{RealField, Scalar, Vector2};
use std::ops::{Add, Deref, DerefMut, Sub};

/// 制約条件を適用していない変位を表す型
#[derive(Clone, Copy, Debug)]
pub(crate) struct Raw<T: Scalar>(Vector2<T>);

impl<T: Scalar> From<Vector2<T>> for Raw<T> {
    fn from(v: Vector2<T>) -> Self {
        Self(v)
    }
}

impl<T: RealField> Add for Raw<T> {
    type Output = Raw<T>;

    fn add(self, rhs: Self) -> Self::Output {
        Raw(self.0 + rhs.0)
    }
}

impl<T: RealField> Sub for Raw<T> {
    type Output = Raw<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        Raw(self.0 - rhs.0)
    }
}

impl<T: RealField> Add<&Raw<T>> for &Raw<T> {
    type Output = Raw<T>;

    fn add(self, rhs: &Raw<T>) -> Self::Output {
        Raw(&self.0 + &rhs.0)
    }
}

impl<T: RealField> Sub<&Raw<T>> for &Raw<T> {
    type Output = Raw<T>;

    fn sub(self, rhs: &Raw<T>) -> Self::Output {
        Raw(&self.0 - &rhs.0)
    }
}

/// 制約条件を満たすような変位を表す型
#[derive(Clone, Copy, Debug)]
pub(super) struct Constrained<T: Scalar>(Vector2<T>);

impl<T: Scalar> Deref for Constrained<T> {
    type Target = Vector2<T>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Scalar> DerefMut for Constrained<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
