use nalgebra::{Point2, RealField, Vector2, convert};
use std::ops::RangeInclusive;

// 境界条件 ω(x) の上/下を区別するためのマーカー
pub(crate) struct Plus; // y = ω(x)
pub(crate) struct Minus; // y = -ω(x)

pub(crate) trait Boundary {
    /// 境界 y = ±ω(x) の符号を表す
    fn sign<T: RealField>() -> T;
}

impl Boundary for Plus {
    #[inline]
    fn sign<T: RealField>() -> T {
        T::one()
    }
}

impl Boundary for Minus {
    #[inline]
    fn sign<T: RealField>() -> T {
        -T::one()
    }
}

/// チャネル境界 ω(x) = sin(2πx) + 0.25sin(4πx) + 1.12
pub(crate) fn omega<B, T>(x: T) -> T
where
    B: Boundary,
    T: RealField + Copy,
{
    let (s, c) = (T::two_pi() * x).sin_cos();
    B::sign::<T>() * (s + convert::<_, T>(0.5) * s * c + convert::<_, T>(1.12))
}

/// チャネル境界の傾き ω'(x) = 2πcos(2πx) + πcos(4πx)
fn omega_prime<B, T>(x: T) -> T
where
    B: Boundary,
    T: RealField + Copy,
{
    let c = (T::two_pi() * x).cos();
    B::sign::<T>() * (T::two_pi() * c * (c + T::one()) - T::pi())
}

/// チャネル境界における内側への単位法線ベクトル
pub(crate) fn normal<B, T>(x: T) -> Vector2<T>
where
    B: Boundary,
    T: RealField + Copy,
{
    Vector2::new(omega_prime::<Plus, T>(x), -B::sign::<T>()).normalize()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Zone {
    Above,  // 上部外側
    Within, // 内部
    Below,  // 下部外側
}

impl Zone {
    pub(crate) fn of_point<T>(range: RangeInclusive<T>, point: &Point2<T>) -> Self
    where
        T: RealField + Copy,
    {
        if *range.end() < point.y {
            Zone::Above
        } else if point.y < *range.start() {
            Zone::Below
        } else {
            Zone::Within
        }
    }
}
