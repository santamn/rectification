use nalgebra::{Point2, RealField, Scalar, Unit, Vector2, convert};
use std::marker::PhantomData;

// 境界条件 ω(x) の上/下を区別するためのマーカー
pub(crate) struct Top<T>(PhantomData<T>); // y = ω(x)
pub(crate) struct Bottom<T>(PhantomData<T>); // y = -ω(x)

pub(crate) trait Boundary<T: Scalar> {
    /// 境界の関数
    fn f(x: &T) -> T;
    /// (x, ±ω(x)) においてベクトル v を反射させたものを返す
    fn reflect_at(x: &T, v: &Vector2<T>) -> Vector2<T>;
    /// 点が境界の外側にあるかどうかを判定する
    fn is_outside(p: &Point2<T>) -> bool;
}

impl<T: RealField + Copy> Boundary<T> for Top<T> {
    #[inline]
    fn f(x: &T) -> T {
        omega(x)
    }

    #[inline]
    fn reflect_at(x: &T, v: &Vector2<T>) -> Vector2<T> {
        let n = Unit::new_normalize(Vector2::new(omega_prime(x), -T::one()));
        v - n.as_ref() * convert::<_, T>(2.0) * n.as_ref().dot(v)
    }

    #[inline]
    fn is_outside(p: &Point2<T>) -> bool {
        p.y > omega(&p.x)
    }
}

impl<T: RealField + Copy> Boundary<T> for Bottom<T> {
    #[inline]
    fn f(x: &T) -> T {
        -omega(x)
    }

    #[inline]
    fn reflect_at(x: &T, v: &Vector2<T>) -> Vector2<T> {
        let n = Unit::new_normalize(Vector2::new(-omega_prime(x), T::one()));
        v - n.as_ref() * convert::<_, T>(2.0) * n.as_ref().dot(v)
    }

    #[inline]
    fn is_outside(p: &Point2<T>) -> bool {
        p.y < -omega(&p.x)
    }
}

/// ω(x) = sin(2πx) + 0.25sin(4πx) + 1.12
fn omega<T>(x: &T) -> T
where
    T: RealField + Copy,
{
    let (s, c) = (T::two_pi() * *x).sin_cos();
    s + convert::<_, T>(0.5) * s * c + convert::<_, T>(1.12)
}

/// チャネル境界の傾き ω'(x) = 2πcos(2πx) + πcos(4πx)
fn omega_prime<T>(x: &T) -> T
where
    T: RealField + Copy,
{
    let c = (T::two_pi() * *x).cos();
    T::two_pi() * c * (c + T::one()) - T::pi()
}
