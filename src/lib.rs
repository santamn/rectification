use nalgebra::{Point2, RealField, Vector2, convert};
use rand::Rng;
use std::cmp::Ordering::{self, *};
use std::ops::ControlFlow::*;

pub mod boundary;
pub mod double;
pub mod single;

pub use double::Diparticle;
pub use single::Monoparticle;

pub trait Particle<T, const C: usize> {
    type Size;

    fn new<R: Rng>(rng: &mut R, size: Self::Size) -> Self;
    fn displacement(&self) -> T;
    fn apply_forces(&mut self, forces: [Vector2<T>; C]);
}

/// 2分探索で方程式 f(x)=0 の根を求める
///
/// 引数には f(a)*f(b) < 0 を満たす a, b を与えること
fn binary_search_root<F, T>(f: F, a: T, b: T) -> T
where
    F: Fn(&T) -> T,
    T: RealField + Copy,
{
    std::iter::successors(
        Some(if f(&a).is_positive() { (a, b) } else { (b, a) }),
        |&(h, l)| {
            let m = (h + l) * convert(0.5);
            Some(if f(&m).is_positive() { (m, l) } else { (h, m) })
        },
    )
    .nth(20)
    .map(|(a, b)| (a + b) * convert(0.5))
    .unwrap()
}

/// 粒子の移動ベクトル dr を壁で反射させたものを返す
fn reflect<F, G, T>(boundary: F, normal_vec: G, current: &Point2<T>, dr: &Vector2<T>) -> Vector2<T>
where
    F: Fn(T) -> T,
    G: Fn(T) -> Vector2<T>,
    T: RealField + Copy,
{
    // 境界との交点を求める
    let x = binary_search_root(
        |&x| current.y + (current.x - x) * dr.y / dr.x - boundary(x),
        current.x + dr.x,
        current.x,
    );
    // 法線ベクトルを使って反射ベクトルを計算
    let n = normal_vec(x);
    // 仮想的な変位 dr を境界で反射したベクトルを返す
    dr - n * n.dot(&(current + dr - Point2::new(x, boundary(x)))) * convert::<_, T>(2.0)
}

/// f,g　の根が (0, 1) に存在するという仮定のもとで、二分探索で関数 f,g の根を大小比較する
///
/// f の根が g の根より大きければ Greater、小さければ Less、十分近ければ Equal を返す
fn binary_cmp_root<F, G, T>(f: F, g: G) -> Ordering
where
    F: Fn(&T) -> T,
    G: Fn(&T) -> T,
    T: RealField + Copy,
{
    (0..20)
        .try_fold(
            if f(&T::zero()).is_positive() {
                (T::zero(), T::one())
            } else {
                (T::one(), T::zero())
            },
            |(h, l), _| {
                let m = (h + l) * convert(0.5);
                match (f(&m).is_positive(), g(&m).is_positive()) {
                    (true, true) => Continue((m, l)),
                    (false, false) => Continue((h, m)),
                    (true, false) => Break(Greater),
                    (false, true) => Break(Less),
                }
            },
        )
        .break_value()
        .unwrap_or(Equal)
}
