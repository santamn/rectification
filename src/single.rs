use crate::{
    Particle, RangeExt,
    RelativePos::*,
    binary_search_root,
    boundary::{Bottom, Boundary, Top},
};
use nalgebra::{Point2, RealField, Scalar, Vector2};
use rand::Rng;
use rand_distr::uniform::SampleUniform;

#[derive(Clone, Copy, Debug)]
pub struct Monoparticle<T: Scalar> {
    initial: Point2<T>, // 初期位置
    current: Point2<T>, // 現在位置
}

impl<T> Particle<T, 1> for Monoparticle<T>
where
    T: RealField + SampleUniform + Copy + Send + Sync,
{
    type Size = ();

    fn new<R: Rng>(rng: &mut R, _size: ()) -> Self {
        let initial = random_point::<R, T>(rng);
        Self {
            initial,
            current: initial,
        }
    }

    #[inline]
    fn displacement(&self) -> T {
        self.current.x - self.initial.x
    }

    fn apply_forces(&mut self, forces: [Vector2<T>; 1]) {
        dbg!(&forces);
        dbg!(&self.current);

        let [dr] = forces;
        let tentative = self.current + dr;

        // 壁との衝突を考慮して変位を適用
        self.current += match (Bottom::<T>::f(&tentative.x)..=Top::<T>::f(&tentative.x))
            .position_of(&tentative.y)
        {
            // 上壁と衝突する場合
            Above => reflected_vector::<Top<_>, T>(&self.current, &dr),
            // 境界内の場合
            Within => dr,
            // 下壁と衝突する場合
            Below => reflected_vector::<Bottom<_>, T>(&self.current, &dr),
        };
    }
}

fn random_point<R, T>(rng: &mut R) -> Point2<T>
where
    R: rand::Rng + ?Sized,
    T: RealField + SampleUniform + Copy + Send + Sync,
{
    let x = rng.random_range(T::zero()..T::one());
    let y = rng.random_range(Bottom::<T>::f(&x)..Top::<T>::f(&x));
    Point2::new(x, y)
}

fn reflected_vector<B, T>(current: &Point2<T>, v: &Vector2<T>) -> Vector2<T>
where
    B: Boundary<T>,
    T: RealField + Copy,
{
    // 粒子と壁の衝突点を求める
    let collision_x = binary_search_root(
        |x| current.y + (current.x - *x) * v.y / v.x - B::f(x),
        &(current.x + v.x),
        &current.x,
    );
    // 現在位置から衝突点までのベクトル
    let to_collision_point = current + v - Point2::new(collision_x, B::f(&collision_x));
    // 反射を計算
    to_collision_point + B::reflect_at(&collision_x, &(v - to_collision_point))
}
