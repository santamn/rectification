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
        let [dr] = forces;
        let tentative = self.current + dr;

        // 壁との衝突を考慮して変位を適用
        self.current += match (Bottom::<T>::f(&tentative.x)..=Top::<T>::f(&tentative.x))
            .position_of(&tentative.y)
        {
            // 上壁と衝突する場合
            Above => {
                let collision_x = binary_search_root(Top::<T>::f, &self.current.x, &tentative.x);
                let to_channel = tentative - Point2::new(collision_x, Top::<T>::f(&collision_x));
                to_channel + Top::<T>::reflect_at(&collision_x, &(dr - to_channel))
            }
            // 境界内の場合
            Within => dr,
            // 下壁と衝突する場合
            Below => {
                let collision_x = binary_search_root(Bottom::<T>::f, &self.current.x, &tentative.x);
                let to_channel = tentative - Point2::new(collision_x, Bottom::<T>::f(&collision_x));
                to_channel + Bottom::<T>::reflect_at(&collision_x, &(dr - to_channel))
            }
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
