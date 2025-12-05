use crate::{
    Particle,
    boundary::{
        Minus, Plus,
        Zone::{self, *},
        normal, omega,
    },
    reflect,
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
    T: RealField + SampleUniform + Copy,
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
        self.current += match Zone::of_point(
            omega::<Minus, T>(tentative.x)..=omega::<Plus, T>(tentative.x),
            &tentative.y,
        ) {
            // 上壁と衝突する場合
            Above => reflect(omega::<Plus, T>, normal::<Plus, T>, &self.current, &dr),
            // 境界内の場合
            Within => dr,
            // 下壁と衝突する場合
            Below => reflect(omega::<Minus, T>, normal::<Minus, T>, &self.current, &dr),
        };
    }
}

fn random_point<R, T>(rng: &mut R) -> Point2<T>
where
    R: rand::Rng + ?Sized,
    T: RealField + SampleUniform + Copy,
{
    let x = rng.random_range(T::zero()..T::one());
    let y = rng.random_range(omega::<Minus, T>(x)..omega::<Plus, T>(x));
    Point2::new(x, y)
}
