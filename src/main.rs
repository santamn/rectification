use nalgebra::{Point2, RealField, Scalar, Vector2, convert};
use rand::{Rng, SeedableRng, rngs::SmallRng};
use rand_distr::{StandardNormal, uniform::SampleUniform};
use rayon::prelude::*;
use std::ops::{Add, Sub};

type Real = f64; // 計算の精度を決める型

const PARTICLES: u64 = 30_000;
const STEPS: u64 = 10_000;
const SQRT_DELTA_T: Real = 0.1;
const F: Real = 1.0;

// 境界条件 ω(x) の上部・下部を区別するためのマーカー
struct Ceiling; // 上部: y = ω(x)
struct Floor; // 下部: y = -ω(x)

trait Boundary {
    /// 境界 y = ±ω(x) の符号を表す
    fn sign<T: RealField>() -> T;
}

impl Boundary for Ceiling {
    #[inline]
    fn sign<T: RealField>() -> T {
        T::one()
    }
}

impl Boundary for Floor {
    #[inline]
    fn sign<T: RealField>() -> T {
        -T::one()
    }
}

#[derive(Copy, Clone, Debug)]
struct Particle<T: Scalar> {
    current: Point2<T>,
    initial: Point2<T>,
}

impl<T> Add<Vector2<T>> for Particle<T>
where
    T: RealField,
{
    type Output = Self;

    fn add(mut self, rhs: Vector2<T>) -> Self::Output {
        self.current += rhs;
        self
    }
}

impl<T> Sub<Vector2<T>> for Particle<T>
where
    T: RealField,
{
    type Output = Self;

    fn sub(mut self, rhs: Vector2<T>) -> Self::Output {
        self.current -= rhs;
        self
    }
}

impl<T> Particle<T>
where
    T: RealField + Copy,
{
    fn new(initial: Point2<T>) -> Self {
        Self {
            current: initial,
            initial,
        }
    }

    fn displacement(&self) -> T {
        self.current.x - self.initial.x
    }
}

fn main() {
    let displacements = (0..PARTICLES)
        .into_par_iter()
        .map(|i| {
            let mut rng = SmallRng::seed_from_u64(i);
            let particle = Particle::new(random_point(&mut rng));

            // 各粒子について 外力F + ブラウン運動 のシミュレーションを行う
            (0..STEPS)
                .map(|_| {
                    // 微小時間に粒子に加わる外力 + ブラウン運動
                    Vector2::new(
                        rng.sample::<Real, _>(StandardNormal) * SQRT_DELTA_T + F,
                        rng.sample::<Real, _>(StandardNormal) * SQRT_DELTA_T,
                    )
                })
                .fold(particle, |acc, dr| {
                    let tentative_pos = acc.current + dr;

                    acc + if omega::<Ceiling, Real>(tentative_pos.x) < tentative_pos.y {
                        // 粒子が天井と衝突する場合
                        reflected_vector::<Ceiling, Real>(&acc.current, &dr)
                    } else if tentative_pos.y < omega::<Floor, Real>(tentative_pos.x) {
                        // 粒子が床と衝突する場合
                        reflected_vector::<Floor, Real>(&acc.current, &dr)
                    } else {
                        // 衝突なし
                        dr
                    }
                })
                .displacement()
        })
        .collect::<Vec<_>>();

    let mean = displacements.iter().sum::<Real>() / PARTICLES as Real;
    let mean_square = displacements.iter().map(|x| x * x).sum::<Real>() / PARTICLES as Real;

    println!(
        "μ(f): {}",
        nonlinear_mobility(mean, F, STEPS as Real, SQRT_DELTA_T * SQRT_DELTA_T)
    );
    println!(
        "D_eff: {}",
        effective_diffusion(
            mean,
            mean_square,
            STEPS as Real,
            SQRT_DELTA_T * SQRT_DELTA_T
        )
    );
}

/// チャネル境界 ω(x) = sin(2πx) + 0.25sin(4πx) + 1.12
fn omega<B, T>(x: T) -> T
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

/// チャネル境界における内側への法線ベクトル
fn normal_vector<B, T>(x: T) -> Vector2<T>
where
    B: Boundary,
    T: RealField + Copy,
{
    Vector2::new(omega_prime::<Ceiling, T>(x), -B::sign::<T>()).normalize()
}

/// チャネル内のランダムな初期位置を生成
fn random_point<R, T>(rng: &mut R) -> Point2<T>
where
    R: rand::Rng + ?Sized,
    T: RealField + SampleUniform + Copy,
{
    let x = rng.random_range(T::zero()..T::one());
    let y = rng.random_range(omega::<Floor, T>(x)..omega::<Ceiling, T>(x));
    Point2::new(x, y)
}

/// ニュートン法で方程式 f(x)=0 の根を求める
// FIXME
fn newton_root<F, G, T>(f: F, df: G, x0: T) -> Point2<T>
where
    T: RealField + Copy,
    F: Fn(T) -> T,
    G: Fn(T) -> T,
{
    std::iter::successors(Some(x0), |&x| Some(x - f(x) / df(x)))
        .find_map(|x| (f(x).abs() <= convert(0.000_1)).then_some(Point2::new(x, f(x))))
        .unwrap()
}

/// 粒子の移動ベクトルを壁で反射させたものを返す
fn reflected_vector<B, T>(current: &Point2<T>, dr: &Vector2<T>) -> Vector2<T>
where
    B: Boundary,
    T: RealField + Copy,
{
    // 粒子と壁の衝突点を求める
    let intersection = newton_root(
        |x| current.y + dr.y * (x - current.x) / dr.x - omega::<B, T>(x),
        |x| dr.y / dr.x - omega_prime::<B, T>(x),
        current.x,
    );
    // 衝突点での法線ベクトルを求める
    let n = normal_vector::<B, T>(intersection.x);

    dr - n * convert::<_, T>(2.0) * n.dot(&(current + dr - intersection))
}

/// 非線形移動度
fn nonlinear_mobility<T>(mean_displacement: T, force: T, steps: T, delta_t: T) -> T
where
    T: RealField,
{
    mean_displacement / force * steps * delta_t
}

/// 有効拡散係数
fn effective_diffusion<T>(
    mean_displacement: T,
    mean_square_displacement: T,
    steps: T,
    delta_t: T,
) -> T
where
    T: RealField + Copy,
{
    (mean_square_displacement - mean_displacement * mean_displacement)
        / (convert::<_, T>(2.0) * steps * delta_t)
}
