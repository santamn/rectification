use nalgebra::{Point2, RealField, Scalar, Vector2, convert};
use rand::{Rng, SeedableRng, rngs::SmallRng};
use rand_distr::{StandardNormal, uniform::SampleUniform};
use rayon::prelude::*;
use std::ops::{Add, Sub};

type Real = f64; // 計算の精度を決める型

const PARTICLES: u64 = 30_000;
const STEPS: u64 = 100_000;
const SQRT_DELTA_T: Real = 0.1;
const F: Real = 1.0;

// 境界条件 ω(x) の上/下を区別するためのマーカー
struct Plus; // y = ω(x)
struct Minus; // y = -ω(x)

trait Boundary {
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

#[derive(Clone, Copy, Debug)]
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
        .into_par_iter() // 各粒子のシミュレーションを並列化
        .map(|i| {
            let mut rng = SmallRng::seed_from_u64(i);
            let particle = Particle::new(random_point(&mut rng));

            (0..STEPS)
                .map(|_| {
                    // 微小時間に粒子に加わる外力F + ブラウン運動
                    Vector2::new(
                        rng.sample::<Real, _>(StandardNormal) * SQRT_DELTA_T + F,
                        rng.sample::<Real, _>(StandardNormal) * SQRT_DELTA_T,
                    )
                })
                .fold(particle, |acc, dr| {
                    let tentative_pos = acc.current + dr;

                    acc + if omega::<Plus, Real>(tentative_pos.x) < tentative_pos.y {
                        // 粒子が天井と衝突する場合
                        reflected_vector::<Plus, Real>(&acc.current, &dr)
                    } else if tentative_pos.y < omega::<Minus, Real>(tentative_pos.x) {
                        // 粒子が床と衝突する場合
                        reflected_vector::<Minus, Real>(&acc.current, &dr)
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
    Vector2::new(omega_prime::<Plus, T>(x), -B::sign::<T>()).normalize()
}

/// チャネル内のランダムな初期位置を生成
fn random_point<R, T>(rng: &mut R) -> Point2<T>
where
    R: rand::Rng + ?Sized,
    T: RealField + SampleUniform + Copy,
{
    let x = rng.random_range(T::zero()..T::one());
    let y = rng.random_range(omega::<Minus, T>(x)..omega::<Plus, T>(x));
    Point2::new(x, y)
}

/// ニュートン法で方程式 f(x)=0 の根を求める
// FIXME: 一定時間内に計算が収束しない
fn newton_root<F, G, T>(f: F, df: G, x0: T) -> Point2<T>
where
    F: Fn(T) -> T,
    G: Fn(T) -> T,
    T: RealField + Copy,
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
    // 移動ベクトルを壁で反射したものを返す
    dr - n * convert::<_, T>(2.0) * n.dot(&(current + dr - intersection))
}

/// 非線形移動度
fn nonlinear_mobility(mean_displacement: Real, force: Real, steps: Real, delta_t: Real) -> Real {
    mean_displacement / force * steps * delta_t
}

/// 有効拡散係数
fn effective_diffusion(
    mean_displacement: Real,
    mean_square_displacement: Real,
    steps: Real,
    delta_t: Real,
) -> Real {
    (mean_square_displacement - mean_displacement * mean_displacement) / (2.0 * steps * delta_t)
}
