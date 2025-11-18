use nalgebra::{Point2, RealField, Scalar, Vector2, convert};
use rand::{Rng, SeedableRng, rngs::SmallRng};
use rand_distr::{Distribution, StandardNormal, uniform::SampleUniform};
use rayon::prelude::*;
use std::f64::consts::SQRT_2;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::{Add, Sub};

type Real = f64; // 計算の精度を決める型

const PARTICLES: u64 = 30_000; // アンサンブル平均に用いる粒子数
const STEPS: u64 = 100_000; // シミュレーションの時間ステップ数
const DELTA_T: Real = 0.000_1; // 時間刻み幅
const TIME: Real = STEPS as Real * DELTA_T; // 総シミュレーション時間

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
    let start = std::time::Instant::now();

    let mut mu_writer = BufWriter::new(File::create("data/mu.dat").unwrap());
    let mut d_writer = BufWriter::new(File::create("data/d_eff.dat").unwrap());
    let mut alpha_writer = BufWriter::new(File::create("data/alpha.dat").unwrap());

    for i in 0..=100 {
        let f_x = i as Real;

        let (mean, mean_square) =
            simulate_brownian_motion(STEPS, PARTICLES, DELTA_T, Vector2::new(f_x, 0.0));
        let (mean_rev, mean_square_rev) =
            simulate_brownian_motion(STEPS, PARTICLES, DELTA_T, Vector2::new(-f_x, 0.0));

        let mu = nonlinear_mobility(mean / TIME, f_x);
        let mu_rev = nonlinear_mobility(mean_rev / TIME, -f_x);
        let d_eff = effective_diffusion(mean, mean_square, TIME);
        let d_eff_rev = effective_diffusion(mean_rev, mean_square_rev, TIME);
        let alpha_val = alpha(mu, mu_rev);

        writeln!(mu_writer, "{} {} {}", f_x, mu, mu_rev).unwrap();
        writeln!(d_writer, "{} {} {}", f_x, d_eff, d_eff_rev).unwrap();
        writeln!(alpha_writer, "{} {}", f_x, alpha_val).unwrap();
    }

    println!("Elapsed: {:.2?}", start.elapsed());
}

/// ブラウン運動のシミュレーションを実行し、粒子の平均変位と平均二乗変位を返す
fn simulate_brownian_motion<T>(steps: u64, particles: u64, delta_t: T, f: Vector2<T>) -> (T, T)
where
    T: RealField + SampleUniform + Copy,
    StandardNormal: Distribution<T>,
{
    let sqrt_delta_t = delta_t.sqrt();
    let displacements = (0..particles)
        .into_par_iter() // 各粒子のシミュレーションを並列化
        .map(|i| {
            let mut rng = SmallRng::seed_from_u64(i);
            let particle = Particle::new(random_point(&mut rng));

            (0..steps)
                .map(|_| {
                    // 微小時間に粒子に加わる外力F + ブラウン運動
                    f * delta_t
                        + Vector2::new(
                            rng.sample::<T, _>(StandardNormal),
                            rng.sample::<T, _>(StandardNormal),
                        ) * convert::<_, T>(SQRT_2)
                            * sqrt_delta_t
                })
                .fold(particle, |acc, dr| {
                    let tentative_pos: Point2<T> = acc.current + dr;

                    acc + if omega::<Plus, T>(tentative_pos.x) < tentative_pos.y {
                        // 粒子が天井と衝突する場合
                        reflected_vector::<Plus, T>(&acc.current, &dr)
                    } else if tentative_pos.y < omega::<Minus, T>(tentative_pos.x) {
                        // 粒子が床と衝突する場合
                        reflected_vector::<Minus, T>(&acc.current, &dr)
                    } else {
                        // 衝突なし
                        dr
                    }
                })
                .displacement()
        })
        .collect::<Vec<_>>();

    (
        displacements.iter().fold(T::zero(), |a, &x| a + x) / convert::<_, T>(particles as f64),
        displacements.iter().fold(T::zero(), |a, &x| a + x * x) / convert::<_, T>(particles as f64),
    )
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

/// チャネル境界における内側への単位法線ベクトル
fn normal_vector<B, T>(x: T) -> Vector2<T>
where
    B: Boundary,
    T: RealField + Copy,
{
    Vector2::new(omega_prime::<Plus, T>(x), -B::sign::<T>()).normalize()
}

/// 0 <= x <= 1 の範囲でチャネル内のランダムな座標を生成
fn random_point<R, T>(rng: &mut R) -> Point2<T>
where
    R: rand::Rng + ?Sized,
    T: RealField + SampleUniform + Copy,
{
    let x = rng.random_range(T::zero()..T::one());
    let y = rng.random_range(omega::<Minus, T>(x)..omega::<Plus, T>(x));
    Point2::new(x, y)
}

/// 2分探索で方程式 f(x)=0 の根を求める
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

/// 粒子の移動ベクトルを壁で反射させたものを返す
fn reflected_vector<B, T>(current: &Point2<T>, dr: &Vector2<T>) -> Vector2<T>
where
    B: Boundary,
    T: RealField + Copy,
{
    // 粒子と壁の衝突点を求める
    let x = binary_search_root(
        |&x| current.y + (current.x - x) * dr.y / dr.x - omega::<B, T>(x),
        current.x + dr.x,
        current.x,
    );
    // 衝突点での法線ベクトルを求める
    let n = normal_vector::<B, T>(x);
    // 移動ベクトルを壁で反射したものを返す
    dr - n * convert::<_, T>(2.0) * n.dot(&(current + dr - Point2::new(x, omega::<B, T>(x))))
}

/// 非線形移動度 μ(f) = ⟨v⟩/f
fn nonlinear_mobility(mean_speed: Real, force: Real) -> Real {
    mean_speed / force
}

/// 有効拡散係数 D_eff = (⟨x²⟩ - ⟨x⟩²)/2t
fn effective_diffusion(
    mean_displacement: Real,
    mean_square_displacement: Real,
    time: Real,
) -> Real {
    (mean_square_displacement - mean_displacement * mean_displacement) / (2.0 * time)
}

/// 整流尺度 α = |μ(f) - μ(-f)| / (μ(f) + μ(-f))
fn alpha(mu: Real, mu_rev: Real) -> Real {
    (mu - mu_rev).abs() / (mu + mu_rev)
}
