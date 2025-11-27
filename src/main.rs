use nalgebra::{RealField, Vector2, convert};
use rand::{SeedableRng, rngs::SmallRng};
use rand_distr::{Distribution, StandardNormal, uniform::SampleUniform};
use rayon::prelude::*;
use rectification::{Diparticle, Particle};
use std::f64::consts::SQRT_2;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::iter::Sum;

type Real = f64; // 計算の精度を決める型

const PARTICLES: u64 = 30_000; // アンサンブル平均に用いる粒子数
const STEPS: u64 = 100_000; // シミュレーションの時間ステップ数
const DELTA_T: Real = 0.000_1; // 時間刻み幅
const TIME: Real = STEPS as Real * DELTA_T; // 総シミュレーション時間
const LENGTH: Real = 0.01; // ディパーティクルの長さ

fn main() {
    let start = std::time::Instant::now();

    let mut mu_writer = BufWriter::new(File::create("data/mu_di.dat").unwrap());
    let mut d_writer = BufWriter::new(File::create("data/d_eff_di.dat").unwrap());
    let mut alpha_writer = BufWriter::new(File::create("data/alpha_di.dat").unwrap());

    for i in 0..=150 {
        let f_x = i as Real;

        let (mean, mean_square) =
            simulate_brownian_motion(STEPS, PARTICLES, DELTA_T, Vector2::new(f_x, 0.0));
        let (mean_rev, mean_square_rev) =
            simulate_brownian_motion(STEPS, PARTICLES, DELTA_T, Vector2::new(-f_x, 0.0));

        let mu = nonlinear_mobility(mean / TIME, f_x);
        let mu_rev = nonlinear_mobility(mean_rev / TIME, -f_x);

        writeln!(mu_writer, "{} {} {}", f_x, mu, mu_rev).unwrap();
        writeln!(
            d_writer,
            "{} {} {}",
            f_x,
            effective_diffusion(mean, mean_square, TIME),
            effective_diffusion(mean_rev, mean_square_rev, TIME)
        )
        .unwrap();
        writeln!(alpha_writer, "{} {}", f_x, alpha(mu, mu_rev)).unwrap();
    }

    println!("Elapsed: {:.2?}", start.elapsed());
}

/// ブラウン運動のシミュレーションを実行し、粒子の平均変位と平均二乗変位を返す
fn simulate_brownian_motion<T>(steps: u64, particles: u64, delta_t: T, f: Vector2<T>) -> (T, T)
where
    T: RealField + SampleUniform + Sum + Copy,
    StandardNormal: Distribution<T>,
{
    let scale = convert::<_, T>(SQRT_2) * delta_t.sqrt();
    let displacements = (0..particles)
        .into_par_iter() // 各粒子のシミュレーションを並列化
        .map(|i| {
            let mut rng = SmallRng::seed_from_u64(i);
            let particle = Diparticle::new(&mut rng, convert(LENGTH));

            (0..steps)
                .map(|_| {
                    // 微小時間に粒子に加わる外力F + ブラウン運動
                    (
                        f * delta_t + noise(&mut rng, scale), // 端1に加わる力
                        f * delta_t + noise(&mut rng, scale), // 端2に加わる力
                    )
                })
                .fold(particle, |mut acc, (dr_1, dr_2)| {
                    acc.apply_forces([dr_1, dr_2]);
                    acc
                })
                .displacement()
        })
        .collect::<Vec<_>>();

    (
        displacements.iter().copied().sum::<T>() / convert(particles as Real),
        displacements.iter().copied().map(|x| x * x).sum::<T>() / convert(particles as Real),
    )
}

/// 正規分布に従う2次元ノイズベクトルを生成
fn noise<T, R>(rng: &mut R, scale: T) -> Vector2<T>
where
    T: RealField + SampleUniform,
    R: rand::Rng + ?Sized,
    StandardNormal: Distribution<T>,
{
    Vector2::new(
        rng.sample::<T, _>(StandardNormal),
        rng.sample::<T, _>(StandardNormal),
    ) * scale
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
