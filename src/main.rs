use nalgebra::{RealField, Vector2, convert};
use rand::{Rng, SeedableRng, rngs::SmallRng};
use rand_distr::{Distribution, StandardNormal, uniform::SampleUniform};
use rayon::prelude::*;
use rectification::{Diparticle, Particle};
use std::f64::consts::SQRT_2;
use std::fs::File;
use std::io::{BufWriter, Write};

type Real = f64; // 計算の精度を決める型

// 単粒子の場合のパラメータ(_01)
// const STEPS: usize = 1_000_000; // シミュレーションの時間ステップ数  10^6
// const DELTA_T: Real = 1e-8; //     時間刻み幅                10^-8

// 双粒子の場合のパラメータ
const PARTICLES: u64 = 30_000; //                アンサンブル平均に用いる粒子数  3×10^4
const STEPS: usize = 100_000; //                 シミュレーションの時間ステップ数  10^5
const TIME: Real = STEPS as Real * DELTA_T; //   総シミュレーション時間         0.001 : 短すぎる
const LENGTH: Real = 0.01; //                    ディパーティクルの長さ         0.01 < チャネルの最狭部の幅 約0.038
const DELTA_T: Real = LENGTH * LENGTH * 1e-4; // 時間刻み幅 (√δt = LENGTH/100 となるように設定)

fn main() {
    let start = std::time::Instant::now();

    let mut mu_writer = BufWriter::new(File::create("data/di/mu_150_001.dat").unwrap());
    let mut d_writer = BufWriter::new(File::create("data/di/d_eff_150_001.dat").unwrap());
    let mut alpha_writer = BufWriter::new(File::create("data/di/alpha_150_001.dat").unwrap());

    for i in 1..=150 {
        let f = Vector2::new(i as Real, 0.0);

        let (mean, mean_square) =
            simulate_brownian_motion::<Diparticle<_>, _, _>(STEPS, PARTICLES, DELTA_T, f, LENGTH);
        let (mean_rev, mean_square_rev) =
            simulate_brownian_motion::<Diparticle<_>, _, _>(STEPS, PARTICLES, DELTA_T, -f, LENGTH);

        let mu = nonlinear_mobility(mean / TIME, f.x);
        let mu_rev = nonlinear_mobility(mean_rev / TIME, -f.x);
        writeln!(mu_writer, "{} {} {}", f.x, mu, mu_rev).unwrap();
        writeln!(
            d_writer,
            "{} {} {}",
            f.x,
            effective_diffusion(mean, mean_square, TIME),
            effective_diffusion(mean_rev, mean_square_rev, TIME)
        )
        .unwrap();
        writeln!(alpha_writer, "{} {}", f.x, alpha(mu, mu_rev)).unwrap();
    }

    println!("Elapsed: {:.2?}", start.elapsed());
}

/// ブラウン運動のシミュレーションを実行し、粒子の平均変位と平均二乗変位を返す
#[allow(dead_code)]
fn simulate_brownian_motion<P, T, const C: usize>(
    steps: usize,
    particles: u64,
    delta_t: T,
    f: Vector2<T>,
    length: P::Size,
) -> (T, T)
where
    P: Particle<T, C>,
    T: RealField + SampleUniform + Copy,
    StandardNormal: Distribution<T>,
{
    let scale = convert::<_, T>(SQRT_2) * delta_t.sqrt();

    (0..particles)
        .into_par_iter() // 各粒子のシミュレーションを並列化
        .map(|i| {
            let mut rng = SmallRng::seed_from_u64(i);
            let particle = P::new(&mut rng, length);

            let delta_x = std::iter::repeat_with(|| {
                // 微小時間に粒子に加わる力 = 外力F + ブラウン運動
                std::array::from_fn(|_| f * delta_t + noise(&mut rng, scale))
            })
            .take(steps)
            .fold(particle, |mut p, forces| {
                p.apply_forces(forces);
                p
            })
            .displacement();

            (delta_x, delta_x * delta_x) // (変位, 二乗変位)
        })
        .reduce_with(|(a, aa), (x, xx)| (a + x, aa + xx))
        .map(|(sum, sq_sum)| {
            (
                sum / convert(particles as f64),    // 平均変位
                sq_sum / convert(particles as f64), // 平均二乗変位
            )
        })
        .unwrap()
}

/// 正規分布に従う2次元ノイズベクトルを生成
fn noise<T, R>(rng: &mut R, scale: T) -> Vector2<T>
where
    T: RealField + SampleUniform,
    R: Rng + ?Sized,
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

#[cfg(test)]
mod test {
    use nalgebra::Vector2;
    use rectification::{Diparticle, Monoparticle};

    #[test]
    fn test_monoparticle() {
        let (_, _): (f64, f64) = super::simulate_brownian_motion::<Monoparticle<_>, _, _>(
            100_000,
            30_000,
            1e-4,
            Vector2::new(1.0, 0.0),
            (),
        );
    }

    #[test]
    fn test_diparticle() {
        // 実行時間: 10^11 [step*pariticles] -> 6810秒 (約1.9時間)
        super::simulate_brownian_motion::<Diparticle<_>, _, _>(
            10_000_000, // 10^7
            30_000,     // 3×10^4
            1e-8,
            Vector2::new(1.0, 0.0),
            0.01, // チャネルの最狭部の幅は 約0.038 なので、それより小さい値にする
        );
    }
}
