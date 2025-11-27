use crate::boundary::{
    Boundary, Minus, Plus,
    Zone::{self, *},
    normal, omega, random_point,
};
use crate::{Particle, binary_cmp_root, binary_search_root, reflect};
use nalgebra::{Point2, RealField, Scalar, Vector2, convert};
use rand_distr::uniform::SampleUniform;
use std::cmp::Ordering::*;

#[derive(Clone, Copy, Debug)]
pub struct Diparticle<T: Scalar> {
    initial: Point2<T>,
    current: Point2<T>,
    length: T,
    angle: T,
}

impl<T> Particle<T, 2> for Diparticle<T>
where
    T: RealField + SampleUniform + Copy,
{
    type Size = T;

    fn new<R>(rng: &mut R, size: T) -> Self
    where
        R: rand::Rng,
    {
        let initial = random_point::<R, T>(rng);
        Self {
            initial,
            current: initial,
            length: size,
            angle: rng.random_range(T::zero()..T::two_pi()),
        }
    }

    #[inline]
    fn displacement(&self) -> T {
        self.current.x - self.initial.x
    }

    // この関数の挙動には次のような仮定が置かれている
    // 1. 粒子の両端が同時に相異なる壁に衝突することはない
    // 2. 粒子の両端は高々1回しか壁に衝突しない
    fn apply_forces(&mut self, forces: [Vector2<T>; 2]) {
        // 両端にある粒子の現在位置を求める
        let (p1, p2) = self.edge_positions();
        // 壁との衝突を考えない場合の両端の変位を求める
        let [dr_1, dr_2] = forces;
        let (dr_1, dr_2) = self.edges_displacements(&dr_1, &dr_2);
        // 壁との衝突を考えない場合の、両端の変位を適用した後の位置を求める
        let (tentative_1, tentative_2) = (p1 + dr_1, p2 + dr_2);

        let (dr_c, dtheta) = match (
            Zone::of_point(
                omega::<Minus, T>(tentative_1.x)..=omega::<Plus, T>(tentative_1.x),
                &tentative_1,
            ),
            Zone::of_point(
                omega::<Minus, T>(tentative_2.x)..=omega::<Plus, T>(tentative_2.x),
                &tentative_2,
            ),
        ) {
            (Above, Above) => {
                // 両端がどちらも上壁に衝突する場合
                // 両端のうちどちらが先に衝突するかを判定
                match binary_cmp_root(
                    |&t| omega::<Plus, T>(p1.x + dr_1.x * t) - (p1.y + dr_1.y * t), // 端1の衝突時間
                    |&t| omega::<Plus, T>(p2.x + dr_2.x * t) - (p2.y + dr_2.y * t), // 端2の衝突時間
                ) {
                    Less => {
                        // 先に端1が衝突する場合
                        let (dr_1, dr_2) = self.advance_to_p1_collision::<Plus>(&dr_1, &dr_2);
                        let (_, p2) = self.edge_positions();

                        // 端2が衝突するかどうかを判定
                        if omega::<Plus, T>(p2.x + dr_2.x) < p2.y + dr_2.y {
                            // 端2も衝突する場合
                            let (dr_1, dr_2) = self.advance_to_p2_collision::<Plus>(&dr_1, &dr_2);
                            self.delta(&dr_1, &dr_2)
                        } else {
                            // 端2は衝突しない場合
                            self.delta(&dr_1, &dr_2)
                        }
                    }
                    Equal => {
                        // 両端が同時に衝突する場合
                        let (p1, p2) = self.edge_positions();
                        self.delta(
                            &reflect(omega::<Plus, T>, normal::<Plus, T>, &p1, &dr_1),
                            &reflect(omega::<Plus, T>, normal::<Plus, T>, &p2, &dr_2),
                        )
                    }
                    Greater => {
                        // 先に端2が衝突する場合
                        let (dr_1, dr_2) = self.advance_to_p2_collision::<Plus>(&dr_1, &dr_2);
                        let (p1, _) = self.edge_positions();

                        // 端1が衝突するかどうかを判定
                        if omega::<Plus, T>(p1.x + dr_1.x) < p1.y + dr_1.y {
                            // 端1も衝突する場合
                            let (dr_1, dr_2) = self.advance_to_p1_collision::<Plus>(&dr_1, &dr_2);
                            self.delta(&dr_1, &dr_2)
                        } else {
                            // 端1は衝突しない場合
                            self.delta(&dr_1, &dr_2)
                        }
                    }
                }
            }
            (Above, Within) => {
                // 端1が上壁に衝突する場合
                let (dr_1, dr_2) = self.advance_to_p1_collision::<Plus>(&dr_1, &dr_2);
                let (_, p2) = self.edge_positions();

                // 端2が衝突するかどうかを判定
                if omega::<Plus, T>(p2.x + dr_2.x) < p2.y + dr_2.y {
                    // 端2も衝突する場合
                    let (dr_1, dr_2) = self.advance_to_p2_collision::<Plus>(&dr_1, &dr_2);
                    self.delta(&dr_1, &dr_2)
                } else {
                    // 端2は衝突しない場合
                    self.delta(&dr_1, &dr_2)
                }
            }
            (Within, Above) => {
                // 端2が上壁に衝突する場合
                let (dr_1, dr_2) = self.advance_to_p2_collision::<Plus>(&dr_1, &dr_2);
                let (p1, _) = self.edge_positions();

                // 端1が衝突するかどうかを判定
                if omega::<Plus, T>(p1.x + dr_1.x) < p1.y + dr_1.y {
                    // 端1も衝突する場合
                    let (dr_1, dr_2) = self.advance_to_p1_collision::<Plus>(&dr_1, &dr_2);
                    self.delta(&dr_1, &dr_2)
                } else {
                    // 端1は衝突しない場合
                    self.delta(&dr_1, &dr_2)
                }
            }
            (Within, Within) => self.delta(&dr_1, &dr_2), // 両端とも壁に衝突しない場合
            (Below, Within) => {
                // 端1が下壁に衝突する場合
                let (dr_1, dr_2) = self.advance_to_p1_collision::<Minus>(&dr_1, &dr_2);
                let (_, p2) = self.edge_positions();

                // 端2が衝突するかどうかを判定
                if p2.y + dr_2.y < omega::<Minus, T>(p2.x + dr_2.x) {
                    // 端2も衝突する場合
                    let (dr_1, dr_2) = self.advance_to_p2_collision::<Minus>(&dr_1, &dr_2);
                    self.delta(&dr_1, &dr_2)
                } else {
                    // 端2は衝突しない場合
                    self.delta(&dr_1, &dr_2)
                }
            }
            (Within, Below) => {
                // 端2が下壁に衝突する場合
                let (dr_1, dr_2) = self.advance_to_p2_collision::<Minus>(&dr_1, &dr_2);
                let (p1, _) = self.edge_positions();

                // 端1が衝突するかどうかを判定
                if p1.y + dr_1.y < omega::<Minus, T>(p1.x + dr_1.x) {
                    // 端1も衝突する場合
                    let (dr_1, dr_2) = self.advance_to_p1_collision::<Minus>(&dr_1, &dr_2);
                    self.delta(&dr_1, &dr_2)
                } else {
                    // 端1は衝突しない場合
                    self.delta(&dr_1, &dr_2)
                }
            }
            (Below, Below) => {
                // 両端がどちらも下壁に衝突する場合
                // 両端のうちどちらが先に衝突するかを判定
                match binary_cmp_root(
                    |&t| omega::<Minus, T>(p1.x + dr_1.x * t) - (p1.y + dr_1.y * t), // 端1の衝突時間
                    |&t| omega::<Minus, T>(p2.x + dr_2.x * t) - (p2.y + dr_2.y * t), // 端2の衝突時間
                ) {
                    Less => {
                        // 先に端1が衝突する場合
                        let (dr_1, dr_2) = self.advance_to_p1_collision::<Minus>(&dr_1, &dr_2);
                        let (_, p2) = self.edge_positions();

                        // 端2が衝突するかどうかを判定
                        if p2.y + dr_2.y < omega::<Minus, T>(p2.x + dr_2.x) {
                            // 端2も衝突する場合
                            let (dr_1, dr_2) = self.advance_to_p2_collision::<Minus>(&dr_1, &dr_2);
                            self.delta(&dr_1, &dr_2)
                        } else {
                            // 端2は衝突しない場合
                            self.delta(&dr_1, &dr_2)
                        }
                    }
                    Equal => {
                        // 両端が同時に衝突する場合
                        let (p1, p2) = self.edge_positions();
                        self.delta(
                            &reflect(omega::<Minus, T>, normal::<Minus, T>, &p1, &dr_1),
                            &reflect(omega::<Minus, T>, normal::<Minus, T>, &p2, &dr_2),
                        )
                    }
                    Greater => {
                        // 先に端2が衝突する場合
                        let (dr_1, dr_2) = self.advance_to_p2_collision::<Minus>(&dr_1, &dr_2);
                        let (p1, _) = self.edge_positions();

                        // 端1が衝突するかどうかを判定
                        if p1.y + dr_1.y < omega::<Minus, T>(p1.x + dr_1.x) {
                            // 端1も衝突する場合
                            let (dr_1, dr_2) = self.advance_to_p1_collision::<Minus>(&dr_1, &dr_2);
                            self.delta(&dr_1, &dr_2)
                        } else {
                            // 端1は衝突しない場合
                            self.delta(&dr_1, &dr_2)
                        }
                    }
                }
            }
            _ => {
                unreachable!("粒子が上下の壁に同時に衝突することはない")
            }
        };

        self.current += dr_c;
        self.angle += dtheta;
    }
}

impl<T> Diparticle<T>
where
    T: RealField + Copy,
{
    /// 粒子の動径方向についての単位ベクトル
    #[inline]
    fn e_r(&self) -> Vector2<T> {
        let (s, c) = self.angle.sin_cos();
        Vector2::new(c, s)
    }

    /// 粒子の偏角方向についての単位ベクトル
    #[inline]
    fn e_theta(&self) -> Vector2<T> {
        let (s, c) = self.angle.sin_cos();
        Vector2::new(-s, c)
    }

    /// 両端の変位から、重心と角度の変化量を求める
    #[inline]
    fn delta(&self, dr_1: &Vector2<T>, dr_2: &Vector2<T>) -> (Vector2<T>, T) {
        (
            (dr_1 + dr_2) * convert::<_, T>(0.5),
            (dr_1 - dr_2).dot(&self.e_theta()) / self.length,
        )
    }

    /// 両端の現在位置を求める
    #[inline]
    fn edge_positions(&self) -> (Point2<T>, Point2<T>) {
        let e_r = self.e_r();
        (
            self.current + e_r * self.length * convert::<_, T>(0.5),
            self.current - e_r * self.length * convert::<_, T>(0.5),
        )
    }

    /// 両端の仮想的な変位から、実現可能な両端の変位を求める
    #[inline]
    fn edges_displacements(
        &self,
        dr_1: &Vector2<T>,
        dr_2: &Vector2<T>,
    ) -> (Vector2<T>, Vector2<T>) {
        let dr_c = (dr_1 + dr_2) * convert::<_, T>(0.5);
        let e_theta = self.e_theta();
        let dv = e_theta * e_theta.dot(&(dr_1 - dr_2)) * convert::<_, T>(0.5);
        (dr_c + dv, dr_c - dv)
    }

    /// 端1が壁に衝突し反射するまで時間を進める
    fn advance_to_p1_collision<B>(
        &mut self,
        v_1: &Vector2<T>,
        v_2: &Vector2<T>,
    ) -> (Vector2<T>, Vector2<T>)
    where
        B: Boundary,
    {
        // 端1が衝突するまでの時間を求める
        let (p1, _) = self.edge_positions();
        let t = binary_search_root(
            |&t| omega::<B, T>(p1.x + v_1.x * t) - (p1.y + v_1.y * t),
            T::zero(),
            T::one(),
        );

        // 端1が衝突するまでの変位を適用
        let (dr_c, dtheta) = self.delta(&(v_1 * t), &(v_2 * t));
        self.current += dr_c;
        self.angle += dtheta;

        // 端1の移動ベクトルを反射させる
        let v_rest = v_1 * (T::one() - t);
        let n = normal::<B, T>(p1.x + v_1.x * t);
        let v_1_reflected = v_rest - n * n.dot(&v_rest) * convert::<_, T>(2.0);

        // 反射後の変位から実際に実現される両端の変位を求める
        self.edges_displacements(&v_1_reflected, v_2)
    }

    // 端2が壁に衝突し反射するまで時間を進める
    fn advance_to_p2_collision<B>(
        &mut self,
        v_1: &Vector2<T>,
        v_2: &Vector2<T>,
    ) -> (Vector2<T>, Vector2<T>)
    where
        B: Boundary,
    {
        // 端2が衝突するまでの時間を求める
        let (_, p2) = self.edge_positions();
        let t = binary_search_root(
            |&t| omega::<B, T>(p2.x + v_2.x * t) - (p2.y + v_2.y * t),
            T::zero(),
            T::one(),
        );

        // 端2が衝突するまでの変位を適用
        let (dr_c, dtheta) = self.delta(&(v_1 * t), &(v_2 * t));
        self.current += dr_c;
        self.angle += dtheta;

        // 端2の移動ベクトルを反射させる
        let v_rest = v_2 * (T::one() - t);
        let n = normal::<B, T>(p2.x + v_2.x * t);
        let v_2_reflected = v_rest - n * n.dot(&v_rest) * convert::<_, T>(2.0);

        // 反射後の変位から実際に実現される両端の変位を求める
        self.edges_displacements(v_1, &v_2_reflected)
    }
}
