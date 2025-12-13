use crate::{
    Particle, binary_cmp_root, binary_search_root,
    boundary::{
        Boundary, Minus, Plus,
        Zone::{self, *},
        normal, omega,
    },
    reflect,
};
use nalgebra::{Point2, RealField, Scalar, Unit, Vector2, convert};
use rand_distr::uniform::SampleUniform;
use std::cmp::Ordering::*;

mod edge_displacement;

use edge_displacement::{ConstrainedEdgeDisplacements, EdgeDisplacements};

#[derive(Clone, Copy, Debug)]
pub struct Diparticle<T: Scalar> {
    initial: Point2<T>, // 初期位置
    current: Point2<T>, // 現在位置
    length: T,          // 粒子間の長さ
    angle: T,           // 粒子間の相対距離ベクトル(r_1 - r_2)がx軸となす偏角
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
        let initial = random_point::<R, T>(rng, size * convert::<_, T>(0.5));
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
    // 2. 粒子の両端はそれぞれ高々1回しか壁に衝突しない
    fn apply_forces(&mut self, forces: [Vector2<T>; 2]) {
        // 両端の現在位置を求める
        let (p1, p2) = self.edge_positions();
        // 壁との衝突を考えない場合の、両端の変位を求める
        let [dr_1, dr_2] = forces;
        let bound_edge_displacements = self.constrained_edge_displacements(&dr_1, &dr_2);
        let (dr_1, dr_2) = bound_edge_displacements.into_inner();
        // 壁との衝突を考えない場合の、両端の変位を適用した後の位置を求める
        let (tentative_1, tentative_2) = (p1 + dr_1, p2 + dr_2);

        // DEBUG: 両端の距離が粒子の長さと等しいことを確認
        // let dist = (tentative_1 - tentative_2).norm();
        // assert!(
        //     dist.ulps_eq(&self.length, T::default_epsilon(), 8),
        //     "edge distance norm={:?}, expected length={:?}",
        //     convert::<_, T>(dist),
        //     convert::<_, T>(self.length)
        // );

        match (
            Zone::of_point(
                omega::<Minus, T>(tentative_1.x)..=omega::<Plus, T>(tentative_1.x),
                &tentative_1.y,
            ),
            Zone::of_point(
                omega::<Minus, T>(tentative_2.x)..=omega::<Plus, T>(tentative_2.x),
                &tentative_2.y,
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

                        // まず端1が衝突するまでの変位を適用する
                        let bound_edge_displacements =
                            self.advance_to_edge1_collision::<Plus>(&bound_edge_displacements);

                        // 続いて端2が衝突するかどうかを判定
                        let (_, p2) = self.edge_positions();
                        self.move_particle(&if omega::<Plus, T>(p2.x + dr_2.x) < p2.y + dr_2.y {
                            // 端2も衝突する場合
                            self.advance_to_edge2_collision::<Plus>(&bound_edge_displacements)
                        } else {
                            // 端2は衝突しない場合
                            bound_edge_displacements
                        })
                    }
                    Equal => {
                        // 両端が同時に衝突する場合
                        // 棒の制約を無視して両端がそれぞれ独立に反射する場合の変位を計算し、結果の変位を適用する
                        self.move_particle(&self.constrained_edge_displacements(
                            &reflect(omega::<Plus, T>, normal::<Plus, T>, &p1, &dr_1),
                            &reflect(omega::<Plus, T>, normal::<Plus, T>, &p2, &dr_2),
                        ))
                    }
                    Greater => {
                        // 先に端2が衝突する場合

                        // まず端2が衝突するまでの変位を適用する
                        let bound_edge_displacements =
                            self.advance_to_edge2_collision::<Plus>(&bound_edge_displacements);

                        // 続いて端1が衝突するかどうかを判定
                        let (p1, _) = self.edge_positions();
                        self.move_particle(&if omega::<Plus, T>(p1.x + dr_1.x) < p1.y + dr_1.y {
                            // 端1も衝突する場合
                            self.advance_to_edge1_collision::<Plus>(&bound_edge_displacements)
                        } else {
                            // 端1は衝突しない場合
                            bound_edge_displacements
                        });
                    }
                }
            }
            (Above, Within) => {
                // 端1が上壁に衝突する場合

                // まず端1が衝突するまでの変位を適用する
                let bound_edge_displacements =
                    self.advance_to_edge1_collision::<Plus>(&bound_edge_displacements);

                // 続いて端2が衝突するかどうかを判定
                let (_, p2) = self.edge_positions();
                self.move_particle(&if omega::<Plus, T>(p2.x + dr_2.x) < p2.y + dr_2.y {
                    // 端2も衝突する場合
                    self.advance_to_edge2_collision::<Plus>(&bound_edge_displacements)
                } else {
                    // 端2は衝突しない場合
                    bound_edge_displacements
                });
            }
            (Within, Above) => {
                // 端2が上壁に衝突する場合

                // まず端2が衝突するまでの変位を適用する
                let bound_edge_displacements =
                    self.advance_to_edge2_collision::<Plus>(&bound_edge_displacements);

                // 続いて端1が衝突するかどうかを判定
                let (p1, _) = self.edge_positions();
                self.move_particle(&if omega::<Plus, T>(p1.x + dr_1.x) < p1.y + dr_1.y {
                    // 端1も衝突する場合
                    self.advance_to_edge1_collision::<Plus>(&bound_edge_displacements)
                } else {
                    // 端1は衝突しない場合
                    bound_edge_displacements
                });
            }
            (Within, Within) => self.move_particle(&bound_edge_displacements), // 両端とも壁に衝突しない場合
            (Below, Within) => {
                // 端1が下壁に衝突する場合

                // まず端1が衝突するまでの変位を適用する
                let bound_edge_displacements =
                    self.advance_to_edge1_collision::<Minus>(&bound_edge_displacements);

                // 続いて端2が衝突するかどうかを判定
                let (_, p2) = self.edge_positions();
                self.move_particle(&if p2.y + dr_2.y < omega::<Minus, T>(p2.x + dr_2.x) {
                    // 端2も衝突する場合
                    self.advance_to_edge2_collision::<Minus>(&bound_edge_displacements)
                } else {
                    // 端2は衝突しない場合
                    bound_edge_displacements
                });
            }
            (Within, Below) => {
                // 端2が下壁に衝突する場合

                // まず端2が衝突するまでの変位を適用する
                let bound_edge_displacements =
                    self.advance_to_edge2_collision::<Minus>(&bound_edge_displacements);

                // 続いて端1が衝突するかどうかを判定
                let (p1, _) = self.edge_positions();
                self.move_particle(&if p1.y + dr_1.y < omega::<Minus, T>(p1.x + dr_1.x) {
                    // 端1も衝突する場合
                    self.advance_to_edge1_collision::<Minus>(&bound_edge_displacements)
                } else {
                    // 端1は衝突しない場合
                    bound_edge_displacements
                });
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

                        // まず端1が衝突するまでの変位を適用する
                        let bound_edge_displacements =
                            self.advance_to_edge1_collision::<Minus>(&bound_edge_displacements);

                        // 続いて端2が衝突するかどうかを判定
                        let (_, p2) = self.edge_positions();
                        self.move_particle(&if p2.y + dr_2.y < omega::<Minus, T>(p2.x + dr_2.x) {
                            // 端2も衝突する場合
                            self.advance_to_edge2_collision::<Minus>(&bound_edge_displacements)
                        } else {
                            // 端2は衝突しない場合
                            bound_edge_displacements
                        })
                    }
                    Equal => {
                        // 両端が同時に衝突する場合
                        // 棒の制約を無視して両端がそれぞれ独立に反射する場合の変位を計算し、結果の変位を適用する
                        self.move_particle(&self.constrained_edge_displacements(
                            &reflect(omega::<Minus, T>, normal::<Minus, T>, &p1, &dr_1),
                            &reflect(omega::<Minus, T>, normal::<Minus, T>, &p2, &dr_2),
                        ));
                    }
                    Greater => {
                        // 先に端2が衝突する場合

                        // まず端2が衝突するまでの変位を適用する
                        let bound_edge_displacements =
                            self.advance_to_edge2_collision::<Minus>(&bound_edge_displacements);

                        // 続いて端1が衝突するかどうかを判定
                        let (p1, _) = self.edge_positions();
                        self.move_particle(&if p1.y + dr_1.y < omega::<Minus, T>(p1.x + dr_1.x) {
                            // 端1も衝突する場合
                            self.advance_to_edge1_collision::<Minus>(&bound_edge_displacements)
                        } else {
                            // 端1は衝突しない場合
                            bound_edge_displacements
                        });
                    }
                }
            }
            _ => {
                unreachable!("粒子が上下の壁に同時に衝突することはない")
            }
        };
    }
}

impl<T> Diparticle<T>
where
    T: RealField + Copy,
{
    /// 粒子の動径方向についての単位ベクトル
    #[inline]
    fn e_r(&self) -> Unit<Vector2<T>> {
        let (s, c) = self.angle.sin_cos();
        Unit::new_unchecked(Vector2::new(c, s))
    }

    /// 粒子の偏角方向についての単位ベクトル
    #[inline]
    fn e_theta(&self) -> Unit<Vector2<T>> {
        let (s, c) = self.angle.sin_cos();
        Unit::new_unchecked(Vector2::new(-s, c))
    }

    /// 粒子の重心位置と偏角の差分を適用する
    #[inline]
    fn move_particle(&mut self, edge_displacements: &ConstrainedEdgeDisplacements<T>) {
        let (dr_1, dr_2) = edge_displacements.into_inner();
        self.current += (dr_1 + dr_2) * convert::<_, T>(0.5);
        self.angle += (dr_1 - dr_2).norm() / self.length;
    }

    /// 両端の現在位置を求める
    #[inline]
    fn edge_positions(&self) -> (Point2<T>, Point2<T>) {
        let e_r = self.e_r().into_inner();
        (
            self.current + e_r * self.length * convert::<_, T>(0.5),
            self.current - e_r * self.length * convert::<_, T>(0.5),
        )
    }

    /// 両端の仮想的な変位から、実現可能な両端の変位を求める
    ///
    /// 棒の長さが変化しないような制約を課す
    #[inline]
    fn constrained_edge_displacements(
        &self,
        dr_1: &Vector2<T>,
        dr_2: &Vector2<T>,
    ) -> ConstrainedEdgeDisplacements<T> {
        (*dr_1, *dr_2).constrain_along(&self.e_r())
    }

    /// 端1が壁に衝突し反射するまでの変位を適用し、それ以降の端の変位を返す
    fn advance_to_edge1_collision<B>(
        &mut self,
        edge_displacements: &ConstrainedEdgeDisplacements<T>,
    ) -> ConstrainedEdgeDisplacements<T>
    where
        B: Boundary,
    {
        // 端1が衝突するまでの時間を求める
        let (p1, _) = self.edge_positions();
        let (v_1, v_2) = edge_displacements.into_inner();
        let t = binary_search_root(
            |&t| omega::<B, T>(p1.x + v_1.x * t) - (p1.y + v_1.y * t),
            T::zero(),
            T::one(),
        );

        // 端1が衝突するまでの変位を適用
        self.move_particle(&(edge_displacements * t));

        // 端1の残りの移動ベクトルを反射させる
        let v1_rest = v_1 * (T::one() - t);
        let v2_rest = v_2 * (T::one() - t);
        let n = normal::<B, T>(p1.x + v_1.x * t); // 衝突点における法線ベクトル
        let v_1_reflected = v1_rest - n * n.dot(&v1_rest) * convert::<_, T>(2.0);

        // 端1の反射後の変位から、棒の両端の変位を求める
        self.constrained_edge_displacements(&v_1_reflected, &v2_rest)
    }

    // 端2が壁に衝突し反射するまで時間を進める
    fn advance_to_edge2_collision<B>(
        &mut self,
        edge_displacements: &ConstrainedEdgeDisplacements<T>,
    ) -> ConstrainedEdgeDisplacements<T>
    where
        B: Boundary,
    {
        // 端2が衝突するまでの時間を求める
        let (_, p2) = self.edge_positions();
        let (v_1, v_2) = edge_displacements.into_inner();
        let t = binary_search_root(
            |&t| omega::<B, T>(p2.x + v_2.x * t) - (p2.y + v_2.y * t),
            T::zero(),
            T::one(),
        );

        // 端2が衝突するまでの変位を適用
        self.move_particle(&(edge_displacements * t));

        // 端2の移動ベクトルを反射させる
        let v1_rest = v_1 * (T::one() - t);
        let v2_rest = v_2 * (T::one() - t);
        let n = normal::<B, T>(p2.x + v_2.x * t); // 衝突点における法線ベクトル
        let v_2_reflected = v2_rest - n * n.dot(&v2_rest) * convert::<_, T>(2.0);

        // 端2の反射後の変位から、棒の両端の変位を求める
        self.constrained_edge_displacements(&v1_rest, &v_2_reflected)
    }
}

fn random_point<R, T>(rng: &mut R, radius: T) -> Point2<T>
where
    R: rand::Rng + ?Sized,
    T: RealField + SampleUniform + Copy,
{
    let x = rng.random_range(T::zero()..T::one());
    let y = rng.random_range((omega::<Minus, T>(x) + radius)..(omega::<Plus, T>(x) - radius));
    Point2::new(x, y)
}
