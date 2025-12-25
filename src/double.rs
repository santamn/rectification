use crate::{
    Particle, RangeExt,
    RelativePos::*,
    binary_cmp_root, binary_search_root,
    boundary::{Bottom, Boundary, Top},
};
use nalgebra::{Point2, RealField, Scalar, Unit, Vector2, convert};
use rand_distr::uniform::SampleUniform;
use std::cmp::Ordering::*;

mod edge_displacement;

use edge_displacement::{ConstrainedEdgeDisplacements, EdgeDisplacements};

#[derive(Clone, Debug)]
pub struct Diparticle<T: Scalar> {
    initial: Point2<T>,           // 初期位置
    current: Point2<T>,           // 現在位置
    length: T,                    // 粒子間の長さ
    angle: T,                     // 棒がx軸となす偏角
    history: Vec<(Point2<T>, T)>, // 履歴
}

impl<T> Particle<T, 2> for Diparticle<T>
where
    T: RealField + SampleUniform + Copy + Send + Sync,
{
    type Size = T;

    fn new<R>(rng: &mut R, size: T) -> Self
    where
        R: rand::Rng,
    {
        let initial_point = random_point::<R, T>(rng, size * convert::<_, T>(0.5));
        let initial_angle = rng.random_range(T::zero()..T::two_pi());
        let mut history = Vec::with_capacity(3_000_000);
        history.push((initial_point, initial_angle));
        Self {
            initial: initial_point,
            current: initial_point,
            length: size,
            angle: initial_angle,
            history,
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
        let bound_edge_displacements = self.constrain_edge_displacements(&dr_1, &dr_2);
        let (dr_1, dr_2) = bound_edge_displacements.into_inner();
        // 壁との衝突を考えない場合の、両端の変位を適用した後の位置を求める
        let (tentative_1, tentative_2) = (p1 + dr_1, p2 + dr_2);

        match (
            (Bottom::<T>::f(&tentative_1.x)..=Top::<T>::f(&tentative_1.x))
                .position_of(&tentative_1.y),
            (Bottom::<T>::f(&tentative_2.x)..=Top::<T>::f(&tentative_2.x))
                .position_of(&tentative_2.y),
        ) {
            // 両端がどちらも上壁に衝突する場合
            (Above, Above) => {
                // 両端のうちどちらが先に衝突するかを判定
                match binary_cmp_root(
                    |&t| Top::<T>::f(&(p1.x + dr_1.x * t)) - (p1.y + dr_1.y * t), // 端1の衝突時間
                    |&t| Top::<T>::f(&(p2.x + dr_2.x * t)) - (p2.y + dr_2.y * t), // 端2の衝突時間
                ) {
                    // 先に端1が衝突する場合
                    Less => {
                        // まず端1が衝突するまでの変位を適用する
                        let dr_after_e1_collided =
                            self.advance_to_edge1_collision::<Top<T>>(&bound_edge_displacements);

                        // 端2の衝突を処理
                        let (_, p2) = self.edge_positions();
                        let (_, dr_2) = dr_after_e1_collided.into_inner();
                        let dr = if Top::<T>::is_outside(&(p2 + dr_2)) {
                            // 端2も衝突する場合
                            self.advance_to_edge2_collision::<Top<T>>(&dr_after_e1_collided)
                        } else {
                            // 端2は衝突しない場合
                            dr_after_e1_collided
                        };
                        self.move_particle(&dr);
                    }
                    // 両端が同時に衝突する場合
                    Equal => {
                        // 棒の制約を無視して両端がそれぞれ独立に反射するとして両端の変位を計算し、結果の変位を適用する
                        let time_to_hit = binary_search_root(
                            |&t| Top::<T>::f(&(p1.x + dr_1.x * t)) - (p1.y + dr_1.y * t),
                            &T::zero(),
                            &T::one(),
                        );

                        let to_boundary = bound_edge_displacements * time_to_hit;
                        let (dr_1_rest, dr_2_rest) =
                            (bound_edge_displacements - to_boundary).into_inner();
                        self.move_particle(
                            &(to_boundary
                                + self.constrain_edge_displacements(
                                    &Top::<T>::reflect_at(&(p1 + dr_1 * time_to_hit).x, &dr_1_rest),
                                    &Top::<T>::reflect_at(&(p2 + dr_2 * time_to_hit).x, &dr_2_rest),
                                )),
                        );
                    }
                    // 先に端2が衝突する場合
                    Greater => {
                        // まず端2が衝突するまでの変位を適用する
                        let dr_after_e2_collided =
                            self.advance_to_edge2_collision::<Top<T>>(&bound_edge_displacements);

                        // 端1の衝突を処理
                        let (p1, _) = self.edge_positions();
                        let (dr_1, _) = dr_after_e2_collided.into_inner();
                        let dr = if Top::<T>::is_outside(&(p1 + dr_1)) {
                            // 端1も衝突する場合
                            self.advance_to_edge1_collision::<Top<T>>(&dr_after_e2_collided)
                        } else {
                            // 端1は衝突しない場合
                            dr_after_e2_collided
                        };
                        self.move_particle(&dr);
                    }
                }
            }
            // 端1が上壁に衝突する場合
            (Above, Within) => {
                // まず端1が衝突するまでの変位を適用する
                let dr_after_e1_collided =
                    self.advance_to_edge1_collision::<Top<T>>(&bound_edge_displacements);

                // 端2の衝突を処理
                let (_, p2) = self.edge_positions();
                let (_, dr_2) = dr_after_e1_collided.into_inner();
                let dr = if Top::<T>::is_outside(&(p2 + dr_2)) {
                    // 端2も衝突する場合
                    self.advance_to_edge2_collision::<Top<T>>(&dr_after_e1_collided)
                } else {
                    // 端2は衝突しない場合
                    dr_after_e1_collided
                };
                self.move_particle(&dr);
            }
            // 端2が上壁に衝突する場合
            (Within, Above) => {
                // まず端2が衝突するまでの変位を適用する
                let dr_after_e2_collided =
                    self.advance_to_edge2_collision::<Top<T>>(&bound_edge_displacements);

                // 端1の衝突を処理
                let (p1, _) = self.edge_positions();
                let (dr_1, _) = dr_after_e2_collided.into_inner();
                let dr = if Top::<T>::is_outside(&(p1 + dr_1)) {
                    // 端1も衝突する場合
                    self.advance_to_edge1_collision::<Top<T>>(&dr_after_e2_collided)
                } else {
                    // 端1は衝突しない場合
                    dr_after_e2_collided
                };
                self.move_particle(&dr);
            }
            // 両端とも壁に衝突しない場合
            (Within, Within) => self.move_particle(&bound_edge_displacements),
            // 端1が下壁に衝突する場合
            (Below, Within) => {
                // まず端1が衝突するまでの変位を適用する
                let dr_after_e1_collided =
                    self.advance_to_edge1_collision::<Bottom<T>>(&bound_edge_displacements);

                // 端2の衝突を処理
                let (_, p2) = self.edge_positions();
                let (_, dr_2) = dr_after_e1_collided.into_inner();
                let dr = if Bottom::<T>::is_outside(&(p2 + dr_2)) {
                    // 端2も衝突する場合
                    self.advance_to_edge2_collision::<Bottom<T>>(&dr_after_e1_collided)
                } else {
                    // 端2は衝突しない場合
                    dr_after_e1_collided
                };
                self.move_particle(&dr);
            }
            // 端2が下壁に衝突する場合
            (Within, Below) => {
                // まず端2が衝突するまでの変位を適用する
                let dr_after_e2_collided =
                    self.advance_to_edge2_collision::<Bottom<T>>(&bound_edge_displacements);

                // 端1の衝突を処理
                let (p1, _) = self.edge_positions();
                let (dr_1, _) = dr_after_e2_collided.into_inner();
                let dr = if Bottom::<T>::is_outside(&(p1 + dr_1)) {
                    // 端1も衝突する場合
                    self.advance_to_edge1_collision::<Bottom<T>>(&dr_after_e2_collided)
                } else {
                    // 端1は衝突しない場合
                    dr_after_e2_collided
                };
                self.move_particle(&dr);
            }
            // 両端がどちらも下壁に衝突する場合
            (Below, Below) => {
                // 両端のうちどちらが先に衝突するかを判定
                match binary_cmp_root(
                    |&t| Bottom::<T>::f(&(p1.x + dr_1.x * t)) - (p1.y + dr_1.y * t), // 端1の衝突時間
                    |&t| Bottom::<T>::f(&(p2.x + dr_2.x * t)) - (p2.y + dr_2.y * t), // 端2の衝突時間
                ) {
                    // 先に端1が衝突する場合
                    Less => {
                        // まず端1が衝突するまでの変位を適用する
                        let dr_after_e1_collided =
                            self.advance_to_edge1_collision::<Bottom<T>>(&bound_edge_displacements);

                        // 端2の衝突を処理
                        let (_, p2) = self.edge_positions();
                        let (_, dr_2) = dr_after_e1_collided.into_inner();
                        let dr = if Bottom::<T>::is_outside(&(p2 + dr_2)) {
                            // 端2も衝突する場合
                            self.advance_to_edge2_collision::<Bottom<T>>(&dr_after_e1_collided)
                        } else {
                            // 端2は衝突しない場合
                            dr_after_e1_collided
                        };
                        self.move_particle(&dr);
                    }
                    // 両端が同時に衝突する場合
                    Equal => {
                        // 棒の制約を無視して両端がそれぞれ独立に反射する場合の変位を計算し、結果の変位を適用する
                        let time_to_hit = binary_search_root(
                            |&t| Bottom::<T>::f(&(p1.x + dr_1.x * t)) - (p1.y + dr_1.y * t),
                            &T::zero(),
                            &T::one(),
                        );

                        let to_boundary = bound_edge_displacements * time_to_hit;
                        let (dr_1_rest, dr_2_rest) =
                            (bound_edge_displacements - to_boundary).into_inner();
                        self.move_particle(
                            &(to_boundary
                                + self.constrain_edge_displacements(
                                    &Bottom::<T>::reflect_at(
                                        &(p1 + dr_1 * time_to_hit).x,
                                        &dr_1_rest,
                                    ),
                                    &Bottom::<T>::reflect_at(
                                        &(p2 + dr_2 * time_to_hit).x,
                                        &dr_2_rest,
                                    ),
                                )),
                        );
                    }
                    // 先に端2が衝突する場合
                    Greater => {
                        // まず端2が衝突するまでの変位を適用する
                        let dr_after_e2_collided =
                            self.advance_to_edge2_collision::<Bottom<T>>(&bound_edge_displacements);

                        // 端1の衝突を処理
                        let (p1, _) = self.edge_positions();
                        let (dr_1, _) = dr_after_e2_collided.into_inner();
                        let dr = if Bottom::<T>::is_outside(&(p1 + dr_1)) {
                            // 端1も衝突する場合
                            self.advance_to_edge1_collision::<Bottom<T>>(&dr_after_e2_collided)
                        } else {
                            // 端1は衝突しない場合
                            dr_after_e2_collided
                        };
                        self.move_particle(&dr);
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
    T: RealField + Copy + Send + Sync,
{
    /// 粒子の動径方向についての単位ベクトル
    #[inline]
    fn e_r(&self) -> Unit<Vector2<T>> {
        let (s, c) = self.angle.sin_cos();
        Unit::new_unchecked(Vector2::new(c, s))
    }

    /// 粒子の重心位置と偏角の差分を適用する
    #[inline]
    fn move_particle(&mut self, edge_displacements: &ConstrainedEdgeDisplacements<T>) {
        let (dr_1, dr_2) = edge_displacements.into_inner();
        self.current += (dr_1 + dr_2) * convert::<_, T>(0.5);
        self.angle += (dr_1 - dr_2).norm() / self.length;
        self.history.push((self.current, self.angle));
    }

    #[inline]
    fn history(&self) -> &Vec<(Point2<T>, T)> {
        &self.history
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

    /// 両端の変位について棒の長さが変化しないような制約を課す
    #[inline]
    fn constrain_edge_displacements(
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
        B: Boundary<T>,
    {
        // 端1が衝突するまでの時間を求める
        let (p1, _) = self.edge_positions();
        let (v_1, _) = edge_displacements.into_inner();
        let time_to_hit = binary_search_root(
            |&t| B::f(&(p1.x + v_1.x * t)) - (p1.y + v_1.y * t),
            &T::zero(),
            &T::one(),
        );

        // 端1が衝突するまでの変位を適用
        let to_boundary = edge_displacements * time_to_hit;
        self.move_particle(&to_boundary);

        // 端1の反射後の変位から、棒の両端の変位を求める
        let (e1_to_boundary, _) = to_boundary.into_inner();
        let (v_1_rest, v_2_rest) = (*edge_displacements - to_boundary).into_inner();
        self.constrain_edge_displacements(
            &B::reflect_at(&(p1 + e1_to_boundary).x, &v_1_rest),
            &v_2_rest,
        )
    }

    // 端2が壁に衝突し反射するまで時間を進める
    fn advance_to_edge2_collision<B>(
        &mut self,
        edge_displacements: &ConstrainedEdgeDisplacements<T>,
    ) -> ConstrainedEdgeDisplacements<T>
    where
        B: Boundary<T>,
    {
        // 端2が衝突するまでの時間を求める
        let (_, p2) = self.edge_positions();
        let (_, v_2) = edge_displacements.into_inner();
        let time_to_hit = binary_search_root(
            |&t| B::f(&(p2.x + v_2.x * t)) - (p2.y + v_2.y * t),
            &T::zero(),
            &T::one(),
        );

        // 端2が衝突するまでの変位を適用
        let to_boundary = edge_displacements * time_to_hit;
        self.move_particle(&to_boundary);

        // 端2の移動ベクトルを反射させる
        let (_, e2_to_boundary) = to_boundary.into_inner();
        let (v_1_rest, v_2_rest) = (*edge_displacements - to_boundary).into_inner();
        self.constrain_edge_displacements(
            &v_1_rest,
            &B::reflect_at(&(p2 + e2_to_boundary).x, &v_2_rest),
        )
    }
}

fn random_point<R, T>(rng: &mut R, radius: T) -> Point2<T>
where
    R: rand::Rng + ?Sized,
    T: RealField + SampleUniform + Copy + Send + Sync,
{
    let x = rng.random_range(T::zero()..T::one());
    let y = rng.random_range((Bottom::<T>::f(&x) + radius)..(Top::<T>::f(&x) - radius));
    Point2::new(x, y)
}
