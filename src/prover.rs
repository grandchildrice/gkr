// src/prover.rs

use ark_bls12_381::Fr as ScalarField;
use ark_ff::Zero;
use rand::Rng;
use crate::ml_extension::{DenseMLE, SparseMLE};
use crate::sumcheck::protocol;
use crate::sumcheck::get_r;

/// Linear GKR の証明メッセージ（フェーズごとに Prover から送られるメッセージ列）
pub struct LinearGKRProof {
    pub phase1_msgs: Vec<Vec<ScalarField>>,
    pub phase2_msgs: Vec<Vec<ScalarField>>,
}

/// Linear GKR Prover（型は固定して ScalarField を利用）
pub struct LinearGKRProver;

impl LinearGKRProver {
    /// f1: 3*l 変数の疎な multilinear extension
    /// f2, f3: それぞれ l 変数の密な multilinear extension
    /// g: 固定ベクトル（長さ l）
    pub fn prove<R: Rng>(
        f1: &SparseMLE<ScalarField>,
        f2: &DenseMLE<ScalarField>,
        f3: &DenseMLE<ScalarField>,
        g: &[ScalarField],
        rng: &mut R,
    ) -> LinearGKRProof {
        let l = g.len();

        // ── Phase 1 ──
        // f1 の最初の l 変数を固定し、 h_g(x) = ∑_y f1(g, x, y) * f3(y) を計算
        let (h_g, f1_fixed_g) = initialize_phase_one(f1, f3, g);
        // P1(x) = h_g(x) * f2(x) の全和（sum-check の対象値）を計算
        let claimed_sum_phase1 = compute_claimed_sum(&h_g, f2);
        let mut prover_state1 = protocol::prover_init(l, claimed_sum_phase1);
        let mut phase1_msgs = Vec::with_capacity(l);
        let mut u = Vec::with_capacity(l);

        for _ in 0..l {
            let msg = protocol::prove_round(&mut prover_state1, rng);
            phase1_msgs.push(msg.clone());
            let r_i: ScalarField = get_r().unwrap();
            u.push(r_i);
            protocol::apply_challenge(&mut prover_state1, r_i);
        }

        // ── Phase 2 ──
        // f1_fixed_g は f1(g, x, y) となっているので，さらに x = u を固定して f1(g, u, y) を得る
        let f1_fixed_gu = initialize_phase_two(&f1_fixed_g, &u);
        let f2_at_u = f2.evaluate(&u);
        // Phase2 の対象は P2(y) = f1(g,u,y) * f3(y) * f2(u) と考える
        let claimed_sum_phase2 = f2_at_u * compute_dense_sum(&f1_fixed_gu, f3);
        let mut prover_state2 = protocol::prover_init(l, claimed_sum_phase2);
        let mut phase2_msgs = Vec::with_capacity(l);
        let mut v = Vec::with_capacity(l);

        for _ in 0..l {
            let msg = protocol::prove_round(&mut prover_state2, rng);
            phase2_msgs.push(msg.clone());
            let r_j: ScalarField = get_r().unwrap();
            v.push(r_j);
            protocol::apply_challenge(&mut prover_state2, r_j);
        }

        LinearGKRProof { phase1_msgs, phase2_msgs }
    }
}

/// f1 の最初の l 変数を固定し、 h_g(x) = ∑_y f1(g,x,y)*f3(y) を計算する
fn initialize_phase_one(
    f1: &SparseMLE<ScalarField>,
    f3: &DenseMLE<ScalarField>,
    g: &[ScalarField],
) -> (DenseMLE<ScalarField>, SparseMLE<ScalarField>) {
    let l = g.len();
    assert_eq!(f1.num_vars, 3 * l);
    let f1_fixed_g = f1.fix_variables(g);
    let size = 1 << l;
    let mut h_evals = vec![ScalarField::zero(); size];
    // f1_fixed_g は 2*l 変数（前半 l が x，後半 l が y）として格納されているとする
    for (&index, &val) in f1_fixed_g.evaluations.iter() {
        if val.is_zero() {
            continue;
        }
        let x_index = index & ((1 << l) - 1);
        let y_index = index >> l;
        h_evals[x_index] += val * f3.evaluations[y_index];
    }
    let h_g = DenseMLE::from_evaluations_vec(l, h_evals);
    (h_g, f1_fixed_g)
}

/// Phase1 で固定した f1 の残りの変数を、u（Phase1 の乱数列）で固定して f1(g,u,y) を得る
fn initialize_phase_two(
    f1_fixed_g: &SparseMLE<ScalarField>,
    u: &[ScalarField],
) -> DenseMLE<ScalarField> {
    let f1_fixed_gu = f1_fixed_g.fix_variables(u);
    f1_fixed_gu.to_dense_multilinear_extension()
}

/// Phase1 の claimed sum の計算：∑_x h_g(x)*f2(x)
fn compute_claimed_sum(h_g: &DenseMLE<ScalarField>, f2: &DenseMLE<ScalarField>) -> ScalarField {
    let l = h_g.num_vars;
    let size = 1 << l;
    let mut sum = ScalarField::zero();
    for i in 0..size {
        sum += h_g.evaluations[i] * f2.evaluations[i];
    }
    sum
}

/// Phase2 の dense sum：∑_y f1(g,u,y)*f3(y)
fn compute_dense_sum(f1_fixed_gu: &DenseMLE<ScalarField>, f3: &DenseMLE<ScalarField>) -> ScalarField {
    let l = f1_fixed_gu.num_vars;
    let size = 1 << l;
    let mut sum = ScalarField::zero();
    for i in 0..size {
        sum += f1_fixed_gu.evaluations[i] * f3.evaluations[i];
    }
    sum
}
