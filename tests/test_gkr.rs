// tests/linear_gkr_test.rs

#[macro_use]
extern crate lazy_static;

use ark_bls12_381::Fr as ScalarField;
use rstest::rstest;
use std::collections::HashMap;

// 各モジュールは src 内の実装（lib.rs 経由で公開）を利用する
use gkr::ml_extension::{DenseMLE, SparseMLE};
use gkr::prover::LinearGKRProver;
use gkr::verifier::LinearGKRVerifier;

lazy_static! {
    // f1: 3 変数の疎な multilinear extension（定数1 の回路）を全評価で定義する
    static ref F1: SparseMLE<ScalarField> = {
        let mut evals = HashMap::new();
        // {0,1}^3 の全 8 点に対して値 1
        for i in 0..(1 << 3) {
            evals.insert(i, 1u32.into());
        }
        SparseMLE { num_vars: 3, evaluations: evals }
    };

    // f2: 1 変数の密な multilinear extension、評価：f2(0)=2, f2(1)=3
    static ref F2: DenseMLE<ScalarField> = {
        DenseMLE::from_evaluations_vec(1, vec![2u32.into(), 3u32.into()])
    };

    // f3: 1 変数の密な multilinear extension、評価：f3(0)=4, f3(1)=5
    static ref F3: DenseMLE<ScalarField> = {
        DenseMLE::from_evaluations_vec(1, vec![4u32.into(), 5u32.into()])
    };

    // 固定ベクトル g（長さ 1）
    static ref G: Vec<ScalarField> = vec![1u32.into()];
}

#[rstest]
fn linear_gkr_test() {
    // 乱数生成器を用意（実際の実装では Fiat–Shamir 等の変換も可能）
    let mut rng = rand::thread_rng();

    // Prover 側：Linear GKR プロトコルの証明を生成
    let proof = LinearGKRProver::prove(&F1, &F2, &F3, &G, &mut rng);

    // 上記の各定義から，Phase1 での claimed sum は以下のように計算できる:
    // f1 は定数 1 で，g = [1] により f1(g,x,y) は {0,1}^2 上の定数 1 となる．
    // したがって h_g(x) = sum_{y in {0,1}} f1(g,x,y) * f3(y) = 4 + 5 = 9（x に依存せず一定）．
    // f2 は [2, 3] なので，
    // claimed_sum = 9 * f2(0) + 9 * f2(1) = 9*2 + 9*3 = 18 + 27 = 45.
    let claimed_sum_phase1: ScalarField = 45u32.into();

    // Verifier 側：Prover から受け取った証明を検証する
    let subclaim = LinearGKRVerifier::verify(1, claimed_sum_phase1, &proof, &mut rng);
    assert!(subclaim.is_ok(), "Linear GKR proof verification failed");
}

