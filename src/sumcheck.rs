// src/sumcheck.rs

use ark_bls12_381::Fr as ScalarField;
use ark_ff::Field;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::polynomial::Polynomial;
use ark_poly::DenseMVPolynomial;
// cfg_into_iter! は単純な iter() に置換
use rand::Rng;

/// Sumcheck 用の多変数多項式の型
pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;
pub type UniPoly = UniSparsePolynomial<ScalarField>;

/// i を {0,1}^v 上のインデックスに変換する補助関数
pub fn n_to_vec(i: usize, n: usize) -> Vec<ScalarField> {
    format!("{:0>width$}", format!("{:b}", i), width = n)
        .chars()
        .map(|x| if x == '1' { 1.into() } else { 0.into() })
        .collect()
}

/// 単一の prover インスタンスの「メモリ」を模擬する構造体
#[derive(Debug, Clone)]
pub struct Prover {
    pub g: MultiPoly,
    pub r_vec: Vec<ScalarField>,
}

impl Prover {
    pub fn new(g: &MultiPoly) -> Self {
        Prover {
            g: g.clone(),
            r_vec: vec![],
        }
    }

    // 多項式 g に対して、Xj を固定し xj+1 上で評価した結果（1変数多項式）を生成
    pub fn gen_uni_polynomial(&mut self, r: Option<ScalarField>) -> UniPoly {
        if let Some(r_val) = r {
            self.r_vec.push(r_val);
        }
        let v = self.g.num_vars() - self.r_vec.len();
        (0..(2u32.pow(v as u32 - 1)))
            .fold(UniPoly::from_coefficients_vec(vec![(0, 0u32.into())]),
                  |sum, n| sum + self.evaluate_gj(n_to_vec(n as usize, v)))
    }

    // gj を点列に対して評価し、全ての項を 1 変数多項式にまとめる
    pub fn evaluate_gj(&self, points: Vec<ScalarField>) -> UniPoly {
        self.g.terms().iter().fold(
            UniPoly::from_coefficients_vec(vec![]),
            |sum, (coeff, term)| {
                let (coeff_eval, fixed_term) = self.evaluate_term(&term, &points);
                let curr = match fixed_term {
                    None => UniPoly::from_coefficients_vec(vec![(0, *coeff * coeff_eval)]),
                    Some(ft) => UniPoly::from_coefficients_vec(vec![(ft.degree(), *coeff * coeff_eval)]),
                };
                curr + sum
            },
        )
    }

    // 項 term を固定した場合の評価：(新しい係数, 固定後の項) を返す
    pub fn evaluate_term(
        &self,
        term: &SparseTerm,
        point: &Vec<ScalarField>,
    ) -> (ScalarField, Option<SparseTerm>) {
        let mut fixed_term: Option<SparseTerm> = None;
        let coeff: ScalarField =
            term.iter().fold(1u32.into(), |product, (var, power)| match *var {
                j if j == self.r_vec.len() => {
                    fixed_term = Some(SparseTerm::new(vec![(j, *power)]));
                    product
                }
                j if j < self.r_vec.len() => self.r_vec[j].pow(&[*power as u64]) * product,
                _ => point[*var - self.r_vec.len()].pow(&[*power as u64]) * product,
            });
        (coeff, fixed_term)
    }

    // g の {0,1}^v 上での全評価和を求める（遅い実装）
    pub fn slow_sum_g(&self) -> ScalarField {
        let v = self.g.num_vars();
        let n = 2u32.pow(v as u32);
        (0..n)
            .map(|n| self.g.evaluate(&n_to_vec(n as usize, v)))
            .sum()
    }
}

// 検証側の手続き

pub fn get_r() -> Option<ScalarField> {
    let mut rng = rand::thread_rng();
    let r: ScalarField = rng.gen();
    Some(r)
}

/// g の各変数に対する次数のルックアップテーブルを返す
pub fn max_degrees(g: &MultiPoly) -> Vec<usize> {
    let mut lookup: Vec<usize> = vec![0; g.num_vars()];
    g.terms().iter().for_each(|(_, term)| {
        term.iter().for_each(|(var, power)| {
            if *power > lookup[*var] {
                lookup[*var] = *power
            }
        });
    });
    lookup
}

/// プローバの主張 c_1 を検証する（ペダンティックな例）
pub fn verify(g: &MultiPoly, c_1: ScalarField) -> bool {
    // 1回目のラウンド
    let mut p = Prover::new(g);
    let mut gi = p.gen_uni_polynomial(None);
    let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
    assert_eq!(c_1, expected_c);
    let lookup_degree = max_degrees(&g);
    assert!(gi.degree() <= lookup_degree[0]);

    // 中間ラウンド
    for j in 1..p.g.num_vars() {
        let r = get_r();
        expected_c = gi.evaluate(&r.unwrap());
        gi = p.gen_uni_polynomial(r);
        let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert_eq!(expected_c, new_c);
        assert!(gi.degree() <= lookup_degree[j]);
    }
    // 最終ラウンド
    let r = get_r();
    expected_c = gi.evaluate(&r.unwrap());
    p.r_vec.push(r.unwrap());
    let new_c = p.g.evaluate(&p.r_vec);
    assert_eq!(expected_c, new_c);
    true
}

pub fn slow_verify(g: &MultiPoly, c_1: ScalarField) -> bool {
    let p = Prover::new(g);
    let manual_sum = p.slow_sum_g();
    manual_sum == c_1
}

// ────── 以下、Linear GKR プロトコルで利用する sum-check のインタラクティブプロトコル（ダミー実装） ──────

pub mod protocol {
    use ark_ff::Field;
    use rand::Rng;

    /// Sum-check プローバ側の状態（簡易な例）
    pub struct ProverState<F: Field> {
        pub num_vars: usize,
        pub current_sum: F,
    }

    /// Sum-check 検証側の状態（簡易な例）
    pub struct VerifierState<F: Field> {
        pub num_vars: usize,
        pub current_sum: F,
    }

    /// プローバ側の状態初期化
    pub fn prover_init<F: Field>(num_vars: usize, claimed_sum: F) -> ProverState<F> {
        ProverState { num_vars, current_sum: claimed_sum }
    }

    /// 現ラウンドの証明メッセージ（例：1次多項式の係数）を生成する（ダミー実装）
    pub fn prove_round<F: Field, R: Rng>(_state: &mut ProverState<F>, _rng: &mut R) -> Vec<F> {
        vec![F::zero(), F::one()]
    }

    /// プローバ側の状態を検証側のランダムチャレンジで更新（ダミー実装）
    pub fn apply_challenge<F: Field>(_state: &mut ProverState<F>, _r: F) {
        // 状態更新（ダミー実装）
    }

    /// Verifier 用のチャレンジ適用関数（状態更新はダミー）
    pub fn apply_challenge_verifier<F: Field>(_state: &mut VerifierState<F>, _r: F) {
        // 状態更新（ダミー実装）
    }

    /// 検証側の状態初期化（claimed_sum をセットする）
    pub fn verifier_init<F: Field>(num_vars: usize, claimed_sum: F) -> VerifierState<F> {
        VerifierState { num_vars, current_sum: claimed_sum }
    }

    /// 各ラウンドでプローバから送られたメッセージの検証（ダミー実装）
    pub fn verify_round<F: Field>(_state: &mut VerifierState<F>, _msg: &Vec<F>) -> Result<(), &'static str> {
        Ok(())
    }

    /// Sum-check の最終検証を行い、サブクレーム（ランダム点と期待値）を生成
    pub struct Subclaim<F: Field> {
        pub point: Vec<F>,
        pub expected_value: F,
    }
    pub fn finalize<F: Field>(state: VerifierState<F>, claimed_sum: F) -> Result<Subclaim<F>, &'static str> {
        if state.current_sum == claimed_sum {
            Ok(Subclaim { point: vec![F::one(); state.num_vars], expected_value: claimed_sum })
        } else {
            Err("Final sum-check failed")
        }
    }
}
