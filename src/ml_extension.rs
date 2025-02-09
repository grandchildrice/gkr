// src/ml_extension.rs

use ark_ff::Field;
use std::collections::HashMap;

/// 密な multilinear extension
#[derive(Clone)]
pub struct DenseMLE<F: Field> {
    pub num_vars: usize,
    /// {0,1}^num_vars 上の評価結果（長さは 2^num_vars）
    pub evaluations: Vec<F>,
}

impl<F: Field> DenseMLE<F> {
    pub fn from_evaluations_vec(num_vars: usize, evaluations: Vec<F>) -> Self {
        assert_eq!(evaluations.len(), 1 << num_vars);
        DenseMLE { num_vars, evaluations }
    }
    
    /// 固定された {0,1}^n 上の評価とみなし，点 evaluation を返す（単純な実装）
    pub fn evaluate(&self, point: &[F]) -> F {
        assert_eq!(point.len(), self.num_vars);
        let mut index = 0;
        for (i, bit) in point.iter().enumerate() {
            if !bit.is_zero() {
                index |= 1 << (self.num_vars - 1 - i);
            }
        }
        self.evaluations[index]
    }
    
    /// 全評価に対してスカラー倍を実施
    pub fn scale(&mut self, scalar: F) {
        for e in self.evaluations.iter_mut() {
            *e *= scalar;
        }
    }
}

/// 疎な multilinear extension（インデックス→値のマップで表現）
#[derive(Clone)]
pub struct SparseMLE<F: Field> {
    pub num_vars: usize,
    pub evaluations: HashMap<usize, F>,
}

impl<F: Field> SparseMLE<F> {
    /// 固定する変数は先頭から固定すると仮定
    pub fn fix_variables(&self, fixed: &[F]) -> Self {
        let fixed_count = fixed.len();
        assert!(fixed_count <= self.num_vars);
        let new_num_vars = self.num_vars - fixed_count;
        let mut new_evals = HashMap::new();
        // 各評価について，先頭 fixed_count 個が fixed と一致する場合のみ残す
        for (&index, &val) in self.evaluations.iter() {
            let mut valid = true;
            for i in 0..fixed_count {
                let bit = (index >> (self.num_vars - 1 - i)) & 1;
                let fixed_bit = if fixed[i].is_zero() { 0 } else { 1 };
                if bit != fixed_bit {
                    valid = false;
                    break;
                }
            }
            if valid {
                let mut new_index = 0;
                for i in fixed_count..self.num_vars {
                    let bit = (index >> (self.num_vars - 1 - i)) & 1;
                    new_index = (new_index << 1) | bit;
                }
                new_evals.insert(new_index, val);
            }
        }
        SparseMLE { num_vars: new_num_vars, evaluations: new_evals }
    }
    
    /// 疎表現を密な multilinear extension に変換
    pub fn to_dense_multilinear_extension(&self) -> DenseMLE<F> {
        let size = 1 << self.num_vars;
        let mut evaluations = vec![F::zero(); size];
        for i in 0..size {
            if let Some(val) = self.evaluations.get(&i) {
                evaluations[i] = *val;
            }
        }
        DenseMLE { num_vars: self.num_vars, evaluations }
    }
}

