// src/verifier.rs

use ark_bls12_381::Fr as ScalarField;
use rand::Rng;
use crate::sumcheck::protocol;
use crate::prover::LinearGKRProof;

/// Linear GKR のサブクレーム。これを次層への入力または最終検証に利用する。
pub struct LinearGKRSubclaim {
    pub u: Vec<ScalarField>,
    pub v: Vec<ScalarField>,
    pub expected_value: ScalarField,
}

/// Linear GKR Verifier
pub struct LinearGKRVerifier;

impl LinearGKRVerifier {
    /// f2_num_vars: f2（および f3）の変数数（l）
    /// claimed_sum: Phase1 で Prover が主張した総和
    /// proof: Prover からの Linear GKR 証明
    pub fn verify<R: Rng>(
        f2_num_vars: usize,
        claimed_sum: ScalarField,
        proof: &LinearGKRProof,
        _rng: &mut R,
    ) -> Result<LinearGKRSubclaim, &'static str> {
        let l = f2_num_vars;
        if proof.phase1_msgs.len() != l || proof.phase2_msgs.len() != l {
            return Err("Invalid proof length");
        }

        // ── Phase 1 の検証 ──
        let mut verifier_state1 = protocol::verifier_init(l, claimed_sum);
        let mut u = Vec::with_capacity(l);
        for msg in proof.phase1_msgs.iter() {
            protocol::verify_round(&mut verifier_state1, msg)?;
            // ダミーの乱数生成（実際は Fiat–Shamir などで生成）
            let r_i: ScalarField = crate::sumcheck::get_r().unwrap();
            u.push(r_i);
            protocol::apply_challenge_verifier(&mut verifier_state1, r_i);
        }
        let subclaim1 = protocol::finalize(verifier_state1, claimed_sum)?;
        let u_point = subclaim1.point;
        let expected_phase1_val = subclaim1.expected_value;

        // ── Phase 2 の検証 ──
        let mut verifier_state2 = protocol::verifier_init(l, expected_phase1_val);
        let mut v = Vec::with_capacity(l);
        for msg in proof.phase2_msgs.iter() {
            protocol::verify_round(&mut verifier_state2, msg)?;
            let r_j: ScalarField = crate::sumcheck::get_r().unwrap();
            v.push(r_j);
            protocol::apply_challenge_verifier(&mut verifier_state2, r_j);
        }
        let subclaim2 = protocol::finalize(verifier_state2, expected_phase1_val)?;
        let v_point = subclaim2.point;
        let expected_phase2_val = subclaim2.expected_value;

        Ok(LinearGKRSubclaim { u: u_point, v: v_point, expected_value: expected_phase2_val })
    }
}
