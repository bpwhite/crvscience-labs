//! Hybrid Continuous Neural Consensus Layer (`neural_tensor.rs`).
//!
//! Combines an 8-feature bioinformatic signal encoder with a 2-layer differentiable MLP,
//! cross-locus voucher/concordance attention bonus, and Bayesian soft-quarantine logit priors.
//! Supports online contrastive weight updates across feedback epochs and persists learned
//! weights to SQLite (`neural_weights_store`) for single-epoch warm-start convergence.
//! Optional Cargo feature `candle-backend` provides direct conversion into `candle_core::Tensor`.

use serde::{Deserialize, Serialize};

/// Number of continuous bioinformatic features fed into each Processing Neuron's consensus unit.
pub const INPUT_DIM: usize = 8;
/// Hidden dimension of the hybrid sequence evaluation MLP.
pub const HIDDEN_DIM: usize = 16;

/// 8-dimensional normalized input vector extracted per `(accession, locus)` candidate:
/// 0. `length_adequacy`: `[0, 1]` position within expected locus bp window
/// 1. `ambig_cleanliness`: `1.0 - min(1.0, ambig_fraction / 0.10)`
/// 2. `reading_frame_integrity`: `1.0` if valid ORF / rRNA structure, `0.0` if internal stops
/// 3. `gc_homogeneity`: `exp(-0.5 * gc_zscore^2)` relative to locus median GC%
/// 4. `kmer_homology_sim`: `[0, 1]` 4-mer cosine similarity to locus centroid
/// 5. `within_locus_k2p_cohesion`: `exp(-0.35 * within_locus_k2p_zscore_pos)`
/// 6. `voucher_cross_locus_link`: `[0, 1]` fraction of other loci sharing the exact specimen voucher
/// 7. `shannon_complexity`: normalized nucleotide Shannon entropy (`H / 2.0` bits)
pub type SeqFeatureVec = [f64; INPUT_DIM];

/// Serializable weights for the 2-layer consensus MLP (`8 -> 16 -> 1`) + adaptive learning state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HybridConsensusWeights {
    pub w1: [[f64; INPUT_DIM]; HIDDEN_DIM],
    pub b1: [f64; HIDDEN_DIM],
    pub w2: [f64; HIDDEN_DIM],
    pub b2: f64,
    /// Gain applied to historical SQLite `skepticism_prior` penalty in logit space.
    pub quarantine_prior_gain: f64,
    pub learning_rate: f64,
}

impl Default for HybridConsensusWeights {
    fn default() -> Self {
        // Calibrated decision boundary centered at weighted feature dot-product ~3.85:
        // - Clean orthologous sequences (dot ~4.7..5.3) activate p ~ 0.88..0.96
        // - Paralog/chimeric/discordant outliers (dot ~2.6..3.4) drop to p ~ 0.15..0.42 (< 0.62 gate)
        // - Previously soft-quarantined records receive `-quarantine_prior_gain * skepticism_prior`,
        //   requiring strong cohort concordance to rehabilitate.
        let base_importance: [f64; INPUT_DIM] = [0.55, 0.70, 1.15, 0.45, 0.95, 0.85, 0.50, 0.40];
        let mut w1 = [[0.0; INPUT_DIM]; HIDDEN_DIM];
        let mut b1 = [0.0; HIDDEN_DIM];
        let mut w2 = [0.0; HIDDEN_DIM];

        for h in 0..HIDDEN_DIM {
            let sign = if h % 2 == 0 { 1.0 } else { -1.0 };
            for i in 0..INPUT_DIM {
                let jitter = ((h * 7 + i * 13 + 3) as f64).sin() * 0.04;
                w1[h][i] = sign * (base_importance[i] * 0.35 + jitter);
            }
            // Center activation hyperplanes around weighted sum 3.85 (`0.35 * 3.85 = 1.3475`)
            b1[h] = -sign * 1.35;
            w2[h] = sign * 0.48;
        }

        Self {
            w1,
            b1,
            w2,
            b2: 0.15,
            quarantine_prior_gain: 1.85,
            learning_rate: 0.04,
        }
    }
}

impl HybridConsensusWeights {
    #[inline]
    fn sigmoid(x: f64) -> f64 {
        if x >= 0.0 {
            1.0 / (1.0 + (-x).exp())
        } else {
            let ex = x.exp();
            ex / (1.0 + ex)
        }
    }

    /// Forward pass returning `(consensus_probability, hidden_activations)`.
    pub fn forward(&self, x: &SeqFeatureVec, skepticism_prior: f64) -> (f64, [f64; HIDDEN_DIM]) {
        let mut hidden = [0.0; HIDDEN_DIM];
        let mut logit = self.b2 - self.quarantine_prior_gain * skepticism_prior;

        for h in 0..HIDDEN_DIM {
            let mut z = self.b1[h];
            for i in 0..INPUT_DIM {
                z += self.w1[h][i] * x[i];
            }
            let act = z.tanh();
            hidden[h] = act;
            logit += self.w2[h] * act;
        }

        (Self::sigmoid(logit), hidden)
    }

    /// Mean binary cross-entropy of the current weights on labelled samples (no update).
    pub fn evaluate_loss(&self, samples: &[(SeqFeatureVec, f64, f64)]) -> f64 {
        if samples.is_empty() {
            return 0.0;
        }
        samples
            .iter()
            .map(|(x, prior, target)| {
                let p = self.forward(x, *prior).0.clamp(1e-7, 1.0 - 1e-7);
                -(target * p.ln() + (1.0 - target) * (1.0 - p).ln())
            })
            .sum::<f64>()
            / samples.len() as f64
    }

    /// Online contrastive feedback step updating weights toward concordant anchors (`1.0`)
    /// and away from discordant/quarantined accessions (`0.0`).
    pub fn online_feedback_step(&mut self, samples: &[(SeqFeatureVec, f64, f64)]) -> f64 {
        if samples.is_empty() {
            return 0.0;
        }
        // Class-balanced update: positives and negatives each contribute half of the gradient.
        // Real cohorts yield many more hard-gate negatives than near-identical positives; an
        // unweighted mean drifts the bias toward rejection every epoch (self-training collapse),
        // and a single-class batch carries no information about the decision boundary at all.
        let n_pos = samples.iter().filter(|s| s.2 >= 0.5).count();
        let n_neg = samples.len() - n_pos;
        if n_pos == 0 || n_neg == 0 {
            return 0.0;
        }
        let w_pos = 0.5 / (n_pos as f64);
        let w_neg = 0.5 / (n_neg as f64);
        let inv_n = 1.0 / (samples.len() as f64);
        let mut grad_w1 = [[0.0; INPUT_DIM]; HIDDEN_DIM];
        let mut grad_b1 = [0.0; HIDDEN_DIM];
        let mut grad_w2 = [0.0; HIDDEN_DIM];
        let mut grad_b2 = 0.0;
        let mut total_loss = 0.0;

        for (x, skep_prior, target) in samples {
            let (pred, hidden) = self.forward(x, *skep_prior);
            let p_clamped = pred.clamp(1e-7, 1.0 - 1e-7);
            total_loss += -(target * p_clamped.ln() + (1.0 - target) * (1.0 - p_clamped).ln());

            let d_logit = (pred - target) * if *target >= 0.5 { w_pos } else { w_neg };
            grad_b2 += d_logit;

            for h in 0..HIDDEN_DIM {
                grad_w2[h] += d_logit * hidden[h];
                let d_hidden = d_logit * self.w2[h] * (1.0 - hidden[h] * hidden[h]);
                grad_b1[h] += d_hidden;
                for i in 0..INPUT_DIM {
                    grad_w1[h][i] += d_hidden * x[i];
                }
            }
        }

        let lr = self.learning_rate;
        let l2_decay = 1e-4;
        for h in 0..HIDDEN_DIM {
            self.w2[h] = self.w2[h] * (1.0 - lr * l2_decay) - lr * grad_w2[h];
            self.b1[h] -= lr * grad_b1[h];
            for i in 0..INPUT_DIM {
                self.w1[h][i] = self.w1[h][i] * (1.0 - lr * l2_decay) - lr * grad_w1[h][i];
            }
        }
        self.b2 -= lr * grad_b2;

        total_loss * inv_n
    }

    /// Optional RuNIE `candle-core` tensor bridge (`--features candle-backend`).
    #[cfg(feature = "candle-backend")]
    pub fn forward_candle_batch(
        &self,
        batch: &[SeqFeatureVec],
        priors: &[f64],
    ) -> candle_core::Result<candle_core::Tensor> {
        use candle_core::{Device, Tensor};
        let dev = Device::Cpu;
        let n = batch.len();
        let flat_x: Vec<f64> = batch.iter().flat_map(|r| r.iter().copied()).collect();
        let x_t = Tensor::from_vec(flat_x, (n, INPUT_DIM), &dev)?;

        let mut w1_t_flat = vec![0.0; INPUT_DIM * HIDDEN_DIM];
        for i in 0..INPUT_DIM {
            for h in 0..HIDDEN_DIM {
                w1_t_flat[i * HIDDEN_DIM + h] = self.w1[h][i];
            }
        }
        let w1_t = Tensor::from_vec(w1_t_flat, (INPUT_DIM, HIDDEN_DIM), &dev)?;
        let b1_t = Tensor::from_vec(self.b1.to_vec(), (1, HIDDEN_DIM), &dev)?;
        let w2_t = Tensor::from_vec(self.w2.to_vec(), (HIDDEN_DIM, 1), &dev)?;

        let hidden = x_t.matmul(&w1_t)?.broadcast_add(&b1_t)?.tanh()?;
        let prior_shift: Vec<f64> = priors
            .iter()
            .map(|p| self.b2 - self.quarantine_prior_gain * p)
            .collect();
        let shift_t = Tensor::from_vec(prior_shift, (n, 1), &dev)?;
        let logits = hidden.matmul(&w2_t)?.broadcast_add(&shift_t)?;
        let neg = logits.neg()?;
        let denom = (neg.exp()? + 1.0)?;
        denom.recip()
    }
}
