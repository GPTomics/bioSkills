#!/usr/bin/env python
# Reference: chrombpnet 0.1.7+, tensorflow 2.13+, tfmodisco-lite 2.0+, shap 0.43+, numpy 1.26+ | Verify API if version differs
# In silico variant-effect prediction with chromBPNet: load trained model,
# encode reference and alternate-allele sequences (2114 bp window), predict
# base-resolution profiles + counts, compute log2_fc for each variant.
# |log2_fc| > 1 indicates strong effect; agreement with EnFormer increases
# confidence (Pampari Avsec 2024/2025).

import argparse
import numpy as np
import tensorflow as tf
import pandas as pd


# === DNA encoding ===
def encode_dna_one_hot(seq, seq_length=2114):
    '''Encode DNA sequence as one-hot for chromBPNet input. Default 2114 bp;
    variant should be at the central position. Ambiguous bases (N) become 0.'''
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    seq = seq.upper()[:seq_length].ljust(seq_length, 'N')
    one_hot = np.zeros((seq_length, 4), dtype=np.float32)
    for i, b in enumerate(seq):
        if b in mapping:
            one_hot[i, mapping[b]] = 1.0
    return one_hot


def get_variant_seq(genome_fasta, chrom, pos_1based, ref, alt, seq_length=2114):
    '''Extract ref and alt sequences centered on variant.
    Requires pyfaidx or pysam for genome sequence retrieval (mocked here).'''
    from pyfaidx import Fasta
    genome = Fasta(genome_fasta)
    half = seq_length // 2
    start = pos_1based - half - 1   # 0-based start
    end = start + seq_length
    ref_seq = str(genome[chrom][start:end])
    # Verify ref allele matches
    pivot = half
    if ref_seq[pivot] != ref:
        raise ValueError(f'Ref mismatch at {chrom}:{pos_1based} ({ref_seq[pivot]} vs {ref})')
    alt_seq = ref_seq[:pivot] + alt + ref_seq[pivot+1:]
    return ref_seq, alt_seq


# === Predict variant effect ===
def predict_variant_effect(model, ref_seq, alt_seq):
    '''Predict reference and alternate-allele profiles + counts. Returns
    log2_fc = log2(alt_counts / ref_counts). chromBPNet model outputs:
    (profile_logits, count_logits) -- counts in log scale.'''
    ref_oh = encode_dna_one_hot(ref_seq)[None, ...]
    alt_oh = encode_dna_one_hot(alt_seq)[None, ...]

    ref_profile, ref_count_logits = model.predict(ref_oh, verbose=0)
    alt_profile, alt_count_logits = model.predict(alt_oh, verbose=0)

    # chromBPNet outputs log-scale counts; exp() to get raw counts
    ref_counts = np.exp(ref_count_logits).item()
    alt_counts = np.exp(alt_count_logits).item()
    log2_fc = np.log2(alt_counts / ref_counts)

    # Profile-level KL divergence captures shape changes (vs total count)
    ref_p = tf.nn.softmax(ref_profile, axis=1).numpy().flatten()
    alt_p = tf.nn.softmax(alt_profile, axis=1).numpy().flatten()
    eps = 1e-10
    kl_div = float(np.sum(ref_p * np.log((ref_p + eps) / (alt_p + eps))))

    return {'log2_fc': float(log2_fc),
            'ref_counts': ref_counts,
            'alt_counts': alt_counts,
            'profile_kl': kl_div}


def classify_effect(log2_fc):
    '''Pampari Avsec 2024 thresholds: |log2_fc| > 1 strong; 0.3-1 moderate; <0.3 weak.'''
    abs_lfc = abs(log2_fc)
    if abs_lfc > 1.0:
        return 'strong'
    elif abs_lfc > 0.3:
        return 'moderate'
    else:
        return 'weak'


# === Ensemble for uncertainty ===
def predict_with_ensemble(model_paths, ref_seq, alt_seq):
    '''Train ensemble of chromBPNet models with different seeds; ensemble mean
    and std flag extrapolation regions where single-model predictions are
    unreliable.'''
    predictions = []
    for mpath in model_paths:
        model = tf.keras.models.load_model(mpath, compile=False)
        pred = predict_variant_effect(model, ref_seq, alt_seq)
        predictions.append(pred['log2_fc'])
    return {'ensemble_mean': float(np.mean(predictions)),
            'ensemble_std': float(np.std(predictions)),
            'ensemble_predictions': predictions,
            'high_confidence': float(np.std(predictions) < 0.2)}


# === Main ===
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='chromBPNet variant effect prediction')
    parser.add_argument('--model', required=True, help='Trained chromBPNet .h5 model')
    parser.add_argument('--ensemble', nargs='+', help='Optional: multiple model paths for ensemble')
    parser.add_argument('--variants', required=True, help='TSV: chrom, pos, ref, alt')
    parser.add_argument('--genome', required=True, help='Genome FASTA (indexed)')
    parser.add_argument('--out', required=True, help='Output TSV with predicted effects')
    args = parser.parse_args()

    model = tf.keras.models.load_model(args.model, compile=False)
    variants = pd.read_csv(args.variants, sep='\t')

    results = []
    for _, row in variants.iterrows():
        try:
            ref_seq, alt_seq = get_variant_seq(args.genome, row['chrom'],
                                                 row['pos'], row['ref'], row['alt'])
            pred = predict_variant_effect(model, ref_seq, alt_seq)
            pred['classification'] = classify_effect(pred['log2_fc'])

            if args.ensemble:
                ensemble = predict_with_ensemble(args.ensemble, ref_seq, alt_seq)
                pred.update(ensemble)

            results.append({**row.to_dict(), **pred})
        except ValueError as e:
            print(f'Skipping {row["chrom"]}:{row["pos"]}: {e}')

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.out, sep='\t', index=False)

    # Summary
    strong = (results_df['classification'] == 'strong').sum()
    moderate = (results_df['classification'] == 'moderate').sum()
    print(f'Total variants: {len(results_df)}')
    print(f'Strong effect (|log2_fc| > 1): {strong}')
    print(f'Moderate effect (0.3 < |log2_fc| < 1): {moderate}')
    # Concordance with EnFormer increases confidence for strong-effect variants.
    # Validate experimentally for high-stakes claims.
