#!/usr/bin/env python3

"""
NMD SHAP Analysis Script - Random Forest Version

This script takes a trained Random Forest model for NMD trigger prediction and applies SHAP analysis
to explain predictions in terms of N-terminal vs C-terminal escape mechanisms.

Usage:
    python nmd_shap_analysis_rf.py model_directory input_data.txt output_data.txt

Model directory should contain:
    - model_config.json
    - scaler.joblib  
    - random_forest_model.joblib
"""

import argparse
import json
import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import shap
import joblib
from pathlib import Path

# Import version information
sys.path.insert(0, str(Path(__file__).parent.parent))
from version import get_table_annotation_lines

def sigmoid(x):
    """Convert log-odds to probability using sigmoid function"""
    return 1 / (1 + np.exp(-x))

def load_model_components(model_dir):
    """Load Random Forest model, scaler, and configuration from directory"""
    model_dir = Path(model_dir)
    
    required_files = ['model_config.json', 'scaler.joblib', 'random_forest_model.joblib']
    missing_files = [f for f in required_files if not (model_dir / f).exists()]
    
    if missing_files:
        raise FileNotFoundError(f"Missing required files in {model_dir}: {missing_files}")
    
    with open(model_dir / 'model_config.json', 'r') as f:
        config = json.load(f)
    
    scaler = joblib.load(model_dir / 'scaler.joblib')
    model = joblib.load(model_dir / 'random_forest_model.joblib')
    
    if not isinstance(model, RandomForestClassifier):
        raise ValueError(f"Expected RandomForestClassifier, got {type(model)}")
    
    categorical_features = config.get('categorical_features', ['50nt_rule', 'has_downstream_inframeAUG'])
    continuous_features = config.get('continuous_features', [])
    
    if not continuous_features and hasattr(scaler, 'feature_names_in_'):
        continuous_features = list(scaler.feature_names_in_)
    
    all_features = categorical_features + continuous_features
    
    print(f"Loaded Random Forest model from {model_dir}")
    print(f"Trees: {model.n_estimators}, Features: {len(all_features)}")
    
    return model, scaler, config, categorical_features, continuous_features

def load_input_data(input_file):
    """Load input data from txt file"""
    try:
        for sep in ['\t', ',', ' ']:
            try:
                df = pd.read_csv(input_file, sep=sep)
                if len(df.columns) > 1:
                    break
            except:
                continue
        else:
            df = pd.read_csv(input_file)
        
        print(f"Loaded {len(df)} samples with {len(df.columns)} features")
        return df
        
    except Exception as e:
        raise Exception(f"Error loading input data from {input_file}: {e}")

def define_feature_groups():
    """Define feature groups for NMD mechanism analysis"""
    return {
        'n_terminal_rescue': {
            'features': ['CDS_position', 'dis_to_first_inframeAUG', 'dis_to_first_outframeAUG', 'downstream_inframeAUG_translationAI'],
            'description': 'N-terminal truncation rescue through downstream translation reinitiation'
        },
        'c_terminal_rescue': {
            'features': ['50nt_rule', 'dis_to_exon_end', 'exon_length', 'distance_to_stop', 'downstream_exons', 'dis_to_3utr_end'],
            'description': 'C-terminal truncation rescue through favorable exon structure and NMD rule escape'
        },
        'general_features': {
            'features': ['CAI_25codon_upstream_diff', 'phyloP', 'upstream_exons', 'AF', 'gc_content', 'LOEUF', 'PTC_translationAI', 'Mean_Expression', 'm6A_CDS_length_normalized_unconstrained', 'm6A_all_length_normalized_unconstrained'],
            'description': 'General sequence and population genetic factors'
        }
    }

def apply_shap_analysis_trigger_space(model, scaler, input_data, feature_groups, categorical_features, continuous_features):
    """Apply SHAP analysis using Random Forest model"""
    print("Applying SHAP analysis...")
    
    all_features = categorical_features + continuous_features
    missing_features = [f for f in all_features if f not in input_data.columns]
    available_features = [f for f in all_features if f in input_data.columns]
    
    if missing_features:
        print(f"Warning: Missing features: {missing_features}")
    
    X = pd.DataFrame(index=input_data.index)
    
    # Process categorical features
    for feat in categorical_features:
        if feat in input_data.columns:
            if input_data[feat].dtype == 'bool':
                X[feat] = input_data[feat].astype(int)
            elif input_data[feat].dtype == 'object':
                X[feat] = input_data[feat].map({'True': 1, 'False': 0, True: 1, False: 0}).fillna(0).astype(int)
            else:
                X[feat] = input_data[feat].fillna(0).astype(int)
        else:
            X[feat] = 0
    
    # Process continuous features
    available_continuous = [f for f in continuous_features if f in input_data.columns]
    if available_continuous:
        continuous_data = input_data[available_continuous].copy()
        
        for col in continuous_data.columns:
            if continuous_data[col].isnull().any():
                median_val = continuous_data[col].median()
                continuous_data[col] = continuous_data[col].fillna(median_val)
        
        try:
            continuous_scaled = scaler.transform(continuous_data)
            continuous_scaled_df = pd.DataFrame(
                continuous_scaled,
                index=continuous_data.index,
                columns=available_continuous
            )
            X = pd.concat([X, continuous_scaled_df], axis=1)
        except Exception as e:
            print(f"Scaling failed: {e}")
            X = pd.concat([X, continuous_data], axis=1)
    
    missing_continuous = [f for f in continuous_features if f not in input_data.columns]
    for feat in missing_continuous:
        X[feat] = 0.0
    
    X = X[all_features]
    
    if X.isnull().any().any():
        X = X.fillna(0)
    
    # Make predictions
    try:
        trigger_predictions = model.predict_proba(X)[:, 1]
    except Exception as e:
        print(f"predict_proba failed: {e}")
        trigger_predictions = model.predict(X)
    
    if trigger_predictions.min() < 0 or trigger_predictions.max() > 1:
        raise ValueError(f"Invalid predictions: range {trigger_predictions.min():.3f} to {trigger_predictions.max():.3f}")
    
    escape_predictions = 1 - trigger_predictions
    
    # Calculate SHAP values
    print("Computing SHAP values...")
    explainer = shap.TreeExplainer(model)
    trigger_shap_values = explainer.shap_values(X)
    
    if isinstance(trigger_shap_values, list) and len(trigger_shap_values) == 2:
        trigger_shap_values = trigger_shap_values[1]
    elif hasattr(trigger_shap_values, 'ndim') and trigger_shap_values.ndim == 3:
        if trigger_shap_values.shape[2] == 2:
            trigger_shap_values = trigger_shap_values[:, :, 1]
    
    trigger_baseline_raw = explainer.expected_value
    
    if isinstance(trigger_baseline_raw, (list, np.ndarray)):
        if len(trigger_baseline_raw) == 2:
            trigger_baseline_raw = trigger_baseline_raw[1]
        elif hasattr(trigger_baseline_raw, 'shape') and trigger_baseline_raw.ndim > 0:
            trigger_baseline_raw = trigger_baseline_raw[1] if len(trigger_baseline_raw) > 1 else trigger_baseline_raw[0]
    
    if trigger_baseline_raw < 0 or trigger_baseline_raw > 1:
        trigger_baseline = sigmoid(trigger_baseline_raw)
        baseline_was_logodds = True
    else:
        trigger_baseline = trigger_baseline_raw
        baseline_was_logodds = False
    
    escape_baseline = 1 - trigger_baseline
    
    # Verify SHAP additivity
    if baseline_was_logodds:
        reconstructed_logits = trigger_baseline_raw + trigger_shap_values.sum(axis=1)
        reconstructed_probs = sigmoid(reconstructed_logits)
        max_error = np.abs(reconstructed_probs - trigger_predictions).max()
    else:
        reconstructed_probs = trigger_baseline + trigger_shap_values.sum(axis=1)
        max_error = np.abs(reconstructed_probs - trigger_predictions).max()
    
    print(f"SHAP additivity check: max error = {max_error:.10f}")
    
    return {
        'trigger_predictions': trigger_predictions,
        'escape_predictions': escape_predictions,
        'trigger_shap_values': trigger_shap_values,
        'trigger_baseline': trigger_baseline,
        'escape_baseline': escape_baseline,
        'feature_names': list(X.columns),
        'selected_features': available_features,
        'missing_features': missing_features,
        'categorical_features': categorical_features,
        'continuous_features': continuous_features,
        'baseline_was_logodds': baseline_was_logodds,
        'raw_baseline': trigger_baseline_raw
    }

def calculate_group_contributions_trigger_space(shap_results, feature_groups):
    """Calculate contributions from each feature group"""
    print("Calculating group contributions...")
    
    trigger_shap_values = shap_results['trigger_shap_values']
    feature_names = shap_results['feature_names']
    n_samples = len(trigger_shap_values)
    
    trigger_contributions = {
        'n_terminal_trigger_contrib': np.zeros(n_samples),
        'c_terminal_trigger_contrib': np.zeros(n_samples),
        'general_trigger_contrib': np.zeros(n_samples)
    }
    
    group_mapping = {
        'n_terminal_rescue': 'n_terminal_trigger_contrib',
        'c_terminal_rescue': 'c_terminal_trigger_contrib', 
        'general_features': 'general_trigger_contrib'
    }
    
    for group_name, group_info in feature_groups.items():
        contrib_key = group_mapping[group_name]
        group_features = group_info['features']
        
        group_indices = []
        found_features = []
        
        for feature in group_features:
            if feature in feature_names:
                group_indices.append(feature_names.index(feature))
                found_features.append(feature)
        
        if group_indices:
            group_contributions = trigger_shap_values[:, group_indices].sum(axis=1)
            trigger_contributions[contrib_key] = group_contributions
            print(f"{group_name}: {len(found_features)} features, mean = {np.mean(group_contributions):+.4f}")
    
    escape_contributions = {
        'n_terminal_contrib': -trigger_contributions['n_terminal_trigger_contrib'],
        'c_terminal_contrib': -trigger_contributions['c_terminal_trigger_contrib'],
        'general_contrib': -trigger_contributions['general_trigger_contrib']
    }
    
    total_trigger_contrib = (trigger_contributions['n_terminal_trigger_contrib'] + 
                           trigger_contributions['c_terminal_trigger_contrib'] + 
                           trigger_contributions['general_trigger_contrib'])
    
    total_escape_contrib = -total_trigger_contrib
    
    escape_contributions['total_contrib'] = total_escape_contrib
    trigger_contributions['total_trigger_contrib'] = total_trigger_contrib
    
    return escape_contributions, trigger_contributions

def calculate_mechanism_probabilities(escape_contributions):
    """Calculate probability of N vs C terminal mechanism using softmax"""
    n_contrib = escape_contributions['n_terminal_contrib']
    c_contrib = escape_contributions['c_terminal_contrib']
    
    # Use raw contributions adjusted by temperature for softmax
    temperature = 0.1
    scores = np.stack([n_contrib / temperature, c_contrib / temperature], axis=1)
    
    # Apply softmax with numerical stability
    exp_scores = np.exp(scores - scores.max(axis=1, keepdims=True))
    probs = exp_scores / exp_scores.sum(axis=1, keepdims=True)
    
    n_prob = probs[:, 0]
    c_prob = probs[:, 1]
    
    # Handle edge case: both mechanisms promote trigger (both negative)
    both_negative = (n_contrib < 0) & (c_contrib < 0)
    n_prob[both_negative] = 0.5
    c_prob[both_negative] = 0.5
    
    return {
        'n_terminal_probability': n_prob,
        'c_terminal_probability': c_prob
    }

def calculate_nt_ct_classification(escape_contributions):
    """Calculate N-terminal vs C-terminal classification"""
    n_contrib = escape_contributions['n_terminal_contrib']
    c_contrib = escape_contributions['c_terminal_contrib']
    
    pos_n = np.maximum(0, n_contrib)
    pos_c = np.maximum(0, c_contrib)
    
    total_nt_ct_escape = pos_n + pos_c
    
    n_terminal_nt_ct_relative = np.where(total_nt_ct_escape > 0, pos_n / total_nt_ct_escape, 0)
    c_terminal_nt_ct_relative = np.where(total_nt_ct_escape > 0, pos_c / total_nt_ct_escape, 0)
    
    # Simple 3-category classification
    n_terminal_dominant = n_terminal_nt_ct_relative > c_terminal_nt_ct_relative
    c_terminal_dominant = c_terminal_nt_ct_relative > n_terminal_nt_ct_relative
    
    classification = np.where(n_terminal_dominant, 'N_terminal',
                            np.where(c_terminal_dominant, 'C_terminal', 'Uncertain'))
    
    return {
        'mechanism_classification': classification,
        'has_nt_ct_mechanisms': total_nt_ct_escape > 0
    }

def create_output_data_minimal(input_data, shap_results, escape_contributions, separate_features=False):
    """Create output dataframe with analysis results
    
    Args:
        input_data: Original input data with all features
        shap_results: SHAP analysis results
        escape_contributions: Escape contribution values
        separate_features: If True, return separate dataframes for predictions and features
        
    Returns:
        If separate_features=False (default): Minimal predictions dataframe only 
            (essential columns + analysis results: CHR, POS, REF_ALLELE, ALT_ALLELE, 
            transcript_id, gene_id, nmd_trigger_probability, mechanism_classification,
            c_terminal_probability, n_terminal_probability)
        If separate_features=True: Tuple of (predictions_df, features_df)
            - predictions_df: Minimal predictions table
            - features_df: All original input features plus SHAP contribution values
    """
    print("Creating output data...")
    
    nt_ct_results = calculate_nt_ct_classification(escape_contributions)
    mech_probs = calculate_mechanism_probabilities(escape_contributions)
    
    # Add mechanism classification and probabilities - only for NMD escape cases
    is_escape = shap_results['escape_predictions'] > 0.5
    
    mechanism_classification = np.where(
        is_escape,
        nt_ct_results['mechanism_classification'],
        None
    )
    
    n_terminal_probability = np.where(
        is_escape,
        mech_probs['n_terminal_probability'],
        np.nan
    )
    
    c_terminal_probability = np.where(
        is_escape,
        mech_probs['c_terminal_probability'],
        np.nan
    )
    
    # Create minimal predictions table
    essential_cols = ['CHR', 'POS', 'REF_ALLELE', 'ALT_ALLELE', 'transcript_id', 'gene_id']
    pred_data = {}
    
    for col in essential_cols:
        if col in input_data.columns:
            pred_data[col] = input_data[col]
    
    pred_data['nmd_trigger_probability'] = shap_results['trigger_predictions']
    pred_data['mechanism_classification'] = mechanism_classification
    pred_data['c_terminal_probability'] = c_terminal_probability
    pred_data['n_terminal_probability'] = n_terminal_probability
    
    predictions_df = pd.DataFrame(pred_data)
    
    if separate_features:
        # Create features table with all original features plus SHAP values
        features_df = input_data.copy()
        features_df['n_terminal_escape_contrib'] = escape_contributions['n_terminal_contrib']
        features_df['c_terminal_escape_contrib'] = escape_contributions['c_terminal_contrib']
        features_df['general_escape_contrib'] = escape_contributions['general_contrib']
        
        print(f"Created predictions table with {len(predictions_df)} samples and {len(predictions_df.columns)} columns")
        print(f"Created features table with {len(features_df)} samples and {len(features_df.columns)} columns")
        
        return predictions_df, features_df
    else:
        # Default behavior: return minimal prediction table only
        print(f"Created minimal predictions table with {len(predictions_df)} samples and {len(predictions_df.columns)} columns")
        
        return predictions_df

def print_summary_statistics_minimal(output_data_or_pred, escape_predictions):
    """Print summary statistics"""
    print("\n" + "="*70)
    print("MECHANISM ANALYSIS SUMMARY")
    print("="*70)
    
    # Handle both single dataframe and tuple of (predictions, features)
    if isinstance(output_data_or_pred, tuple):
        pred_data = output_data_or_pred[0]
    else:
        pred_data = output_data_or_pred
    
    n_samples = len(pred_data)
    
    print(f"\nNMD Predictions (n={n_samples}):")
    trigger_probs = pred_data['nmd_trigger_probability']
    
    print(f"  Trigger: Mean={trigger_probs.mean():.4f}, Range={trigger_probs.min():.4f}-{trigger_probs.max():.4f}")
    print(f"  Escape: Mean={escape_predictions.mean():.4f}, Range={escape_predictions.min():.4f}-{escape_predictions.max():.4f}")
    
    # Count escape cases
    n_escape = (escape_predictions > 0.5).sum()
    print(f"\nNMD Escape Cases (escape probability > 0.5): {n_escape}/{n_samples} ({n_escape/n_samples*100:.1f}%)")
    
    if n_escape > 0:
        print(f"\nMechanism Classification (escape cases only):")
        escape_classifications = pred_data.loc[escape_predictions > 0.5, 'mechanism_classification']
        class_counts = escape_classifications.value_counts()
        for mechanism, count in class_counts.items():
            percentage = (count / n_escape) * 100
            print(f"  {mechanism}: {count} ({percentage:.1f}%)")
        
        print(f"\nMechanism Probabilities (escape cases only):")
        n_probs = pred_data.loc[escape_predictions > 0.5, 'n_terminal_probability']
        c_probs = pred_data.loc[escape_predictions > 0.5, 'c_terminal_probability']
        print(f"  N-terminal: Mean={n_probs.mean():.3f}, Range={n_probs.min():.3f}-{n_probs.max():.3f}")
        print(f"  C-terminal: Mean={c_probs.mean():.3f}, Range={c_probs.min():.3f}-{c_probs.max():.3f}")

def save_output_data(output_data_or_tuple, output_file, features_output_file=None, command=None):
    """Save output data to file(s) with version annotations
    
    Args:
        output_data_or_tuple: Either a single dataframe or tuple of (predictions_df, features_df)
        output_file: Path to save predictions/main output
        features_output_file: Optional path to save features separately
    """
    sep = ',' if output_file.endswith('.csv') else '\t'
    
    if isinstance(output_data_or_tuple, tuple):
        # Separate predictions and features
        predictions_df, features_df = output_data_or_tuple
        
        # Save predictions with version annotations
        with open(output_file, 'w') as f:
            # Write version annotation lines
            for line in get_table_annotation_lines(command=command):
                f.write(line + '\n')
            
            # Write the data
            predictions_df.to_csv(f, sep=sep, index=False, float_format='%.6f')
        print(f"Predictions saved to {output_file}")
        
        # Save features if output file specified
        if features_output_file:
            sep_features = ',' if features_output_file.endswith('.csv') else '\t'
            with open(features_output_file, 'w') as f:
                # Write version annotation lines
                for line in get_table_annotation_lines(command=command):
                    f.write(line + '\n')
                
                # Write the data
                features_df.to_csv(f, sep=sep_features, index=False, float_format='%.6f')
            print(f"Features saved to {features_output_file}")
    else:
        # Single combined output with version annotations
        with open(output_file, 'w') as f:
            # Write version annotation lines
            for line in get_table_annotation_lines(command=command):
                f.write(line + '\n')
            
            # Write the data
            output_data_or_tuple.to_csv(f, sep=sep, index=False, float_format='%.6f')
        print(f"Results saved to {output_file}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Apply SHAP analysis to NMD Random Forest model'
    )
    
    parser.add_argument('model_directory', help='Directory containing model files')
    parser.add_argument('input_file', help='Input data file (.txt or .csv)')
    parser.add_argument('output_file', help='Output file (.txt or .csv)')
    parser.add_argument('--features-output', help='Optional separate file for features table (all original features + SHAP values)')
    parser.add_argument('--command', type=str, default=None, help='Command line used to run the software (for annotation purposes)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    try:
        print("Loading model...")
        model, scaler, config, categorical_features, continuous_features = load_model_components(args.model_directory)
        
        print("\nLoading input data...")
        input_data = load_input_data(args.input_file)
        
        print("\nDefining feature groups...")
        feature_groups = define_feature_groups()
        
        print("\nApplying SHAP analysis...")
        shap_results = apply_shap_analysis_trigger_space(model, scaler, input_data, feature_groups, categorical_features, continuous_features)
        
        print("\nCalculating group contributions...")
        escape_contributions, trigger_contributions = calculate_group_contributions_trigger_space(shap_results, feature_groups)
        
        print("\nCreating output...")
        # Determine if we need separate features output
        separate_features = args.features_output is not None
        output_data = create_output_data_minimal(input_data, shap_results, escape_contributions, 
                                                  separate_features=separate_features)
        
        print("\nSaving results...")
        # Reconstruct command from sys.argv if not provided
        command = args.command if args.command else ' '.join(sys.argv)
        save_output_data(output_data, args.output_file, args.features_output, command=command)
        
        print_summary_statistics_minimal(output_data, shap_results['escape_predictions'])
        
        print("\n" + "="*70)
        print("ANALYSIS COMPLETED SUCCESSFULLY")
        print("="*70)
        print(f"\nResults saved to: {args.output_file}")
        if args.features_output:
            print(f"Features saved to: {args.features_output}")
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
