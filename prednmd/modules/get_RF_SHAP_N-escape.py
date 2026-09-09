#!/usr/bin/env python3

import argparse
import json
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import shap
import joblib
from pathlib import Path

# Import version information when the script is installed inside the predNMD
# package.  Keep a small standalone fallback for downloaded copies.
sys.path.insert(0, str(Path(__file__).parent.parent))
try:
    from version import get_table_annotation_lines
except ImportError:
    def get_table_annotation_lines(command=None):
        lines = []
        if command:
            lines.append(f"# Command: {command}")
        return lines


CDS_REGROUP_BOUNDARY = 0.5          # relative CDS position at/after which CDS_position -> C group
CDS_POSITION_FEATURE = 'CDS_position'
AUG_DISTANCE_COLUMN = 'dis_to_first_inframeAUG'
AUG_BOOLEAN_COLUMN = 'has_downstream_inframeAUG'
AUG_SENTINEL = 100000.0             # value of AUG_DISTANCE_COLUMN meaning "no downstream in-frame AUG"

# Missing values for these features have explicit biological/sentinel meanings in
# the training data.  Apply these raw values before StandardScaler.transform().
DEFAULT_SPECIAL_IMPUTATION_VALUES = {
    'AF': 1e-6,
    'downstream_inframeAUG_translationAI': 1e-4,
    'PTC_translationAI': 1e-4,
    'dis_to_first_inframeAUG': 100000.0,
    'dis_to_first_outframeAUG': 100000.0,
}

# Structural annotations that every legitimate stop-gain variant has by
# construction.  A missing value here means the annotation step failed or the row
# is not a real PTC, so it is an error rather than something to impute: filling it
# with a training median turns a broken row into a plausible-looking prediction.
REQUIRED_CONTINUOUS_FEATURES = {
    'CDS_position',
    'distance_to_stop',
    'exon_length',
    'dis_to_exon_end',
    'downstream_exons',
    'upstream_exons',
    'dis_to_3utr_end',
    'gc_content',
}

# Same reasoning on the categorical side.  0 is not a neutral default for
# 50nt_rule; it is the positive claim that the PTC fails the rule.
REQUIRED_CATEGORICAL_FEATURES = {'50nt_rule'}

# Common column names used by VEP/predNMD input tables. The minimal output
# always uses the standardized names on the left. Exact standardized names are
# preferred when more than one alias is present.
OUTPUT_IDENTIFIER_ALIASES = {
    'CHR': ['CHR', 'CHROM', '#CHROM', 'chrom', 'chr'],
    'POS': ['POS', 'Position', 'position', 'pos'],
    'REF_ALLELE': ['REF_ALLELE', 'REF', 'ref'],
    'ALT_ALLELE': ['ALT_ALLELE', 'ALT', 'alt'],
    'transcript_id': [
        'transcript_id', 'TRANSCRIPT_ID', 'Feature', 'FEATURE',
        'Transcript', 'TRANSCRIPT', 'transcript',
    ],
    'gene_id': ['gene_id', 'GENE_ID', 'Gene', 'GENE', 'gene'],
}

def sigmoid(x):
    """Convert log-odds to probability using sigmoid function"""
    return 1 / (1 + np.exp(-x))

def _load_training_medians(model_dir, config, continuous_features):
    """Load raw-value training medians from joblib, with JSON fallback."""
    configured_name = config.get('continuous_medians_file', 'continuous_medians.joblib')
    median_path = model_dir / str(configured_name)

    inline_medians = config.get('continuous_medians')
    used_joblib = median_path.is_file()

    if used_joblib:
        loaded = joblib.load(median_path)
        if isinstance(loaded, pd.Series):
            medians = loaded.copy()
        elif isinstance(loaded, dict):
            medians = pd.Series(loaded)
        else:
            try:
                medians = pd.Series(loaded, index=continuous_features)
            except Exception as exc:
                raise ValueError(
                    f"Unsupported training-median artifact in {median_path}: "
                    f"{type(loaded).__name__}"
                ) from exc
        median_source = str(median_path)
    elif isinstance(inline_medians, dict):
        medians = pd.Series(inline_medians)
        median_source = f"{model_dir / 'model_config.json'} (continuous_medians block)"
    else:
        raise FileNotFoundError(
            "The model bundle does not contain training-time continuous-feature "
            "medians. Expected continuous_medians.joblib (or a "
            "continuous_medians mapping in model_config.json). Retrain/save the "
            "model with the updated training script before inference."
        )

    if medians.index.has_duplicates:
        duplicated = medians.index[medians.index.duplicated()].unique().tolist()
        raise ValueError(
            "Saved training medians contain duplicate feature entries: "
            + ', '.join(str(feature) for feature in duplicated)
        )

    missing = [feature for feature in continuous_features if feature not in medians.index]
    if missing:
        raise ValueError(
            "Saved training medians are missing continuous features: "
            + ', '.join(missing)
        )

    medians = pd.to_numeric(
        medians.reindex(continuous_features), errors='coerce'
    ).astype(float)
    invalid = medians.index[~np.isfinite(medians.to_numpy(dtype=float))].tolist()
    if invalid:
        raise ValueError(
            "Saved training medians are non-finite for: " + ', '.join(invalid)
        )

    # The export utility writes the medians twice (joblib + an inline JSON block).
    # The joblib wins, so say so when the two copies have drifted apart.
    if used_joblib and isinstance(inline_medians, dict):
        inline_series = pd.to_numeric(
            pd.Series(inline_medians).reindex(continuous_features), errors='coerce'
        ).astype(float)
        difference = (medians - inline_series).abs().max()
        if not np.isfinite(difference) or difference > 1e-9:
            print(
                f"Warning: {median_path.name} and the continuous_medians block in "
                "model_config.json disagree; the joblib file is being used"
            )

    return medians, median_source


def load_model_components(model_dir):
    """Load Random Forest, scaler, feature manifest, and training medians."""
    model_dir = Path(model_dir).expanduser().resolve()

    required_files = ['model_config.json', 'scaler.joblib', 'random_forest_model.joblib']
    missing_files = [name for name in required_files if not (model_dir / name).is_file()]
    if missing_files:
        raise FileNotFoundError(f"Missing required files in {model_dir}: {missing_files}")

    with open(model_dir / 'model_config.json', 'r') as handle:
        config = json.load(handle)

    scaler = joblib.load(model_dir / 'scaler.joblib')
    model = joblib.load(model_dir / 'random_forest_model.joblib')
    if not isinstance(model, RandomForestClassifier):
        raise ValueError(f"Expected RandomForestClassifier, got {type(model)}")

    continuous_features = list(config.get('continuous_features', []))
    if not continuous_features and hasattr(scaler, 'feature_names_in_'):
        continuous_features = [str(value) for value in scaler.feature_names_in_]
    if not continuous_features:
        raise ValueError(
            "model_config.json does not define continuous_features and the saved "
            "scaler does not expose feature_names_in_."
        )

    # No hardcoded fallback list: whatever is not continuous is read off the saved
    # model, so a bundle can never be scored against a guessed feature manifest.
    categorical_features = config.get('categorical_features')
    if categorical_features is None:
        if not hasattr(model, 'feature_names_in_'):
            raise ValueError(
                "model_config.json does not define categorical_features and the "
                "saved Random Forest does not expose feature_names_in_ to derive "
                "them from."
            )
        continuous_set = set(continuous_features)
        categorical_features = [
            str(value) for value in model.feature_names_in_
            if str(value) not in continuous_set
        ]
        print(
            "model_config.json does not define categorical_features; derived "
            f"{categorical_features} from the saved model"
        )
    categorical_features = list(categorical_features)

    all_features = categorical_features + continuous_features
    if len(set(all_features)) != len(all_features):
        raise ValueError("The saved feature manifest contains duplicate feature names")

    training_medians, median_source = _load_training_medians(
        model_dir, config, continuous_features
    )

    special_imputation_values = dict(DEFAULT_SPECIAL_IMPUTATION_VALUES)
    configured_special = config.get('special_imputation_values', {})
    if isinstance(configured_special, dict):
        unknown = [
            str(feature) for feature in configured_special
            if str(feature) not in continuous_features
        ]
        if unknown:
            raise ValueError(
                "special_imputation_values in model_config.json names feature(s) "
                "that are not continuous model features, so they would be silently "
                "ignored: " + ', '.join(unknown)
            )
        for feature, value in configured_special.items():
            try:
                special_imputation_values[str(feature)] = float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid special imputation value for {feature!r}: {value!r}"
                ) from exc

    if hasattr(scaler, 'n_features_in_') and int(scaler.n_features_in_) != len(continuous_features):
        raise ValueError(
            f"Saved scaler expects {int(scaler.n_features_in_)} continuous features, "
            f"but model_config.json defines {len(continuous_features)}"
        )
    if hasattr(scaler, 'feature_names_in_'):
        scaler_features = [str(value) for value in scaler.feature_names_in_]
        if scaler_features != continuous_features:
            raise ValueError(
                "Saved scaler feature order does not match model_config.json. "
                f"Scaler={scaler_features}; config={continuous_features}"
            )
    if hasattr(model, 'n_features_in_') and int(model.n_features_in_) != len(all_features):
        raise ValueError(
            f"Saved Random Forest expects {int(model.n_features_in_)} total features, "
            f"but the manifest defines {len(all_features)}"
        )
    if hasattr(model, 'feature_names_in_'):
        model_features = [str(value) for value in model.feature_names_in_]
        if model_features != all_features:
            raise ValueError(
                "Saved Random Forest feature order does not match model_config.json. "
                f"Model={model_features}; config={all_features}"
            )

    print(f"Loaded Random Forest model from {model_dir}")
    print(f"Trees: {model.n_estimators}, Features: {len(all_features)}")
    print(f"Loaded training medians from: {median_source}")

    return (
        model,
        scaler,
        config,
        categorical_features,
        continuous_features,
        training_medians,
        special_imputation_values,
    )

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

def _coerce_binary_feature(series, feature, required=False):
    """Match training behavior: valid binary values, with missing -> 0.

    When ``required`` is set, a missing value is an error instead: 0 is a
    substantive claim about the variant, not a neutral placeholder.
    """
    if pd.api.types.is_bool_dtype(series.dtype):
        numeric = series.astype('Int64').astype(float)
    else:
        normalized = series.copy()
        if pd.api.types.is_object_dtype(normalized.dtype) or isinstance(
            normalized.dtype, pd.StringDtype
        ):
            text_values = normalized.astype('string').str.strip().str.lower()
            mapped = text_values.map({
                'true': 1, 'false': 0, 'yes': 1, 'no': 0, 'y': 1, 'n': 0,
            })
            normalized = normalized.where(mapped.isna(), mapped)
        numeric = pd.to_numeric(normalized, errors='coerce')

    if required and numeric.isna().any():
        rows = numeric.index[numeric.isna()]
        raise ValueError(
            f"Required binary feature {feature!r} has no usable value for "
            f"{len(rows)} row(s) (first: {rows[:5].tolist()}). Every stop-gain "
            "variant must have this feature; check the annotation step rather "
            "than defaulting it to 0."
        )

    numeric = numeric.fillna(0)
    invalid = ~numeric.isin([0, 1])
    if invalid.any():
        examples = series.loc[invalid].drop_duplicates().head(10).tolist()
        raise ValueError(
            f"Binary feature {feature!r} contains values other than 0/1: {examples}"
        )
    return numeric.astype(int)


def prepare_model_features(
    input_data,
    scaler,
    categorical_features,
    continuous_features,
    training_medians,
    special_imputation_values,
):
    """Reproduce training-time raw imputation and scaling exactly."""
    all_features = list(categorical_features) + list(continuous_features)
    missing_columns = [feature for feature in all_features if feature not in input_data.columns]
    if missing_columns:
        raise ValueError(
            "Input data is missing required model feature columns: "
            + ', '.join(missing_columns)
        )

    X = pd.DataFrame(index=input_data.index)
    for feature in categorical_features:
        X[feature] = _coerce_binary_feature(
            input_data[feature],
            feature,
            required=feature in REQUIRED_CATEGORICAL_FEATURES,
        )

    continuous_data = input_data[list(continuous_features)].copy()
    for feature in continuous_features:
        continuous_data[feature] = pd.to_numeric(
            continuous_data[feature], errors='coerce'
        )

    # Structural annotations are checked before any imputation runs: a gap here is
    # an upstream failure, not a value to fill in.
    required_gaps = {}
    for feature in continuous_features:
        if feature not in REQUIRED_CONTINUOUS_FEATURES:
            continue
        missing = continuous_data[feature].isna()
        if missing.any():
            required_gaps[feature] = continuous_data.index[missing]
    if required_gaps:
        details = '; '.join(
            f"{feature}: {len(rows):,} row(s) (first: {rows[:5].tolist()})"
            for feature, rows in required_gaps.items()
        )
        raise ValueError(
            "Required feature(s) have no usable value for some rows. Every "
            "stop-gain variant must have these by construction, so check the "
            "annotation step rather than imputing them. " + details
        )

    special_counts = {}
    for feature, fill_value in special_imputation_values.items():
        if feature not in continuous_data.columns:
            continue
        missing = continuous_data[feature].isna()
        if missing.any():
            continuous_data.loc[missing, feature] = float(fill_value)
            special_counts[feature] = int(missing.sum())

    # Any remaining missing continuous values receive the corresponding raw-value
    # median learned from the training table, never a median from this input batch.
    median_counts = {}
    for feature in continuous_features:
        missing = continuous_data[feature].isna()
        if missing.any():
            median = float(training_medians.loc[feature])
            continuous_data.loc[missing, feature] = median
            median_counts[feature] = int(missing.sum())

    remaining = continuous_data.columns[continuous_data.isna().any()].tolist()
    if remaining:
        raise ValueError(
            "Missing values remain after special-value and training-median "
            "imputation: " + ', '.join(remaining)
        )

    if special_counts:
        print("Applied project-specific missing-value replacements:")
        for feature in continuous_features:
            if feature in special_counts:
                print(
                    f"  {feature}: {special_counts[feature]:,} row(s) -> "
                    f"{special_imputation_values[feature]:g}"
                )
    if median_counts:
        print("Applied saved training medians to remaining missing values:")
        for feature in continuous_features:
            if feature in median_counts:
                print(
                    f"  {feature}: {median_counts[feature]:,} row(s) -> "
                    f"{float(training_medians.loc[feature]):g}"
                )

    # Failure is explicit: never fall back to unscaled data.
    continuous_scaled = scaler.transform(continuous_data[list(continuous_features)])
    continuous_scaled_df = pd.DataFrame(
        continuous_scaled,
        index=continuous_data.index,
        columns=list(continuous_features),
    )
    X = pd.concat([X, continuous_scaled_df], axis=1)
    X = X[all_features]

    if not X.index.equals(input_data.index):
        raise AssertionError("Row order changed during inference preprocessing")
    if X.columns.tolist() != all_features:
        raise AssertionError("Inference feature names or order changed")
    if X.isna().any().any():
        bad = X.columns[X.isna().any()].tolist()
        raise AssertionError(
            "Unexpected missing values after inference preprocessing: "
            + ', '.join(bad)
        )
    if not np.isfinite(X.to_numpy(dtype=float)).all():
        raise AssertionError("Non-finite model values remain after preprocessing")

    report = {
        'special_imputation_counts': special_counts,
        'training_median_imputation_counts': median_counts,
    }
    return X, continuous_data, report


def apply_shap_analysis_trigger_space(
    model,
    scaler,
    input_data,
    feature_groups,
    categorical_features,
    continuous_features,
    training_medians,
    special_imputation_values,
):
    """Apply RF prediction and SHAP using training-consistent preprocessing."""
    print("Applying SHAP analysis...")

    X, imputed_continuous_data, preprocessing_report = prepare_model_features(
        input_data=input_data,
        scaler=scaler,
        categorical_features=categorical_features,
        continuous_features=continuous_features,
        training_medians=training_medians,
        special_imputation_values=special_imputation_values,
    )

    trigger_predictions = model.predict_proba(X)[:, 1]
    if trigger_predictions.min() < 0 or trigger_predictions.max() > 1:
        raise ValueError(
            f"Invalid predictions: range {trigger_predictions.min():.3f} "
            f"to {trigger_predictions.max():.3f}"
        )
    escape_predictions = 1 - trigger_predictions

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
            trigger_baseline_raw = (
                trigger_baseline_raw[1]
                if len(trigger_baseline_raw) > 1
                else trigger_baseline_raw[0]
            )

    if trigger_baseline_raw < 0 or trigger_baseline_raw > 1:
        trigger_baseline = sigmoid(trigger_baseline_raw)
        baseline_was_logodds = True
    else:
        trigger_baseline = trigger_baseline_raw
        baseline_was_logodds = False
    escape_baseline = 1 - trigger_baseline

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
        'selected_features': list(X.columns),
        'missing_features': [],
        'categorical_features': list(categorical_features),
        'continuous_features': list(continuous_features),
        'baseline_was_logodds': baseline_was_logodds,
        'raw_baseline': trigger_baseline_raw,
        'imputed_continuous_data': imputed_continuous_data,
        'preprocessing_report': preprocessing_report,
        'special_imputation_values': dict(special_imputation_values),
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

    grouped_features = {
        feature for info in feature_groups.values() for feature in info['features']
    }
    unassigned = [name for name in feature_names if name not in grouped_features]
    if unassigned:
        print(
            "Warning: model feature(s) belong to no mechanism group, so their SHAP "
            "contributions are excluded from the N/C/general decomposition: "
            + ', '.join(unassigned)
        )
    
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

def no_downstream_aug_mask(input_data, special_imputation_values=None):
    """Identify rows with no downstream in-frame AUG.

    Known values in an explicit boolean column take priority.  Missing/unknown
    boolean values fall back to the AUG-distance column.  Missing AUG distance is
    interpreted using the same 100000 sentinel used for inference imputation.
    """
    n_rows = len(input_data)
    no_aug = np.zeros(n_rows, dtype=bool)
    resolved = np.zeros(n_rows, dtype=bool)

    if AUG_BOOLEAN_COLUMN in input_data.columns:
        # Parse to a number rather than matching literal spellings: a float64
        # column (which any missing value forces) renders as '1.0'/'0.0', which
        # string matching would reject, silently discarding a usable column.
        text_values = input_data[AUG_BOOLEAN_COLUMN].astype('string').str.strip().str.lower()
        mapped = text_values.map(
            {'true': 1.0, 'false': 0.0, 'yes': 1.0, 'no': 0.0, 'y': 1.0, 'n': 0.0}
        )
        numeric_values = mapped.fillna(pd.to_numeric(text_values, errors='coerce'))
        known = numeric_values.isin([0, 1]).to_numpy()
        no_aug[known] = (numeric_values == 0).to_numpy()[known]
        resolved[known] = True

    if AUG_DISTANCE_COLUMN in input_data.columns:
        distance = pd.to_numeric(
            input_data[AUG_DISTANCE_COLUMN], errors='coerce'
        )
        fill_value = AUG_SENTINEL
        if special_imputation_values is not None:
            fill_value = float(
                special_imputation_values.get(AUG_DISTANCE_COLUMN, AUG_SENTINEL)
            )
        distance = distance.fillna(fill_value).to_numpy(dtype=float)
        distance_no_aug = distance >= AUG_SENTINEL
        unresolved = ~resolved
        no_aug[unresolved] = distance_no_aug[unresolved]
        resolved[unresolved] = True

    if not resolved.all():
        print(
            f"Warning: {int((~resolved).sum()):,} row(s) lack a usable "
            f"{AUG_BOOLEAN_COLUMN!r} or {AUG_DISTANCE_COLUMN!r}; the downstream-AUG "
            "gate is left inactive for those rows"
        )
    return no_aug


def apply_mechanism_rules(escape_contributions, shap_results, input_data):
    """Apply the two attribution rules to the group contributions.

    Both are repartitions of a fixed total: a group contribution is the sum of its
    features' SHAP values, so moving a feature between groups is addition and
    n + c + general is unchanged. Neither rule re-runs the model.

    1. CDS_position moves from the N group to the C group for PTCs at or past
       CDS_REGROUP_BOUNDARY of the CDS.
    2. Where no downstream in-frame AUG exists, the whole remaining N contribution
       moves to the general group and the N contribution becomes zero.

    The gate mask is returned inside escape_contributions so that the probability
    and classification functions can honour it.
    """
    print("Applying mechanism attribution rules...")

    n_contrib = np.asarray(escape_contributions['n_terminal_contrib'], dtype=float).copy()
    c_contrib = np.asarray(escape_contributions['c_terminal_contrib'], dtype=float).copy()
    g_contrib = np.asarray(escape_contributions['general_contrib'], dtype=float).copy()

    # Rule 1: CDS_position to the C group past the midpoint of the CDS.
    feature_names = list(shap_results['feature_names'])

    # Read the post-imputation values the model was actually given, so the rule and
    # the SHAP value it is repartitioning describe the same variant.
    imputed = shap_results.get('imputed_continuous_data')

    def _position_values(column):
        if imputed is not None and column in imputed.columns:
            return imputed[column].to_numpy(dtype=float)
        if column in input_data.columns:
            return pd.to_numeric(input_data[column], errors='coerce').to_numpy(dtype=float)
        return None

    cds_pos = _position_values('CDS_position')
    dist_stop = _position_values('distance_to_stop')
    have_position = cds_pos is not None and dist_stop is not None

    if CDS_POSITION_FEATURE in feature_names and have_position:
        shap_cds = np.asarray(shap_results['trigger_shap_values'], dtype=float)[
            :, feature_names.index(CDS_POSITION_FEATURE)
        ]
        total = cds_pos + dist_stop

        with np.errstate(invalid='ignore', divide='ignore'):
            relative = np.where(total > 0, cds_pos / total, np.nan)

        # An undefined relative position keeps the original grouping rather than
        # being moved on the strength of a missing value.
        moved = relative >= CDS_REGROUP_BOUNDARY
        shift = np.where(moved & np.isfinite(shap_cds), shap_cds, 0.0)

        # The contributions are escape-space (minus the trigger-space SHAP sum), so
        # removing a feature from N adds its trigger-space value back to N.
        n_contrib = n_contrib + shift
        c_contrib = c_contrib - shift
        print(f"  CDS_position -> C group for {int(np.sum(moved)):,} variant(s) "
              f"at relative CDS position >= {CDS_REGROUP_BOUNDARY}")
        undefined = int(np.sum(~np.isfinite(relative)))
        if undefined:
            print(f"  {undefined:,} variant(s) had an undefined relative CDS position "
                  "and kept their original grouping")
    else:
        print("  CDS_position regrouping skipped: CDS_position/distance_to_stop not available")
        moved = np.zeros(len(input_data), dtype=bool)

    # Rule 2: no downstream in-frame AUG means no reinitiation, so nothing can
    # remain in the N group.
    gated = no_downstream_aug_mask(
        input_data, shap_results.get('special_imputation_values')
    )
    if np.any(gated):
        g_contrib = np.where(gated, g_contrib + n_contrib, g_contrib)
        n_contrib = np.where(gated, 0.0, n_contrib)
        print(f"  N contribution zeroed for {int(np.sum(gated)):,} variant(s) "
              "with no downstream in-frame AUG")

    escape_contributions = dict(escape_contributions)
    escape_contributions['n_terminal_contrib'] = n_contrib
    escape_contributions['c_terminal_contrib'] = c_contrib
    escape_contributions['general_contrib'] = g_contrib
    escape_contributions['no_downstream_aug'] = gated

    for group, values in (('n_terminal_rescue', n_contrib),
                          ('c_terminal_rescue', c_contrib),
                          ('general_features', g_contrib)):
        print(f"  {group} after rules: mean = {np.mean(values):+.4f}")

    return escape_contributions


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

    # No downstream in-frame AUG means N-terminal rescue is impossible. The softmax
    # compares n - c, so a zeroed N contribution would beat any negative C
    # contribution and report N as likely for exactly those variants; the
    # probabilities are therefore set directly rather than left to the softmax.
    gated = escape_contributions.get('no_downstream_aug')
    if gated is not None and np.any(gated):
        n_prob = np.where(gated, 0.0, n_prob)
        c_prob = np.where(gated, 1.0, c_prob)

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

    # With reinitiation ruled out, C-terminal is the only mechanism available to a
    # variant that escapes, so the call is forced rather than left as the
    # 'Uncertain' this rule would return where C evidence is also non-positive.
    gated = escape_contributions.get('no_downstream_aug')
    if gated is not None and np.any(gated):
        classification = np.where(gated, 'C_terminal', classification)

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
    
    # Create a stable minimal predictions table. VEP-style input commonly uses
    # Feature/Gene instead of transcript_id/gene_id, so map known aliases to
    # standardized output names rather than silently dropping those columns.
    pred_data = {}
    identifier_mapping = {}
    for output_column, aliases in OUTPUT_IDENTIFIER_ALIASES.items():
        source_column = next(
            (candidate for candidate in aliases if candidate in input_data.columns),
            None,
        )
        if source_column is None:
            pred_data[output_column] = pd.Series(
                pd.NA, index=input_data.index, dtype='object'
            )
        else:
            pred_data[output_column] = input_data[source_column]
            identifier_mapping[output_column] = source_column

    for output_column in OUTPUT_IDENTIFIER_ALIASES:
        source_column = identifier_mapping.get(output_column)
        if source_column is None:
            print(
                f"Warning: no input column was found for {output_column!r}. "
                f"Checked aliases: {OUTPUT_IDENTIFIER_ALIASES[output_column]}. "
                "The output column is present but contains missing values."
            )
        elif source_column != output_column:
            print(f"Mapped input column {source_column!r} -> output column {output_column!r}")

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
        (
            model,
            scaler,
            config,
            categorical_features,
            continuous_features,
            training_medians,
            special_imputation_values,
        ) = load_model_components(args.model_directory)
        
        print("\nLoading input data...")
        input_data = load_input_data(args.input_file)
        
        print("\nDefining feature groups...")
        feature_groups = define_feature_groups()
        
        shap_results = apply_shap_analysis_trigger_space(
            model,
            scaler,
            input_data,
            feature_groups,
            categorical_features,
            continuous_features,
            training_medians,
            special_imputation_values,
        )
        
        escape_contributions, trigger_contributions = calculate_group_contributions_trigger_space(shap_results, feature_groups)

        escape_contributions = apply_mechanism_rules(escape_contributions, shap_results, input_data)
        
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
