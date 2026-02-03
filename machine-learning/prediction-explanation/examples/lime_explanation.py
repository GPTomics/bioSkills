'''LIME explanations for individual omics predictions'''

import pandas as pd
import numpy as np
from lime.lime_tabular import LimeTabularExplainer
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T
y = meta.loc[X.index, 'condition'].values
class_names = np.unique(y).tolist()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1))
])
pipe.fit(X_train, y_train)

# LIME needs training data for background distribution
# mode='classification' for classifiers
explainer = LimeTabularExplainer(
    X_train.values,
    feature_names=X_train.columns.tolist(),
    class_names=class_names,
    mode='classification',
    discretize_continuous=True
)

# num_features=20: Number of features to show in explanation
# Explain sample with highest disease probability
probs = pipe.predict_proba(X_test)[:, 1]
top_idx = np.argmax(probs)

exp = explainer.explain_instance(
    X_test.iloc[top_idx].values,
    pipe.predict_proba,
    num_features=20,
    top_labels=1
)

exp.save_to_file(f'lime_explanation_sample{top_idx}.html')
print(f'Saved LIME explanation for sample {top_idx} (p={probs[top_idx]:.3f})')

# Extract as dataframe
exp_list = exp.as_list(label=1)
lime_df = pd.DataFrame(exp_list, columns=['feature_rule', 'weight'])
lime_df.to_csv(f'lime_features_sample{top_idx}.csv', index=False)
lime_df.head(10)
