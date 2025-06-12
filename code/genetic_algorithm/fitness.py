import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.utils import resample

from imblearn.under_sampling import RandomUnderSampler

# Function to evaluate the fitness of a feature subset using Random Forest  
# subset: Binary vector (1 indicates the feature is selected, 0 indicates it is not)  
# data: Dataset in pandas.DataFrame format, where columns represent features and the last column is the label.  
# target: vector representing the class for each sample.
# n_trees: number of trees of the forest
def evaluate_fitness(subset, data, target, n_trees, algorithm="rf"):
    if sum(subset) == 0:
        return 0  # Evitar subconjuntos vacíos
    X = data.iloc[:, subset == 1].values
    y = target.values if hasattr(target, 'values') else target
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    rus = RandomUnderSampler(random_state=1234567890)
    X_resampled, y_resampled = rus.fit_resample(X_train, y_train)

    if algorithm == "rf":
        clf = RandomForestClassifier(n_estimators=n_trees)
        clf.fit(X_resampled, y_resampled)
        result = clf.score(X_test, y_test)

    elif algorithm == "svm":
        clf = SVC(kernel='rbf')
        clf.fit(X_resampled, y_resampled)
        result = clf.score(X_test, y_test)

    else:
        print("No se ha ingresado un algoritmo de fitness válido.")
    
    return result
    