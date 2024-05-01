import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=42
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X, y, n_jobs):
        self.classes_ = sorted(np.unique(y))
        ## ENTER YOUR CODE HERE (/¯◡ ‿ ◡)/¯☆*##
        for i in range(self.n_estimators):
            np.random.seed(self.random_state + i)
            feat_ids = np.random.choice(X.shape[1], size=self.max_features, replace=False)
            if feat_ids.isinstance(int):
                feat_ids = [feat_ids]
            self.feat_ids_by_tree.append(feat_ids)
            bootstrap_indices = np.random.choice(X.shape[0], size=X.shape[0], replace=True)
            X_bootstrap = X[bootstrap_indices]
            y_bootstrap = y[bootstrap_indices]
            tree = DecisionTreeClassifier(
                max_depth=self.max_depth, max_features=self.max_features, random_state=self.random_state
            )
            tree.fit(np.reshape(X_bootstrap[:,feat_ids], (-1,len(feat_ids))), y_bootstrap)
            self.trees.append(tree)
            
        return self

    def predict_proba(self, X, n_jobs):
        ## ENTER YOUR CODE HERE (/¯◡ ‿ ◡)/¯☆*##
        probas = np.zeros((X.shape[0], len(self.classes_)))
        for tree, feat_ids in zip(self.trees, self.feat_ids_by_tree):
            probas += tree.predict_proba(np.reshape(X[:, feat_ids],(-1,len(feat_ids))))
        return probas / len(self.trees)

    def predict(self, X, n_jobs):
        probas = self.predict_proba(X)
        predictions = np.argmax(probas, axis=1)
        return predictions
