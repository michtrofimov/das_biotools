import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from concurrent.futures import ProcessPoolExecutor


class RandomForestClassifierCustom(BaseEstimator):
    """
    Custom implementation of Random Forest classifier using multiple processes for parallelization.

    Parameters:
    - n_estimators (int): Number of trees in the forest.
    - max_depth (int): Maximum depth of the trees.
    - max_features (int): Maximum number of features considered for splitting a node.
    - random_state (int): Seed used by the random number generator.

    Attributes:
    - trees (list): List of trained decision trees.
    - feat_ids_by_tree (list): List of feature indices used by each tree.

    Methods:
    - fit(X, y, n_jobs=1): Fit the model to the training data.
    - _fit_process(n_estimators, X, y, workers_id): Helper function for fitting trees in parallel.
    - predict_proba(X, n_jobs): Predict class probabilities for the input samples.
    - _predict_proba_process(n_estimators, X, workers_id): Helper function for predicting probabilities in parallel.
    - predict(X, n_jobs=1): Predict class labels for the input samples.
    """

    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=42
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X, y, n_jobs=1):
        """
        Fit the Random Forest classifier to the training data.

        Parameters:
        - X (array-like): Training data.
        - y (array-like): Target values.
        - n_jobs (int): Number of processes to use for parallelization.
        """
        self.classes_ = sorted(np.unique(y))
        n_estimators_subset_size = self.n_estimators // n_jobs
        workers_id = list(range(n_jobs))
        with ProcessPoolExecutor(n_jobs) as pool:
            results = pool.map(
                self._fit_process,
                [n_estimators_subset_size] * n_jobs,
                [X] * n_jobs,
                [y] * n_jobs,
                workers_id,
            )
            for result, feat_ids_list in results:
                self.trees.extend(list(result))
                self.feat_ids_by_tree.extend(feat_ids_list)

    def _fit_process(self, n_estimators, X, y, workers_id):
        """
        Helper function to fit trees in parallel.

        Parameters:
        - n_estimators (int): Number of trees to fit.
        - X (array-like): Training data.
        - y (array-like): Target values.
        - workers_id (int): Identifier for the current worker process.

        Returns:
        - trees (list): List of trained decision trees.
        - feat_ids_list (list): List of feature indices used by each tree.
        """
        trees = []
        feat_ids_list = []

        for i in range(n_estimators):
            np.random.seed(self.random_state + i + (n_estimators * workers_id))
            feat_ids = np.random.choice(
                X.shape[1], size=self.max_features, replace=False
            )

            if isinstance(feat_ids, int):
                feat_ids = [feat_ids]

            bootstrap_indices = np.random.choice(
                X.shape[0], size=X.shape[0], replace=True
            )

            X_bootstrap = X[bootstrap_indices]
            y_bootstrap = y[bootstrap_indices]
            tree = DecisionTreeClassifier(
                max_depth=self.max_depth,
                max_features=self.max_features,
                random_state=self.random_state,
            )
            tree.fit(
                np.reshape(X_bootstrap[:, feat_ids], (-1, len(feat_ids))),
                y_bootstrap,
            )
            trees.append(tree)
            feat_ids_list.append(feat_ids)

        return trees, feat_ids_list

    def predict_proba(self, X, n_jobs):
        """
        Predict class probabilities for the input samples.

        Parameters:
        - X (array-like): Input samples.
        - n_jobs (int): Number of processes to use for parallelization.

        Returns:
        - y_probas (array-like): Predicted probabilities for each class.
        """
        n_estimators_subset_size = self.n_estimators // n_jobs
        workers_id = list(range(n_jobs))

        y_probas = np.zeros((X.shape[0], len(self.classes_)))
        with ProcessPoolExecutor(n_jobs) as pool:
            results = pool.map(
                self._predict_proba_process,
                [n_estimators_subset_size] * n_jobs,
                [X] * n_jobs,
                workers_id,
            )
            for result in results:
                y_probas += result
        return y_probas / n_jobs

    def _predict_proba_process(self, n_estimators, X, workers_id):
        """
        Helper function to predict probabilities in parallel.

        Parameters:
        - n_estimators (int): Number of trees to use for predictions.
        - X (array-like): Input samples.
        - workers_id (int): Identifier for the current worker process.

        Returns:
        - probas (array-like): Predicted probabilities for each class.
        """
        probas = np.zeros((X.shape[0], len(self.classes_)))
        tree_indices = n_estimators * workers_id
        for tree_idx in range(tree_indices, (tree_indices + n_estimators)):
            tree = self.trees[tree_idx]
            feat_ids = self.feat_ids_by_tree[tree_idx]

            probas += tree.predict_proba(
                np.reshape(X[:, feat_ids], (-1, len(feat_ids)))
            )
        return probas / len(self.trees)

    def predict(self, X, n_jobs=1):
        """
        Predict class labels for the input samples.

        Parameters:
        - X (array-like): Input samples.
        - n_jobs (int): Number of processes to use for parallelization.

        Returns:
        - predictions (array-like): Predicted class labels.
        """
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions
