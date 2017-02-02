import numpy as np 
import matplotlib.pyplot as plt

from sklearn import cross_validation, datasets, linear_model

diabetes = datasets.load_diabetes()
X = diabetes.data[:150]
print(X)
y = diabetes.target[:150]
print(y)

lasso = linear_model.Lasso()
alphas = np.logspace(-4, 0.5, 30)

scores = list()
scores_std = list()

for alpha in alphas:
	lasso.alpha = alpha
	# Default: Split into three folds and compute the score three
	# times. Then we take the average score.
	this_scores = cross_validation.cross_val_score(lasso, X, y)
	scores.append(np.mean(this_scores))
	scores_std.append(np.std(this_scores))

plt.figure(figsize=(4,3))
plt.semilogx(alphas, scores)
plt.semilogx(alphas, np.array(scores) + np.array(scores_std / np.sqrt(len(X))), 'b--')
plt.semilogx(alphas, np.array(scores) - np.array(scores_std / np.sqrt(len(X))), 'b--')
plt.ylabel('CV Score')
plt.xlabel('alpha')
plt.axhline(np.max(scores), linestyle='--', color='.5')

######################################################
# Bonus: How much can we trust the selection of alpha?

lasso_cv = linear_model.LassoCV(alphas=alphas)
k_fold   = cross_validation.KFold(len(X), 3)

print("Alpha parameters maximizing the generalizaiton score on different subsets of data: ")
for k, (train, test) in enumerate(k_fold):
	lasso_cv.fit(X[train], y[train])
	print("[fold {0}] alpha: {1:.5f}, score: {2:.5f}".
		  format(k, lasso_cv.alpha_, lasso_cv.score(X[test], y[test])))


plt.show()





