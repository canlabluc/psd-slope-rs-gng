""" Randomized Lasso, intead of t-tests.
Which electrodes predict age?
"""

import numpy as np
import scipy as sp
import pandas as pd
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn import cross_validation, datasets
mpl.rcParams['figure.figsize'] = (16, 10)

data = pd.read_csv('../data/pipeline-full/ya-oa-full-linreg-02-24.csv')

cols = list(data.columns.values)
cols.remove('SUBJECT')
cols.remove('CLASS')
cols.remove('AGE')
cols.remove('SEX')

X = data[cols]
y = data.AGE

alpha = -15.4
resamplings = 8

rlasso = linear_model.RandomizedLasso(alpha=alpha, n_resampling=resamplings)
rlasso.fit(X, y)

print("Features sorted by score, using {} resamplings: ".format(resamplings))
feature_list = sorted(zip(map(lambda x: round(x, 4), rlasso.scores_), cols), reverse=True)
for f in feature_list:
	print(f)















# lasso = linear_model.Lasso(max_iter=10000)
# alphas = np.logspace(-4, 1, 30)

# scores = list()
# scores_std = list()

# for alpha in alphas:
# 	lasso.alpha = alpha
# 	this_scores = cross_validation.cross_val_score(lasso, X, y)
# 	scores.append(np.mean(this_scores))
# 	scores_std.append(np.std(this_scores))

# plt.semilogx(alphas, scores)
# plt.semilogx(alphas, np.array(scores) + np.array(scores_std / np.sqrt(len(X))), 'b--')
# plt.semilogx(alphas, np.array(scores) - np.array(scores_std / np.sqrt(len(X))), 'b--')
# plt.ylabel('CV Score')
# plt.xlabel('alpha')
# plt.axhline(np.max(scores), linestyle='--', color='.5')

# print(np.max(scores))

# # plt.show()

# # Nicely prints coefficients of linear models [0].
# # [0]: http://blog.datadive.net/selecting-good-features-part-ii-linear-models-and-regularization/
# def prettyprint(coefs, names=None, sort=False, n_coefs=20):
# #     if names == None:
# #         names = ["X%s" % x for x in range(len(coefs))]
#     lst = zip(coefs, names)
#     if sort:
#         lst = sorted(lst, key = lambda x:-np.abs(x[0]))
#     return " + \n".join("%s * %s" % (round(coef, 3), name) for coef, name in lst)

# pred_train, pred_test, tar_train, tar_test = cross_validation.train_test_split(X, y, test_size=.3, random_state=123)

# model = linear_model.LassoLarsCV(cv=10, precompute=False).fit(pred_train, tar_train)

# print(prettyprint(model.coef_, X.columns, sort=True))

















