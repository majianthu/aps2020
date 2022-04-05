# aps2020
This is the code for the paper 'Variable Selection with Copula Entropy' published on Chinese Journal of Applied Probability and Statistics. The preprint paper is available at [here](https://arxiv.org/abs/1910.12389) on ArXiv.

* Ma, Jian. “Variable Selection with Copula Entropy.” Chinese Journal of Applied Probability and Statistics, 2021, 37(4): 405-420. See also arXiv preprint arXiv:1910.12389 (2019).

In the paper, three methods for variable selection are compared on the UCI [heart disease data](http://archive.ics.uci.edu/ml/datasets/heart+disease):
* Copula Entropy [1],
* Hilbert-Schimdt Independence Criterion (HSIC) [2,3],
* Distance Correlation [4].

 The following additional independence measures are also considered in this version of comparison experiment:
* Heller-Heller-Gorfine Tests of Independence [5],
* Hoeffding's D test [6],
* Bergsma-Dassios T* sign covariance [7],
* Ball correlation [8],
* BET: Binary Expansion Testing [9],
* qad: Quantification of Asymmetric Dependence [10],
* MixedIndTests [11],
* NNS: Nonlinear Nonparametric Statistics [12],
* subcopula based dependence measures [13],
* MDM: Mutual Independence Measure [14].

Copula Entropy does better than all the others measures in terms of predictibility and interpretability.

#### References
1. Ma, J., & Sun, Z. (2011). Mutual Information Is Copula Entropy. Tsinghua Science & Technology, 16(1), 51–54. See also arXiv preprint arXiv:0808.0845 (2008).
2. Gretton, A., Fukumizu, K., Teo, C. H., Song, L., Schölkopf, B., & Smola, A. J. (2007). A Kernel Statistical Test of Independence. In Advances in Neural Information Processing Systems 20 (Vol. 20, pp. 585–592).
3. Pfister, N., Bühlmann, P., Schölkopf, B., & Peters, J. (2018). Kernel-based Tests for Joint Independence. Journal of The Royal Statistical Society Series B-Statistical Methodology, 80(1), 5–31.
4. Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances. Annals of Statistics, 35(6), 2769–2794.
5. Heller, R., Heller, Y., Kaufman, S., Brill, B., & Gorfine, M. (2016). Consistent distribution-free K-sample and independence tests for univariate random variables. Journal of Machine Learning Research, 17(1), 978–1031.
6. Hoeffding, W. (1948). A Non-Parametric Test of Independence. Annals of Mathematical Statistics, 19(4), 546–557.
7. Bergsma, W., & Dassios, A. (2014). A consistent test of independence based on a sign covariance related to Kendall’s tau. Bernoulli, 20(2), 1006–1028.
8. Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019). Ball Covariance: A Generic Measure of Dependence in Banach Space. Journal of the American Statistical Association, 115, 307-317.
9. Zhang, K. (2019).BET on Independence. Journal of the American Statistical Association, Taylor & Francis, 114, 1620-1637.
10. Junker, R. R.; Griessenberger, F. & Trutschnig, W. (2021). Estimating scale-invariant directed dependence of bivariate distributions. Computational Statistics & Data Analysis, 153, 107058.
11. Genest, C.; Nešlehová, J. G.; Rémillard, B. & Murphy, O. A. Testing for independence in arbitrary distributions. Biometrika, 2019, 106, 47-68.
12. Viole, Fred and Nawrocki, David N., Deriving Nonlinear Correlation Coefficients from Partial Moments (September 18, 2012). Available at SSRN: https://ssrn.com/abstract=2148522 or http://dx.doi.org/10.2139/ssrn.2148522
13. Arturo Erdely. A subcopula based dependence measure. Kybernetika, 53(2), 231-243, 2017.
14. Ze Jin, David S. Matteson. Generalizing Distance Covariance to Measure and Test Multivariate Mutual Dependence. arXiv preprint arXiv:1709.02532, 2017.
