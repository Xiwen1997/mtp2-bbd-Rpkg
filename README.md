# mtp2-bbd-Rpkg
R implementation on how to use bridge-block decomposition to acclerate the learning of large-scale sparse MTP2 Gaussian graphical models, formulated as
$$
\mathsf{minimize}  -\log\det\left(\boldsymbol{\Theta}\right)+\left\langle \boldsymbol{\Theta},\mathbf{S}\right\rangle +\sum_{i\neq j}\Lambda_{ij}\left|\Theta_{ij}\right|, 
$$

subject to  

$$ 
	\boldsymbol{\Theta}\succ\mathbf{0}, \text{ and } \Theta_{ij}\leq0,\forall i\neq j
$$ 

using the methods proposed in [1]. The codes contain following procedures.

(1) Generating the data.

(2) Computing thresholded graph and bridge-block decomposition.

(3) Solving Sub-problems individually using FPN solver [2].

(4) Obtaining optimal solution using methods in [1].

Please skip first step if you have data matrix or sample covariance matrix provided. The methods could significantly accelerate the convergence of existing algroithms and reduce memory cost when the thresholded graphs are sparse. 
