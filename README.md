# OrdinalEffects

Title: A Latent Causal Inference Framework for Ordinal Variables

Date: 

Authors: Scauda Martina, Giusi Moffa, Jack Kuipers

Abstract: Ordinal data, including Likert scales, economic status, and education levels are commonly encountered in applied research. Yet, existing causal methods often fail to account for the inherent order among categories, as they are primarily developed either for nominal data or for continuous data where relative magnitudes are well-defined. Hence, there is a pressing need for an order-preserving methodology to compute interventional effects between ordinal variables. Presuming a latent Gaussian Directed Acyclic Graph (DAG) model as data-generating mechanism provides one possible solution [1]. Precisely, the model assumes that ordinal variables originate from marginally discretizing at given thresholds a set of Gaussian variables, whose latent covariance matrix is constrained to satisfy the conditional independencies inherent in a DAG. Conditionally on a given latent covariance matrix and thresholds, this model leads to a closed-form function for ordinal causal effects in terms of interventional distributions in the latent space. For binary variables, this approach reduces to classical methods for causal effect estimation. When the underlying DAG is unknown, one can use the Ordinal Structural EM (OSEM) algorithm [2] to learn both a plausible latent DAG, up to an equivalence class, and the model’s parameters from observational data. Simulations demonstrate the performance of the proposed approach in estimating ordinal causal effects both for known and unknown structures of the latent graph. As an illustration of a real world use case, the method is applied to survey data of 408 patients from a study [3]  on the functional relationships between symptoms of obsessive-compulsive disorder and depression.  

Description: 

References:
       [1] Silva, R., Ghahramani, Z. (2009).  The Hidden Life of Latent Variables: Bayesian Learning with Mixed Graph Models. *Journal of Machine Learning Research*, v. 10, n. 41, 1187--1238.   
       
[2] Luo, X. G., Moffa, G., and Kuipers, J. (2021).  Learning Bayesian Networks from Ordinal Data. *Journal of Machine Learning Research*, v. 22, n. 266, 1--44.

[3] McNally, R. J., Mair, P., Mugno, B. L., and Riemann, B. C. (2017).  Co-morbid obsessive-compulsive disorder and depression: a Bayesian network approach.
   *Psychological Medicine*, v. 47, n. 7, 1204--1214.

        
