

This document briefly describes the code that implement the REL and the boosting-type greedy algorithm in my paper

["Econometric Estimation with High-Dimensional Moment Equalities"](http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2491102), 2016, *Journal of Econometrics*, 195, 104-119

 via a demonstration with the Chinese trade data.


## Prerequisite
To run the code, we need to install [**MOSEK**](https://www.mosek.com/) and hook it up with Matlab. **MOSEK** is a commercial convex problem solver. It provides free academic license.


## Data and Files

* **trade_data.mat**: data file for the Chinese empirical application in my paper.
* **master_REL.m**: the master function for the relaxed empirical likelihood.
* **master_boosting:m**: the master function for the boosting procedure in my paper.
* **EKK**: generates the predicted market entry and the sale. It contains the EKK paper's micro model as a nested function.

## Environment
The master files work on the author's PC in Matlab 2014b
