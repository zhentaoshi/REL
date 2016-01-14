---
title: "REL and Boosting Package"
author: Zhentao Shi
date: January 2016
output: pdf_document
---
This document briefly describes the code that implement the REL and the boosting-type greedy algorithm in my paper ["Econometric Estimation with High-Dimensional Moment Equalities"](http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2491102) via a demonstration with the Chinese trade data.


## Prerequisite
To run the code, we need to install [**MOSEK**](https://www.mosek.com/) and hook it up with Matlab. **MOSEK** is a commericial convex problem solver. It provides free academic license.


## The Master Files
* **main_REL.m**: the master function for the relaxed empirical likelihood.
* **main_boosting:m**: the master function for the boosting procedure in my paper.
* Both master files depend on **trade_data.mat**, the data file for the Chinese empirical application in my paper.

## Environment
The master files work on the author's PC in Matlab 2014b


