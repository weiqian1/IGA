# *IGA*: An Interactive Greedy Approach to Group Sparsity Learning in High Dimensions

### Introduction

Interactive Grouped Greedy Algorithm (*IGA*) is a forward-backward stepwise variable/group selection algorithm for sparsity learning in high dimensional regression problems. It also introduces an interactive feature to potentially incorporate human operator's expert opinion in variable/group selection.

Detailed description of *IGA* algorithm is given in a technical report shown below.

This repo contains a MATLAB demo and implementation of IGA algorithm for sparse linear regression and logistic regression. 

### License

*IGA* is released under the GPL-3 License (refer to the LICENSE file for details).

### Citing *IGA*

If you find *IGA* useful in your study, please cite:

    @article{qian2017interactive,
    title={An Interactive Greedy Approach to Group Sparsity in High Dimension},
    author={Qian, Wei and Li, Wending and Sogawa, Yasuhiro and Fujimaki, Ryohei and Yang, Xitong and Liu, Ji},
    journal={arXiv preprint arXiv:1707.02963},
    year={2017}
    }


### Running Demo:
0.	Add the folder `./MATLAB_fcn` to MATLAB search path
0.	Use `demo1.m` to run a high-dimensional linear regression example with IGA.
0.	Use `demo2.m` to run a high-dimensional logistic regression example with IGA.



