# *IGA*: An Interactive Greedy Approach to Group Sparsity Learning in High Dimensions

### Introduction

Interactive Grouped Greedy Algorithm (*IGA*) is a forward-backward stepwise variable/group selection algorithm for sparsity learning in high dimensional regression problems. It also introduces an interactive feature to potentially incorporate human operator's expert opinion in variable/group selection.

The description of *IGA* algorithm is given in the papers shown below.

This repo contains a MATLAB demo and implementation of IGA algorithm for sparse linear regression and logistic regression. 

### License

*IGA* is released under the GPL-3 License (refer to the LICENSE file for details).

### Papers regarding *IGA*

    @article{qian2022adapative,
    title={Adaptive algorithm for multi-armed bandit problem with high-dimensional covariates},
    author={Qian, Wei and Ing, Ching-Kang and Liu, Ji},
    journal={Journal of the American Statistical Association, in press},
    year={2022+}
    }
    
    @article{qian2019interactive,
    title={An interactive greedy approach to group sparsity in high dimensions},
    author={Qian, Wei and Li, Wending and Sogawa, Yasuhiro and Fujimaki, Ryohei and Yang, Xitong and Liu, Ji},
    journal={Technometrics},
    year={2019},
    volume={61},
    issue={3}
    }
    


### Running Demo:
0.	Add the folder `./MATLAB_fcn` to MATLAB search path
0.	Use `demo1.m` to run a high-dimensional linear regression example with IGA.
0.	Use `demo2.m` to run a high-dimensional logistic regression example with IGA.

