## Workflow for simulations in PO2PLS

Please also read the first few lines of each script for extra information. All scripts are slightly modified to allow for running on a modest laptop. The full simulations should be performed on a computing server.

### For main, rank, distribution and case-control

1. Edit the file simuPO2PLS.R to choose the desired scenario and parameters
2. Run/source simuPO2PLS_aggr.R to generate plots


### For CPU times and memory usage

Run/source simuPO2PLS_CPUmem.R to evaluate computational times and plots


### For type I error and power 

1. Edit the file simuPO2PLS_asymp.R to choose the desired parameters
2. Run/source simuPO2PLS_asymp_aggr.R to generate plots

