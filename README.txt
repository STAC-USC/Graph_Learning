Graph Laplacian Learning (GLL) Package v2.1

This MATLAB package includes implementations of graph learning algorithms presented in [1]-[2].

[1] H. E. Egilmez, E. Pavez, and A. Ortega, "Graph learning from data under Laplacian and structural constraints," IEEE Journal of Selected Topics in Signal Processing, 2017.
   
    Arxiv version:
        H. E. Egilmez, E. Pavez, and A. Ortega, "Graph learning from data under structural and Laplacian constraints," CoRR, vol. abs/1611.05181v2,2016. 
    [Online]. Available: https://arxiv.org/abs/1611.05181

[2] H. E. Egilmez, E. Pavez, and A. Ortega, "Graph Learning from Filtered Signals: Graph System and Diffusion Kernel Identification," IEEE Transactions on Signal and Information Processing over Networks, 2018.
 
    Arxiv version:
        H. E. Egilmez, E. Pavez, and A. Ortega, "Graph Learning from Filtered Signals: Graph System and Diffusion Kernel Identification," CoRR, vol. abs/1803.02553,2018. 
    [Online]. Available: https://arxiv.org/abs/1803.02553


To install the package:
(1) Download the source files.
(2) Run the script 'start_graph_learning.m'

The demo script 'demo_animals.m' shows the usage of functions used to estimate three different graph Laplacian matrices discussed in [1].

The demo script 'demo_us_temperature.m' shows the usage of functions used to estimate combinatorial Laplacian matrices from smooth signals discussed in [2]. The code regenerates Fig.7(e) in [2].

The demo script 'demo_artificial_data_on_grid.m' implements the following steps: 
   - generates a grid graph with random edges and a $\beta$-hop filter
   - generates an artificial dataset based on the graph system, which is specified by the generated graph and filter
   - runs the iterative algorithm described in [2].
   - shows ground truth and estimated graphs as well as returns the estimated $\beta$





