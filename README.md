[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10102872.svg)](https://doi.org/10.5281/zenodo.10102872)

# Three phases of quantum annealing: Fast, slow, and very slow
Artur Soriani, Pierre Nazé, Marcus V. S. Bonança, Bartłomiej Gardas, and Sebastian Deffner

## Description
This is the code required to reproduce all results presented in the paper [Three phases of quantum annealing: Fast, slow, and very slow](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.105.042423) published in the journal Physical Review A. The code is written by Artur Soriani.     
Additionally, this code can be used to reproduce results presented in [Assessing the performance of quantum annealing with nonlinear driving](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.105.052442), published in the journal Physical Review A, by the same authors.

## How to run
1. To compile this code, use `g++ -std=c++14 TI_comented.cpp -o TI.x`
2. When compiled, type `./TI.x 100 log -2 2` to run the program.

## Keywords
Adiabatic quantum optimization, Quantum criticality, Quantum thermodynamics

## Citation
If you find this research useful, please cite it under:
@article{PhysRevA.105.042423,                     
  title = {Three phases of quantum annealing: Fast, slow, and very slow},                 
  author = {Soriani, Artur and Naz\'e, Pierre and Bonan\ifmmode \mbox{\c{c}}\else \c{c}\fi{}a, Marcus V. S. and Gardas, Bart\l{}omiej and Deffner, Sebastian},               
  journal = {Phys. Rev. A},             
  volume = {105},            
  issue = {4},          
  pages = {042423},              
  numpages = {12},            
  year = {2022},               
  month = {Apr},               
  publisher = {American Physical Society},              
  doi = {10.1103/PhysRevA.105.042423},              
  url = {https://link.aps.org/doi/10.1103/PhysRevA.105.042423}              
}

or 

@article{PhysRevA.105.052442,
  title = {Assessing the performance of quantum annealing with nonlinear driving},              
  author = {Soriani, Artur and Naz\'e, Pierre and Bonan\ifmmode \mbox{\c{c}}\else \c{c}\fi{}a, Marcus V. S. and Gardas, Bart\l{}omiej and Deffner, Sebastian},        
  journal = {Phys. Rev. A},                 
  volume = {105},             
  issue = {5},            
  pages = {052442},            
  numpages = {9},             
  year = {2022},              
  month = {May},             
  publisher = {American Physical Society},        
  doi = {10.1103/PhysRevA.105.052442},       
  url = {https://link.aps.org/doi/10.1103/PhysRevA.105.052442}          
}
