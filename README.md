# ABXTree
How to give ABX in ED among patients with positive chest image 

This package fit the optimal tree-shaped regime 

**For  OS X 10.11 and higher, R version 4.3.0 and higher, install gfortran-12.2-universal.pkg  from https://cran.r-project.org/bin/macosx/tools/**

```
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("yizhenxu/ABXTree")

library(ABXTree)

nb = 100 # number of burn-in samples
nd = 1000 # number of post-burn-in posterior draws 
nt = 100 # nubmer of trees
```

