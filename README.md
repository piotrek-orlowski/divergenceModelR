# divergenceModelR
Pricing and estimation algorithms for affine jump diffusion models for divergence contracts (see Schneider and Trojani, ``Divergence and the Price of Uncertainty'')

# Installation

Building this package requires the user to first install `affineModelR` and `ukfRcpp` and include linking information in the `Makevars` or `Makevars.win` files. For example, if your R libraries are in `c:/Libs/R` on your Windows box, use the following `Makevars.win` file:
```
PKG_LIBS = $(shell $(R_HOME)/bin/Rscript.exe -e "Rcpp:::LdFlags()") $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -Lc:/Libs/R/affineModelR/libs/x64 -laffineModelR -Lc:/Libs/R/ukfRcpp/libs/x64 -lukfRcpp
PKG_CXXFLAGS = -fpermissive -Ic:/Libs/R/affineModelR/include -Ic:/Libs/R/ukfRcpp/include
```
