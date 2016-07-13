# divergenceModelR
Pricing and estimation algorithms for affine jump diffusion models for divergence contracts (see Schneider and Trojani, ``Divergence and the Price of Uncertainty'').

This package accompanies my paper ``Modeling Divergence Swap Rates''.

# Divergence

Power divergence swaps are a generalization of the VIX-based variance swap. The power divergence family also encompasses Simple Variance Swaps and Gamma Swaps.

Power divergence swaps also allow for defining skewness and quarticitiy swaps. All such contracts are tradable in financial markets: they require a single option portfolio position and dynamic trading in the underlying.

# In development

For the time being, all filtering/likelihood routines in the package are subject to change, including major function signature changes. The divergence *pricing* algorithms, as well as conditional moment calculations in affine models can be treated as stable.

# Installation

Building this package requires the user to first install `affineModelR` and `ukfRcpp` and include linking information in the `Makevars` or `Makevars.win` files. For example, if your R libraries are in `c:/Libs/R` on your Windows box, use the following `Makevars.win` file:
```
PKG_LIBS = $(shell $(R_HOME)/bin/Rscript.exe -e "Rcpp:::LdFlags()") $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -Lc:/Libs/R/affineModelR/libs/x64 -laffineModelR -Lc:/Libs/R/ukfRcpp/libs/x64 -lukfRcpp
PKG_CXXFLAGS = -fpermissive -Ic:/Libs/R/affineModelR/include -Ic:/Libs/R/ukfRcpp/include
```

# Usage

There are two use cases for this package:

(1) To evaluate the prices of divergence and higher-order swaps in affine-model settings, with known parameters.
(2) To estimate affine jump diffusion models based on data about divergence and higher-order swaps.

In order to obtain data on power divergence swaps, users have to form appropriate option portfolios from data available to them.

# Author

This package is being developed by Piotr Orłowski from a code base started together with Andras Sali (https://github.com/andrewsali/)
