```@meta
CurrentModule = MultiSpacecraftAnalysis
```

# MultiSpacecraftAnalysis

```@index
```

```@autodocs
Modules = [MultiSpacecraftAnalysis]
```

## Installation

```julia
using Pkg
Pkg.add("MultiSpacecraftAnalysis")
```

## Reciprocal vectors

[paschmannMultispacecraftAnalysisMethods2008; Chapter 4](@citet), [paschmannAnalysisMethodsMultispacecraft2000; Chapter 14](@citet), 

```@docs
reciprocal_vector
reciprocal_vectors
position_tensor
volumetric_tensor
tetrahedron_quality
```

### Estimation of spatial gradients

```@docs
lingradest
```

> Since $\tilde{g}$ and $\tilde{\boldsymbol{V}}$ are linear functions, the calculation of spatial derivatives, such as the gradient of some scalar function or the divergence or curl of a vector function, can be done quite easily. The results are:

```math
\begin{aligned}
\nabla g \simeq \nabla \tilde{g} & =\sum_{\alpha=0}^3 \boldsymbol{k}_\alpha g_\alpha \\
\hat{\boldsymbol{e}} \cdot \nabla g \simeq \hat{\boldsymbol{e}} \cdot \nabla \tilde{g} & =\sum_{\alpha=0}^3\left(\hat{\boldsymbol{e}} \cdot \boldsymbol{k}_\alpha\right) g_\alpha \\
\nabla \cdot \boldsymbol{V} \simeq \nabla \cdot \tilde{\boldsymbol{V}} & =\sum_{\alpha=0}^3 \boldsymbol{k}_\alpha \cdot \boldsymbol{V}_\alpha \\
\nabla \times \boldsymbol{V} \simeq \nabla \times \tilde{\boldsymbol{V}} & =\sum_{\alpha=0}^3 \boldsymbol{k}_\alpha \times \boldsymbol{V}_\alpha
\end{aligned}
```

## Multi-spacecraft timing

```@docs
ConstantVelocityApproach
```


```@bibliography
```