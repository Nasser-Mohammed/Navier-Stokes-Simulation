## Interactive Simulation

[Live simulation](https://nasser-mohammed.github.io/simulations/programs/Navier-Stokes/index.html)

[![Simulation Preview](fluidImg.png)](https://nasser-mohammed.github.io/simulations/programs/Navier-Stokes/index.html)

## Navier-Stokes Model
$$
\begin{aligned}
&\text{Mass:} && \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0, \\[6pt]
&\text{Momentum:} && \frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{f}, \\[6pt]
&\text{Energy:} && \frac{\partial E}{\partial t} + \nabla \cdot \big((E + p)\mathbf{u}\big) = \nabla \cdot (\boldsymbol{\tau} \cdot \mathbf{u}) + \nabla \cdot (k \nabla T) + \rho \mathbf{f} \cdot \mathbf{u},
\end{aligned}
$$

where:

$$
\boldsymbol{\tau} = \mu \big(\nabla \mathbf{u} + (\nabla \mathbf{u})^T\big) - \tfrac{2}{3}\mu(\nabla \cdot \mathbf{u})\mathbf{I}
$$

is the viscous stress tensor,  
\( \rho \) is density, \( \mathbf{u} \) velocity, \( p \) pressure, \( E \) total energy,  
\( \mu \) viscosity, \( k \) thermal conductivity, \( T \) temperature, and \( \mathbf{f} \) body force (e.g. gravity).
