
<!-- badges: start -->
[![](https://img.shields.io/badge/lifecycle-stable-yellow.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/github/last-commit/psarkhosh/Stoker_solution.svg)](https://github.com/psarkhosh/Stoker_solution/commits/main)
[![License: GPL (&gt;=
3)](https://img.shields.io/badge/license-GPL%20(%3E=%203)-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=%203))
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5598374.svg)](https://doi.org/10.5281/zenodo.5598374)
<!-- badges: end -->

# Introduction
The classical dam-break problem in a rectangular prismatic channel in the absence of the source terms is recognized as the "Reimann problem" (Toro 2001). Ritter (1892) and Stoker (1957) pioneered analytical solutions to solve Saint-Venant, also known as one-dimensional shallow water, equations for the ideal dam-break problem over horizontal dry and wet beds, respectively. These solutions are widely applied as the benchmarks for validating numerical schemes to understand the model's capability to deal with continuous and discontinuous waves (Castro and Hager, 2019). For the sake of convenience, a computer code in MATLAB, titled Stoker_solution, has been developed. This source code is straightforward to simulate ideal dam-break problems over dry- and wet-beds.

![alt text](https://github.com/psarkhosh/Stoker_solution/blob/main/Initial%20condition.jpg)
![alt text](https://github.com/psarkhosh/Stoker_solution/blob/main/Output.jpg)


# Codes

- `Stoker_solution.m`: Original Matlab code by Sarkhosh (2021).
- `StokerSol.py`: Python version of Sarkhosh's code.
- `StokerSol_V3.py`: Enhanced Python code for dam-break problems.

# References
- Castro-Orgaz, O. and Hager, W.H. (2019). Shallow water hydraulics. Springer International Publishing.
- Ritter, A. (1892). Die Fortpflanzung von Wasserwellen [Propagation of water waves]. Zeitschrift Verein Deutscher Ingenieure, 36(2), 947–954 (in German).
- Stoker, J. J. (1957). Water waves: The mathematical theories with applications. New York: Wiley.
- Toro, E.F. (2001). Shock-capturing methods for free-surface shallow flows. Wiley-Blackwell.

# License 
This project is licensed under the GPL-3.0 License - see the LICENSE.md file for more details.

If using these scripts in your data analyses pipelines, please cite as:

**Sarkhosh, P. (2021). Stoker solution package (1.0.0). Zenodo. https://doi.org/10.5281/zenodo.5598374**


# Methodology of Stoker's Dam-Break Solution (Enhanced Python Code)

## Problem Setup

### Initial Conditions

Consider a rectangular channel with a dam at position $x = 0$:

- **Upstream** ($x < 0$): Initial water depth $h_0$
- **Downstream** ($x > 0$): Initial water depth $h_1$ (can be zero for dry bed)
- **Initial velocity**: $u = 0$ everywhere

### Governing Equations

The shallow water equations in conservative form:

$\frac{\partial h}{\partial t} + \frac{\partial (hu)}{\partial x} = 0 \quad \text{(Continuity)}$

$\frac{\partial (hu)}{\partial t} + \frac{\partial \left(hu^2 + \frac{1}{2}gh^2\right)}{\partial x} = 0 \quad \text{(Momentum)}$

Where:
- $h(x,t)$ = water depth
- $u(x,t)$ = depth-averaged velocity  
- $g$ = gravitational acceleration

## Method of Characteristics

### Characteristic Equations

The shallow water equations can be written in characteristic form along curves $\frac{dx}{dt} = \lambda$:

**Positive characteristics** ($\frac{dx}{dt} = u + c$):
$du + \frac{g}{c}dh = 0 \quad \text{along} \quad \frac{dx}{dt} = u + c$

**Negative characteristics** ($\frac{dx}{dt} = u - c$):
$du - \frac{g}{c}dh = 0 \quad \text{along} \quad \frac{dx}{dt} = u - c$

Where $c = \sqrt{gh}$ is the local wave speed.

### Riemann Invariants

Along characteristics, the following quantities remain constant:

- $R_+ = u + 2c$ (constant along positive characteristics)
- $R_- = u - 2c$ (constant along negative characteristics)

## Solution Structure

### Case 1: Dry Bed Downstream ($h_1 = 0$)

The solution consists of three regions:

#### Region I: Undisturbed Upstream
- **Location**: $x < -c_0 t$
- **Conditions**: $h = h_0$, $u = 0$
- **Boundary**: Negative characteristic from dam

#### Region II: Rarefaction Wave (Simple Wave)
- **Location**: $-c_0 t < x < 2c_0 t$
- **Self-similar solution**: Variables depend only on $\xi = x/t$
- **Depth**: $h = \frac{1}{9g}\left(2c_0 - \frac{x}{t}\right)^2$
- **Velocity**: $u = \frac{2}{3}\left(\frac{x}{t} + c_0\right)$

#### Region III: Dry Bed
- **Location**: $x > 2c_0 t$
- **Conditions**: $h = 0$, $u = 0$

Where $c_0 = \sqrt{gh_0}$ is the initial upstream wave speed.

### Case 2: Wet Bed Downstream ($h_1 > 0$)

The solution involves both a rarefaction wave and a shock wave:

#### Region I: Undisturbed Upstream
- Same as dry bed case

#### Region II: Rarefaction Wave
- **Location**: $-c_0 t < x < x_2(t)$
- Similar structure to dry bed case but modified by downstream conditions

#### Region III: Uniform Flow Behind Shock
- **Location**: $x_2(t) < x < x_3(t)$
- **Depth**: $h = h_A$ (intermediate depth)
- **Velocity**: $u = u_A$ (intermediate velocity)

#### Region IV: Undisturbed Downstream
- **Location**: $x > x_3(t)$
- **Conditions**: $h = h_1$, $u = 0$

## Shock Wave Analysis

### Rankine-Hugoniot Relations

For the shock wave, mass and momentum conservation give:

$[h(u - s)] = 0 \quad \text{(Mass conservation)}$

$[hu^2 + \tfrac{1}{2}gh^2 - s \cdot hu] = 0 \quad \text{(Momentum conservation)}$

Where $s$ is the shock speed and $[\cdot]$ denotes the jump across the shock.

### Shock Speed Calculation

The shock speed $C_B$ is determined by solving the nonlinear equation:

$F(C_B) = C_B h_1 - h_1 \left(\sqrt{\frac{8C_B^2}{c_1^2} + 1} - 1\right) \left(\frac{C_B}{2} - c_0 + \sqrt{\frac{gh_1}{2}\left(\sqrt{\frac{8C_B^2}{c_1^2} + 1} - 1\right)}\right) = 0$

This nonlinear equation is solved using the Newton-Raphson method.

### Intermediate State

The flow state immediately behind the shock:

$h_A = \frac{h_1}{2}\left(\sqrt{1 + \frac{8C_B^2}{c_1^2}} - 1\right)$

$u_A = 2c_0 - 2\sqrt{gh_A}$

where $c_1 = \sqrt{gh_1}$ is the downstream wave speed.

## Characteristic Structure

### Characteristic Lines in x-t Plane

The solution structure is clearly revealed in the characteristic diagram:

1. **Positive characteristics** from the dam: $x = 2c_0 t$
2. **Negative characteristics** from the dam: $x = -c_0 t$  
3. **Internal characteristics** in rarefaction: $x = \left(2c_0 - 3\sqrt{gh}\right)t$
4. **Shock trajectory**: $x = C_B t$

### Physical Interpretation

- **Rarefaction wave fan**: Expands the disturbance upstream and downstream
- **Shock wave**: Compresses the flow in wet bed cases
- **Wave interaction**: Determines the final flow configuration

## Numerical Implementation

### Newton-Raphson Iteration

For each time step, solve for shock speed using the iterative scheme:

$C_B^{(n+1)} = C_B^{(n)} - \frac{F(C_B^{(n)})}{F'(C_B^{(n)})}$

where the derivative is:

$F'(C_B) = h_1 - h_1 \left(\frac{2C_B g h_1}{c_1^2 \sqrt{\frac{8C_B^2}{c_1^2} + 1} \sqrt{\frac{gh_1}{2}\left(\sqrt{\frac{8C_B^2}{c_1^2} + 1} - 1\right)}} + \frac{1}{2}\right) \left(\sqrt{\frac{8C_B^2}{c_1^2} + 1} - 1\right)$
$- \frac{8C_B h_1}{c_1^2 \sqrt{\frac{8C_B^2}{c_1^2} + 1}} \left(\frac{C_B}{2} - c_0 + \sqrt{\frac{gh_1}{2}\left(\sqrt{\frac{8C_B^2}{c_1^2} + 1} - 1\right)}\right)$

### Spatial Discretization

For each grid point at time $t$:

1. Determine which region the point lies in based on characteristic boundaries
2. Apply appropriate analytical formula
3. Handle boundary conditions and discontinuities properly

## Key Dimensionless Parameters

### Froude Number
$\text{Fr} = \sqrt{\frac{h_0}{h_1}}$

Characterizes the relative importance of upstream to downstream conditions.

### Depth Ratio
$\eta = \frac{h_1}{h_0}$

Determines whether shock formation occurs ($\eta > 0$) or dry bed expansion ($\eta = 0$).

### Local Froude Number
$\text{Fr}_{\text{local}} = \frac{u}{\sqrt{gh}}$

Indicates subcritical ($\text{Fr} < 1$), critical ($\text{Fr} = 1$), or supercritical ($\text{Fr} > 1$) flow conditions.

## Validation and Verification

### Analytical Benchmarks

1. **Mass conservation**: Total volume preserved
   $\int_{-\infty}^{\infty} h(x,t) \, dx = \text{constant}$

2. **Energy considerations**: Monotonic energy decrease
   $E(t) = \int_{-\infty}^{\infty} \left(\frac{1}{2}hu^2 + \frac{1}{2}gh^2\right) dx \leq E(0)$

3. **Characteristic compatibility**: Riemann invariants preserved along characteristics

4. **Shock relations**: Rankine-Hugoniot conditions satisfied

### Limiting Cases

- **$h_1 \to 0$**: Recovers dry bed solution
- **$h_1 \to h_0$**: Minimal disturbance, weak shock  
- **Small times**: Initial value problem consistency
- **Large times**: Asymptotic behavior

## Applications and Extensions

### Practical Applications

1. **Dam safety analysis**: Flood wave prediction and inundation mapping
2. **Tsunami modeling**: Initial wave generation from submarine landslides
3. **Laboratory experiments**: Validation studies and physical modeling
4. **Numerical method verification**: Benchmark solutions for CFD codes

### Theoretical Extensions

1. **Non-rectangular channels**: Variable cross-sections and complex geometries
2. **Bottom friction**: Energy dissipation effects using Manning's equation
3. **Non-hydrostatic effects**: Dispersive corrections for finite amplitude waves
4. **Multiple dam-breaks**: Interaction problems and wave superposition

## Advantages and Limitations

### Advantages

- **Exact analytical solution**: No numerical approximation errors
- **Physical insight**: Clear wave structure understanding through characteristics
- **Computational efficiency**: Instantaneous evaluation at any point
- **Benchmark quality**: Reference standard for numerical codes

### Limitations

- **Idealized conditions**: Instantaneous dam removal, rectangular channel
- **Shallow water assumption**: Vertical acceleration neglected ($w \ll u$)
- **Inviscid flow**: Friction effects ignored
- **One-dimensional**: Lateral variations not considered

## Implementation Notes

### Computational Considerations

1. **Singularity at $t = 0$**: Requires small initial time offset $t_0 = \epsilon > 0$
2. **Newton-Raphson convergence**: Requires good initial guess $C_B^{(0)} \approx \sqrt{g(h_0 + h_1)/2}$
3. **Characteristic tracking**: Proper region identification using characteristic boundaries
4. **Shock capturing**: Discontinuity representation without numerical oscillations

### Accuracy and Stability

- Solution is **unconditionally stable** (analytical, no CFL restriction)
- **Machine precision accuracy** achieved (limited only by floating-point arithmetic)
- **Real-time computation** possible for engineering applications
- **Parameter sensitivity** well-controlled through dimensionless analysis

## Conclusion

Stoker's analytical solution provides a complete, exact description of the dam-break problem under idealized conditions. The method of characteristics reveals the underlying wave physics, while the shock wave analysis handles the nonlinear steepening effects. This solution serves as both a fundamental theoretical result and a practical computational tool for dam-break analysis.

The methodology demonstrates the power of analytical techniques in fluid mechanics, providing physical insight that purely numerical methods cannot easily achieve. Modern implementations combine this analytical framework with advanced visualization techniques to create powerful tools for research and engineering applications.

## References

1. Stoker, J.J. (1957). *Water Waves: The Mathematical Theory with Applications*. Interscience Publishers.
2. Toro, E.F. (2001). *Shock-Capturing Methods for Free-Surface Shallow Flows*. Wiley.
3. Chanson, H. (2004). *The Hydraulics of Open Channel Flow*. Butterworth-Heinemann.
4. Ritter, A. (1892). Die Fortpflanzung der Wasserwellen. *Zeitschrift des Vereines Deutscher Ingenieure*, 36, 947-954.
5. Whitham, G.B. (1974). *Linear and Nonlinear Waves*. Wiley-Interscience.

## Mathematical Appendix

### Complete Solution Formulas

#### Dry Bed Case ($h_1 = 0$)

**Rarefaction region** ($-c_0 t \leq x \leq 2c_0 t$):
$h(x,t) = \frac{1}{9g}\left(2c_0 - \frac{x}{t}\right)^2$

$u(x,t) = \frac{2}{3}\left(\frac{x}{t} + c_0\right)$

$c(x,t) = \frac{1}{3}\left(2c_0 - \frac{x}{t}\right)$

**Wave boundaries**:
- Head: $x = 2c_0 t$ (fastest characteristic)
- Tail: $x = -c_0 t$ (slowest characteristic)

#### Wet Bed Case ($h_1 > 0$)

**Intermediate state behind shock**:
$h_A = \frac{h_1}{2}\left(\sqrt{1 + \frac{8C_B^2}{c_1^2}} - 1\right)$

$u_A = 2c_0 - 2\sqrt{gh_A}$

**Shock speed equation** (implicit):
$F(C_B) = C_B h_1 - h_1 \left(\sqrt{\frac{8C_B^2}{c_1^2} + 1} - 1\right) \left(\frac{C_B}{2} - c_0 + \sqrt{\frac{gh_1}{2}\left(\sqrt{\frac{8C_B^2}{c_1^2} + 1} - 1\right)}\right) = 0$

**Newton-Raphson iteration**:
$C_B^{(n+1)} = C_B^{(n)} - \frac{F(C_B^{(n)})}{F'(C_B^{(n)})}$

### Characteristic Equations Summary

| Wave Type | Characteristic Speed | Riemann Invariant | Physical Meaning |
|-----------|---------------------|-------------------|------------------|
| Positive (C₊) | $\frac{dx}{dt} = u + c$ | $R_+ = u + 2c$ | Rightward information |
| Negative (C₋) | $\frac{dx}{dt} = u - c$ | $R_- = u - 2c$ | Leftward information |
| Shock | $\frac{dx}{dt} = C_B$ | Rankine-Hugoniot | Nonlinear steepening |

### Dimensionless Analysis

**Dimensionless variables**:
- Time: $\tau = \frac{tc_0}{L}$ (convective time scale)
- Space: $\xi = \frac{x}{L}$ (normalized position)  
- Velocity: $U = \frac{u}{c_0}$ (Mach number analog)
- Depth: $H = \frac{h}{h_0}$ (relative depth)

**Dimensionless parameters**:
- Froude number: $\text{Fr} = \frac{u}{\sqrt{gh}}$ (local flow parameter)
- Depth ratio: $\eta = \frac{h_1}{h_0}$ (boundary condition parameter)
- Aspect ratio: $\delta = \frac{h_0}{L}$ (shallow water parameter, $\delta \ll 1$)

### Energy and Momentum Analysis

**Total energy per unit width**:
$E = \int_{-\infty}^{\infty} \left(\frac{1}{2}\rho hu^2 + \frac{1}{2}\rho gh^2\right) dx$

**Energy dissipation at shock**:
$\Delta E = \frac{\rho g (h_A - h_1)^3}{4h_A h_1} > 0$

**Momentum conservation** (exact for inviscid flow):
$\frac{d}{dt}\int_{-\infty}^{\infty} \rho hu \, dx = 0$

This demonstrates that while momentum is conserved globally, energy is irreversibly lost at the shock, consistent with the second law of thermodynamics.
