

## Euler equations ##

1D Euler equations -

<p align="center"><img src="/tex/742f2730e9d10f5118ae00c688cf9e93.svg?invert_in_darkmode&sanitize=true" align=middle width=101.86482734999998pt height=33.81208709999999pt/></p>

where 

<p align="center"><img src="/tex/281128f2928eb2a8cc77be057538e022.svg?invert_in_darkmode&sanitize=true" align=middle width=96.39820905pt height=59.1786591pt/></p>

<p align="center"><img src="/tex/cd9718e24b48364512658c6234c7dc06.svg?invert_in_darkmode&sanitize=true" align=middle width=147.04689779999998pt height=59.1786591pt/></p>

<p align="center"><img src="/tex/be8cda4d064bc4fd70c8ad63969a638d.svg?invert_in_darkmode&sanitize=true" align=middle width=178.54961024999997pt height=40.11819404999999pt/></p>

And <p align="center"><img src="/tex/de29eb8e83da5d35e219754ed00dc2db.svg?invert_in_darkmode&sanitize=true" align=middle width=52.3461345pt height=13.789957499999998pt/></p> for air. Also important is the use of the sound speed

<p align="center"><img src="/tex/fc4385f538bc5d880531aa28c6a4661e.svg?invert_in_darkmode&sanitize=true" align=middle width=68.68477109999999pt height=39.452455349999994pt/></p>

## MHD equations ##

The ideal MHD equations can be given by - 

The 1D ideal MHD equations can be given as - 

## Approximate fluxes at Riemmann interfaces ##

## Discontinuous Galerkin Method ##
### 1D set of conservation laws ###


We seek to find weak solutions to this pde such that for a function <p align="center"><img src="/tex/75d4b214207743cc3c4ee3414bf5faed.svg?invert_in_darkmode&sanitize=true" align=middle width=9.794543549999998pt height=14.611878599999999pt/></p>,

<p align="center"><img src="/tex/c631fc118b491f8ad5efdce4d233f3d0.svg?invert_in_darkmode&sanitize=true" align=middle width=186.08054355pt height=39.452455349999994pt/></p>

In DG method we assume the solution at element <p align="center"><img src="/tex/28b0b71e05b371ee11b28542314966d1.svg?invert_in_darkmode&sanitize=true" align=middle width=9.07536795pt height=11.4155283pt/></p> to take the form

<p align="center"><img src="/tex/5441dd1f2bd361ce4e7659705cf67bb9.svg?invert_in_darkmode&sanitize=true" align=middle width=159.99420689999997pt height=50.37229065pt/></p>

Where <p align="center"><img src="/tex/4cfc2128c77b2339681dd23b9aefefd0.svg?invert_in_darkmode&sanitize=true" align=middle width=38.90136855pt height=17.031940199999998pt/></p> are the basis functions we use to represent the solution in element <p align="center"><img src="/tex/28b0b71e05b371ee11b28542314966d1.svg?invert_in_darkmode&sanitize=true" align=middle width=9.07536795pt height=11.4155283pt/></p> and <p align="center"><img src="/tex/8b7972d36baa9abc56e39ec1c328cfa1.svg?invert_in_darkmode&sanitize=true" align=middle width=31.817979599999997pt height=15.9817185pt/></p> are the expansion coefficients for element <img src="/tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode&sanitize=true" align=middle width=9.075367949999992pt height=22.831056599999986pt/> for the <img src="/tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode&sanitize=true" align=middle width=7.710416999999989pt height=21.68300969999999pt/>th order
polynomial.

Rewriting the weak form of the pde for one element <p align="center"><img src="/tex/28b0b71e05b371ee11b28542314966d1.svg?invert_in_darkmode&sanitize=true" align=middle width=9.07536795pt height=11.4155283pt/></p>,
<p align="center"><img src="/tex/ebc96065050be1dbc48d685c0eb10476.svg?invert_in_darkmode&sanitize=true" align=middle width=198.73947719999998pt height=40.55745375pt/></p>

Which can be written as 

<p align="center"><img src="/tex/f2eae65aadab22661e32b78904375726.svg?invert_in_darkmode&sanitize=true" align=middle width=340.9418529pt height=39.25784444999999pt/></p>

Which can be written as 

<p align="center"><img src="/tex/76edbc54d4f9782f93892fe70f84cb82.svg?invert_in_darkmode&sanitize=true" align=middle width=134.74387739999997pt height=33.81208709999999pt/></p>

Which is what we will solve for.

