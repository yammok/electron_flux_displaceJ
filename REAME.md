# Maxwell's displacement current density (for H2 molecule)

This code is the one to calculate and visualize the current density for one-body electronic density from the data obtained by electronic structure calculation. The electronic current density is expressed by Maxwell's displacement current density, namely, 

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?{\bf&space;j}&space;({\bf&space;x},t)&space;\equiv&space;\frac{1}{4\pi}&space;\frac{\partial}{\partial&space;t}{\bf&space;E}&space;({\bf&space;x},t)" /> ,
</div>

where, 

<div align="center">
<img src="https://latex.codecogs.com/gif.latex?{\bf&space;E}&space;({\bf&space;x},t)&space;\equiv&space;-&space;\int&space;d{\bf&space;r}&space;\frac{{\bf&space;x}-{\bf&space;r}}{|{\bf&space;x}-{\bf&space;r}|^3}&space;\rho&space;_e&space;({\bf&space;r},t)" /> . 
</div>

where, one-body electronic density can be obtained by electronic structure calculation like Gaussian, GAMESS and so on. At present, the eletronic current density can be calculated within and Hydrogen molecule and the restricted Hartree-Fock (RHF) method with STO-NG basis set. However, the above electornic current density can be extended to polyatomic molecules and arbitrary methods and basis sets within electronic structure calculation. 
