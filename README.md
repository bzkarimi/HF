# HF

A simple HF code

Copyright &copy; 2016 Borna Zandkarimi

## **General comments regarding the code**:

1. In the "makefile" please change LIB to the directory in which the LAPACK is installed.  
2. All calculations were done in atomic units.  
3. 1s primitive Gaussian functions were used.  
4. The basis sets are in a file named "Basis.txt".  
5. The coordinates of atoms are in a file named "input.xyz".  
6. The final results will be saved in a file named "Output.txt".  
7. The description of each variable can be found as comment in the header of "HF.f90".  

## **Basis.txt**:

Title  
Basis set  
atom1name   &nbsp; &nbsp;    atomcharge      &nbsp; &nbsp;   slater-exponent(for scaling)  
Gaussian-coefficient  &nbsp; &nbsp;  Gaussian-exponent  
Gaussian-coefficient  &nbsp; &nbsp;  Gaussian-exponent  
.  
.  
.  
atom2name  &nbsp; &nbsp;             atomcharge  &nbsp; &nbsp;         slater-exponent(for scaling)  
Gaussian-coefficient &nbsp; &nbsp;  Gaussian-exponent  
Gaussian-coefficient &nbsp; &nbsp;  Gaussian-exponent  
.  
.  
.  

## **Input.xyz**:

atom1 &nbsp; &nbsp;  x  &nbsp; &nbsp;  y &nbsp; &nbsp;  z  
atom2 &nbsp; &nbsp;  x  &nbsp; &nbsp;  y &nbsp;&nbsp;   z  
.   &nbsp; &nbsp;    .  &nbsp; &nbsp; . &nbsp; &nbsp;  .  
.   &nbsp; &nbsp;    .  &nbsp; &nbsp; . &nbsp; &nbsp;  .  
.   &nbsp; &nbsp;    .  &nbsp; &nbsp; . &nbsp; &nbsp;  .  
