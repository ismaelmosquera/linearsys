
## **Solving Linear Systems With C**  

### *Matrix and Vector*  


This Project is an implementation suitable to efficiently solve NxN systems of linear equations;  
four different skills are presented to do that. We first declare a Matrix and a Vector type,  
and we use them in the algorithms to solve linear systems.  
In addition, we declare and implement functions to perform the most common operations applied to Matrix and Vector, and support  
to storage too.  
The tecniques implemented to solve linear systems are:  
>  
> - Cramer�s rule.  
> - Gaussian elimination.  
> - LU decomposition.  
> - QR factorization.  
>  
  
where we use numerical analisys skills in the developed algorithms.  
QR factorization is implemented just for square matrices.  
  
Computing eigenvalues and eigenvectors:  
An Eigen is a pair eigenvalue/eigenvector.  
find the maximum Eigen usign the Power Method.  
Find the eigensystem for a square matrix using the QR Algorithm.  
  
Perform SVD factorization for a MxN matrix.  
The result of a SVD factorization is as follows:  
M= UDV^t  
where  
>  
> - U is an orthogonal matrix having the left singular vectors of M.  
> - Sigma is a diagonal matrix having the singular values of M.  
> - V is an orthogonal matrix having the right singular vectors of M.  
>  
  
Compute the pseudoinverse of a MxN matrix using SVD.  
Get the nearest orthogonal matrix to a NxN matrix using SVD.  
  
The code can be compiled under any platform having a C compiler, since we use only standard headers  
which come with all C compilers, and the rest is pure C code.  
  
All the headers in the include folder are fully documented about what each funcion does.  
