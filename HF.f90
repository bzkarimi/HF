	Program HF

!====================================================
!                  SOME GENERAL COMMENTS
!----------------------------------------------------
! 1. All calculations were done in atomic units.
! 2. 1s primitive Gaussian functions were used.
! 3. The basis sets are in a file named "Basis.txt".
! 4. The coordinates of atoms are in a file named "input.xyz".
! 5. The final results will be saved in a file named "Output.txt".

!====================================================

!====================================================
!                  VARIABLES(ALPHABETICAL ORDER)
!----------------------------------------------------
!
! alpha(i,j) = Guassian exponent for atom "i" and primitive "j"
! alpha1 = Scaling factor for DGEMM subroutine
! atomname = name of the atom
! beta1 = Scaling factor for DGEMM subroutine
! coeff(i,j) = Contraction coefficient for atom "i" and primitive "j"
! Converge = whenever (delta < Converge) this means that our result is in the desired interval
! d,d1,d2,d3 = distances between centers
! delta = standard deviation for change in density matrix
! Eigen_C/Eigen_C_prime = Eigen vectors of Fock/transformed Fock Matrix
! Eigen_value() = Eigenvalues of matrix
! EN = Electronic energy
! EN_t = total energy
! F(),F_Prime() = Fock Matrix before/after transformation
! G() = 2e part of Fock Matrix
! H_core() = core H matrix
! Integral = output of 2e-integral
! KE,KE_t = Kinetic energy(primitive)/total
! M_2e() = Total matrix regarding to 2e-integrals
! max_it = maximum iteration
! natom = number of atom
! nbas = number of basis
! P_new() = New Density Matrix
! P_old() = Old Density Matrix
! PE,PE_t = Potential energy(primitive)/total
! R(i,j) = Distance between the nuclei "i" and "j"
! S,S_M = Overlap Integral : each primitive/Matrix
! step = the number of iterations
! title = title of the basis
! U() = Unitary Matrix
! X() = Matrix Transformation
! xp,yp,zp = Coordinates of P
! xq,yq,zq = Coordinates of Q
! Z(i) = Nuclear charge for atom i
! Zeta(i) = Slater exponent for atom i
!
!====================================================

!====================================================
!             SUBROUTINES(ALPHABETICAL ORDER)
!----------------------------------------------------
!

! DGEMM('N','N',M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) = Lapack Subroutine
! Calculates C = alpha*A*B + beta*C

! Distance(x1,y1,z1,x2,y2,z2,d) = Calculates the distance between atom 1 and 2

! DSYEV('V', 'U', natom, U, natom, Eigen_value, WORK, LWORK, INFO) = Lapack Subroutine
! Calculates the eigenvalues and eigenvectors of U

! F0(t,F_0) = A function related to error function to calculate integrals

! Input(natom,nbas) = Reads number of atoms and basis

! Integral_2e(a,b,c,d,r1,r2,r3,Integral) = Calculates 2e integrals using 4 exponents: a,b,c,d and 3 distances:
! r1 = R(A,B), r2 = R(C,D), r3 = R(P,Q) see fig below

! Kinetic(a,b,d,KE) = Calculates Kinetic Energy integral using 2 exponents: a,b, and 1 distance:
! r1 = R(A,B) see fig below

! Overlap(a,b,d,S) = Calculates the overlap integral using exponents a,b and the distance d

! Potential(a,b,Z,d1,d2,PE) = Calculates Potential Energy integral using 2 exponents: a,b, and 2 distances:
! d1 = R(A,B), d2 = R(P,C) see fig below


!====================================================
!                  FIGURE OF THE PROBLEM
!----------------------------------------------------
!        A------------P----B
!        C----------Q------D
!
!     A,B,C,D are centers for integration
!====================================================
  IMPLICIT NONE

  Integer :: i, j, k, l, m, n, o, p, max_it, natom, nbas, step, lwork, info

  Double Precision :: alpha1, beta1, Converge, delta, d, d1, d2, d3
  Double Precision :: EN, EN_t, Integral, KE, PE, pi
  Double Precision :: S, xp, yp, zp, xq, yq, zq

  Double Precision,Allocatable :: alpha(:,:), coeff(:,:), KE_t(:,:), R(:,:)
  Double Precision,Allocatable :: S_M(:,:), Z(:), Zeta(:),X(:,:)
  Double Precision,Allocatable :: H_core(:,:), G(:,:), F(:,:), Work(:)
  Double Precision,Allocatable :: F_Prime(:,:), P_new(:,:), P_old(:,:)
  Double Precision,Allocatable :: M_2e(:,:,:,:), PE_t(:,:),Eigen_value(:)
  Double Precision,Allocatable :: Eigen_C(:,:), Eigen_C_Prime(:,:),U(:,:)

  Character(70) :: title,atomname

  pi = 3.1415926535898d0

  Converge = 1.0d-4
  max_it = 20
  alpha1 = 1.0d0
  beta1 = 1.0d0


!====================================================
!                  MAIN PROGRAM
!====================================================

	call Input(natom,nbas)

	lwork = 3*natom-1

	Allocate(S_M(natom,natom))
	Allocate(U(natom,natom))
	Allocate(KE_t(natom,natom))
	Allocate(PE_t(natom,natom))
	Allocate(H_core(natom,natom))
	Allocate(P_new(natom,natom))
	Allocate(P_old(natom,natom))
	Allocate(F_Prime(natom,natom))
	Allocate(F(natom,natom))
	Allocate(G(natom,natom))
	Allocate(X(natom,natom))
	Allocate(Eigen_C_Prime(natom,natom))
	Allocate(Eigen_C(natom,natom))
	Allocate(M_2e(natom,natom,natom,natom))
	Allocate(Z(natom))
	Allocate(Eigen_value(natom))
	Allocate(Work(lwork))
	Allocate(Zeta(natom))
	Allocate(R(natom,3))
	Allocate(alpha(natom,nbas))
	Allocate(coeff(natom,nbas))

!====================================================
!                  INPUT BASIS
!====================================================



!----------------Reading basis from file----------------

	open(Unit = 3,File ="Basis.txt")
	read(3,*)
	read(3,"(A)") title

	do i=1,natom
	  read(3,*) atomname,Z(i),Zeta(i)
	  do j=1,nbas
	    read (3,*) coeff(i,j),alpha(i,j)
	  end do
	end do

	close(3)


!====================================================
!                  INPUT GEOMETRY
!====================================================


!----------------Reading coordinates from file----------------

	open(Unit = 4,File ="Input.xyz")
	do i=1,natom
	  read(4,*) atomname,R(i,1),R(i,2),R(i,3)
	end do

	close(4)


!====================================================
!         SCALING of Coefficients & Exponents
!----------------------------------------------------

	do i=1,natom
	  do j=1,nbas
	    alpha(i,j) = alpha(i,j)*(zeta(i)**2)
	    coeff(i,j) = coeff(i,j)*((2.0d0*alpha(i,j)/pi)**0.75d0)
	  end do
	end do

!====================================================
!                  OVERLAP INTEGRAL
!----------------------------------------------------

	open(Unit = 2,File ="Output.txt" )
	write(2,*) title
	write(2,*) "Overlap Matrix S"
!-------------loop over atoms------------------------
	do i=1,natom
	  do j=1,natom
	    S_M(i,j) = 0.0d0
!-------------loop over primitives-------------------
	    do k=1,nbas
	      do l=1,nbas
		call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d)
	        call Overlap(alpha(i,k),alpha(j,l),d,S)
		S_M(i,j) = S_M(i,j) + S*coeff(i,k)*coeff(j,l)
	      end do
	    end do
!----------------------------------------------------

	    write(2,"(f16.10)",advance = "No") S_M(i,j)
	  end do
	write(2,*)
	end do

!====================================================
!                  KINETIC ENERGY
!----------------------------------------------------
	write(2,*) "Kinetic Matrix"
!-------------loop over atoms------------------------
	do i=1,natom
	  do j=1,natom
	    KE_t(i,j) = 0.0d0
!-------------loop over primitives-------------------
	    do k=1,nbas
	      do l=1,nbas
		call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d)
	        call Kinetic(alpha(i,k),alpha(j,l),d,KE)
		KE_t(i,j) = KE_t(i,j) + KE*coeff(i,k)*coeff(j,l)
	      end do
	    end do
!----------------------------------------------------
	    H_core(i,j) = H_core(i,j) + KE_t(i,j)
	    write(2,"(f16.10)",advance = "No") KE_t(i,j)
	  end do
	write(2,*)
	end do

!====================================================
!                  POTENTIAL ENERGY
!----------------------------------------------------
	write(2,*) "Potential Matrix"
!----------------------------------------------------
!                (A|sum(Zc/(r-R))|B)
!----------------------------------------------------

!-------------loop over atoms------------------------

	do i=1,natom
	  do j=1,natom
	    PE_t(i,j) = 0.0d0
!-------------loop over primitives-------------------

	    do l=1,nbas
	      do m=1,nbas

!-----------------Calculating the P Coordinates------

		xp = (alpha(i,l)*R(i,1)+alpha(j,m)*R(j,1) )/(alpha(i,l)+alpha(j,m))
		yp = (alpha(i,l)*R(i,2)+alpha(j,m)*R(j,2) )/(alpha(i,l)+alpha(j,m))
		zp = (alpha(i,l)*R(i,3)+alpha(j,m)*R(j,3) )/(alpha(i,l)+alpha(j,m))

!----------------loop over all nuclei(for Zc)--------

	        do k=1,natom

!----------------------------------------------------
!    d1 is the distance between P and C
!    d is the distance between A and B
!----------------------------------------------------

		  call Distance(xp,yp,zp,R(k,1),R(k,2),R(k,3),d1)
		  call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d)

	          call Potential(alpha(i,l),alpha(j,m),Z(k),d,d1,PE)
		  PE_t(i,j) = PE_t(i,j) + PE*coeff(i,l)*coeff(j,m)

	        end do

	      end do
	    end do
!----------------------------------------------------
	    H_core(i,j) = H_core(i,j) + PE_t(i,j)
	    write(2,"(f16.10)",advance = "No") PE_t(i,j)
	  end do
	  write(2,*)
	end do

!====================================================
!                 2e-Integrals
!----------------------------------------------------


  do i=1,natom
    do j=1,natom
      do k=1,natom
        do l=1,natom

!-------------loop over primitives-------------------
          do m=1,nbas
            do n=1,nbas
	      do o=1,nbas
	        do p=1,nbas

!-----------------Calculating the P Coordinates------

		  xp = (alpha(i,m)*R(i,1)+alpha(j,n)*R(j,1) )/(alpha(i,m)+alpha(j,n))
		  yp = (alpha(i,m)*R(i,2)+alpha(j,n)*R(j,2) )/(alpha(i,m)+alpha(j,n))
		  zp = (alpha(i,m)*R(i,3)+alpha(j,n)*R(j,3) )/(alpha(i,m)+alpha(j,n))

!-----------------Calculating the Q Coordinates------

		  xq = (alpha(k,o)*R(k,1)+alpha(l,p)*R(l,1) )/(alpha(k,o)+alpha(l,p))
		  yq = (alpha(k,o)*R(k,2)+alpha(l,p)*R(l,2) )/(alpha(k,o)+alpha(l,p))
		  zq = (alpha(k,o)*R(k,3)+alpha(l,p)*R(l,3) )/(alpha(k,o)+alpha(l,p))

!--------------------------------------------------------------------------------------
!      d1 = distance(A,B) , d2 = distance(C,D) , d3 = distance(P,Q)
!--------------------------------------------------------------------------------------

		  call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d1)
		  call Distance(R(k,1),R(k,2),R(k,3),R(l,1),R(l,2),R(l,3),d2)
		  call Distance(xp,yp,zp,xq,yq,zq,d3)

		  call Integral_2e(alpha(i,m),alpha(j,n),alpha(k,o),alpha(l,p),d1,d2,d3,Integral)
		  M_2e(i,j,k,l) = M_2e(i,j,k,l) + Integral*coeff(i,m)*coeff(j,n)*coeff(k,o)*coeff(l,p)
                end do
	      end do
	    end do
	  end do
        end do
      end do
    end do
  end do

!===========================================================================
!                   ASSIGNING THE WHOLE 2e-MATRIX
!---------------------------------------------------------------------------

	write(2,*) "All 2e Matrix"
	do i=1,natom
	  do j=1,natom
	    do k=1,natom
	      do l=1,natom
	        write(2, "(A1,I1,I1,A1,I1,I1,A1,f16.10)") '(',i,j,'|',k,l,')',M_2e(i,j,k,l)
	      end do
	    end do
	  end do
	end do

	write(2,*) "**********************************************"


!---------------------Core Hamiltonian-----------------

	write(2,*) "H_core Matrix"
	do i=1,natom
	  do j=1,natom
	    write(2, "(f16.10)", advance = "No") H_core(i,j)
	  end do
	write(2,*)
	end do


!=====================================================
!             RHF PROGRAM (CLOSED SHELL)
!=====================================================

!=====================================================
!        ORTHOGONALIZATION OF THE BASIS(CANONICAL)
!        X = Us^-0.5, where s = diagonalized S
!-----------------------------------------------------

!-------------Diagonalizing S and producing s^-0.5----

	do i=1,natom
	  do j=1,natom
	    U(i,j) = S_M(i,j)
	  end do
	end do

!-------------Calculating eigenvalues & eigenvectors ----

	call DSYEV('V', 'U', natom, U, natom, Eigen_value, WORK, LWORK, INFO)

	do i=1,natom
	  Eigen_value(i) = (Eigen_value(i))**(-0.5d0)
	end do

!-----------------Transformation----------------------

	do i=1,natom
	  do j=1,natom
	    X(i,j) = Eigen_value(j)*U(i,j)
	  end do
	end do

	write(2,*) "X Matrix"
	do i=1,natom
	  do j=1,natom
	    write(2, "(f16.10)", advance = "No") X(i,j)
	  end do
	write(2,*)
	end do

!====================================================
!                 SCF ITERATION
!====================================================

	step = 0

	write(2,*) "Initial Density Matrix"
	do i=1,natom
	  do j=1,natom
	    P_old(i,j) = 0.0d0
	    P_new(i,j) = 0.0d0
	    write(2, "(f16.10)", advance = "No") P_old(i,j)
	  end do
	write(2,*)
	end do

	write(2,*) "**********************************************"

	20 continue


	write(2,*) "Step", step
	step = step+1

!---------------------Form the 2e part of Fock Matrix-----------

	write(2,*) "G Matrix"
	do i=1,natom
	  do j=1,natom
	    G(i,j) = 0.0d0
	    do k=1,natom
	      do l=1,natom
	        G(i,j) = G(i,j)+P_new(k,l)*(M_2e(i,j,k,l)-0.5d0*M_2e(i,l,k,j))
	      end do
	    end do
	    write(2, "(f16.10)", advance = "No") G(i,j)
	  end do
	write(2,*)
	end do

!------------------Constructing Fock Matrix--------------

	write(2,*) "Fock Matrix"
	do i=1,natom
	  do j=1,natom
	    F(i,j) = H_core(i,j)+G(i,j)
	    write(2, "(f16.10)", advance = "No") F(i,j)
	  end do
	write(2,*)
	end do

	write(2,*)

!-----------------------Electronic Energy----------------

	EN = 0.0d0
	do i=1,natom
	  do j=1,natom
	    EN = EN + 0.5d0*P_new(i,j)*(H_core(i,j)+F(i,j))
	  end do
	end do
	write(2,*) "Electronic Energy = ",EN
	write(2,*)

!----------------------Transform and Diagonalize---------

!----------------------G = F*X----------------------------

	G = 0.0d0
        call DGEMM('N','N',natom,natom,natom,alpha1,F,natom,X,natom,beta1,G,natom)

!----------------------F_Prime = XT*G = XT*F*X--------------------

	F_Prime = 0.0d0
        call DGEMM('T','N',natom,natom,natom,alpha1,X,natom,G,natom,beta1,F_Prime,natom)

	do i=1,natom
	  do j=1,natom
	    Eigen_C_prime(i,j) = F_Prime(i,j)
	  end do
	end do

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

	call DSYEV('V', 'U', natom, Eigen_C_prime, natom, Eigen_value, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------

	Eigen_C = 0.0d0
	call DGEMM('N','N',natom,natom,natom,alpha1,X,natom,Eigen_C_prime,natom,beta1,Eigen_C,natom)

!------------------Creating New Density Matrix----------

	P_old = P_new

	do i=1,natom
	  do j=1,natom
	    P_new(i,j) = 0.0d0
	    do k=1,natom/2
	      P_new(i,j) = P_new(i,j) + 2.0d0*Eigen_C(i,k)*Eigen_C(j,k)
	    end do
	  end do
	end do

!---------------------Writing F', C', Energy matrix, C, P--------------------

	write(2,*) "Transformed Fock Matrix"
	do i=1,natom
	  do j=1,natom
	    write(2, "(f16.10)", advance = "No") F_Prime(i,j)
	  end do
	write(2,*)
	end do

	write(2,*) "Transformed C Matrix"
	do i=1,natom
	  do j=1,natom
	    write(2, "(f16.10)", advance = "No") Eigen_C_Prime(i,j)
	  end do
	write(2,*)
	end do

	write(2,*) "E Matrix"
	do i=1,natom
	  do j=1,natom
	    if (i.eq.j) then
	      write(2, "(f16.10)", advance = "No") Eigen_value(i)
	    else
	      write(2, "(f16.10)", advance = "No") 0.0d0
	    end if
	  end do
	write(2,*)
	end do

	write(2,*) "C Matrix"
	do i=1,natom
	  do j=1,natom
	    write(2, "(f16.10)", advance = "No") Eigen_C(i,j)
	  end do
	write(2,*)
	end do

	write(2,*) "Density Matrix"
	do i=1,natom
	  do j=1,natom
	    write(2, "(f16.10)", advance = "No") P_new(i,j)
	  end do
	write(2,*)
	end do

!---------------------Convergence of Density Matrix(Using Standard Deviation)---------

	delta = 0.0d0
	do i=1,natom
	  do j=1,natom
	    delta = delta+(P_new(i,j)-P_old(i,j))**2
	  end do
	end do

	delta = dsqrt(delta/4.0d0)
	write(2,*) "Delta = ",delta

	write(2,*) "*********************************************"

!-----------------Check for Convergence-------------

	if (delta .lt. Converge) then
	  EN_t = EN
	  do i=1,natom-1
	    do j=i+1,natom
	      call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d)
	      EN_t = EN_t + Z(i)*Z(j)/d
	    end do
	  end do

	  write(2,*) "Congratulations! Calculation Converged"
	  write(2,*) "Electronic Energy = ",EN
	  write(2,*) "Total Energy = ",EN_t

	else if(step .lt. max_it) then
	  goto 20
	else
	  write(2,*) "It does not converge in 20 steps!"
	end if


  End Program

!====================================================
!                 END OF THE MAIN PROGRAM
!====================================================

!====================================================
!                   SUBROUTINES
!----------------------------------------------------

!====================================================
!                      INPUT
!====================================================

  Subroutine Input(natom,nbas)
	Integer :: nbas,natom


	write(*,*)
	write(*,*) "|--------------Restricted HF Calculations-------------|"
	write(*,*)
	write(*,*) "Please specify the number of atoms:"
	read(*,*) natom
	write(*,*) "Please specify the number of Gaussian primitive for each atom:"
	read(*,*) nbas

  return
  end



!====================================================
!                  DISTANCE
!====================================================

  Subroutine Distance(x1,y1,z1,x2,y2,z2,d)
	Double Precision :: d,x1,y1,z1,x2,y2,z2



!-------------Calculating distance between points "1" and "2"------------

	d = dsqrt( (x1-x2) **2.0d0+(y1-y2) **2.0d0+(z1-z2) **2.0d0)

  return
  end

!====================================================
!                  OVERAP INTEGRAL (UNNORMALIZED)
!====================================================

  Subroutine Overlap(a,b,d,S)
	Double Precision :: a,b,d,pi,S
	pi = 3.1415926535898d0

	S = ( ( pi/(a+b) )**1.5d0)*exp((-a*b/(a+b))*d**2.0d0)
  return
  end
!====================================================
!                  KINETIC PART (UNNORMALIZED)
!====================================================

  Subroutine Kinetic(a,b,d,KE)
	Double Precision :: a,b,d,KE,pi
	pi = 3.1415926535898d0

	KE = a*b/(a+b)*(3.0d0-2.0d0*a*b/(a+b)*d**2)*((pi/(a+b))**1.5d0)*exp(-a*b/(a+b)*d**2.0d0)
  return
  end

!====================================================
!                  F0 FUNCTION
!====================================================

  Subroutine F0(t,F_0)
  	Double Precision :: pi,t,F_0
	pi = 3.1415926535898d0

!----------prevent from dividing by zero-----------

	if (t .gt. 1.0d-6) then
	  F_0 = 0.5d0*(pi/t)**0.5d0 *erf(t**0.5d0)
	else
	  F_0 = 1.0d0-t/3.0d0
	end if
  return
  end

!====================================================
!                  POTENTIAL PART (UNNORMALIZED)
!====================================================

  Subroutine Potential(a,b,Z,d1,d2,PE)
	Double Precision :: a,b,d1,d2,F_0,PE,pi,Z
	pi = 3.1415926535898d0

	call F0((a+b)*d2**2.0d0,F_0)
	PE = -2.0d0*pi/(a+b)*Z*exp(-a*b/(a+b)*d1**2.0d0)*F_0
  return
  end
!====================================================
!              TWO ELECTRON INTEGRAL (UNNORMALIZED)
!====================================================

  Subroutine Integral_2e(a,b,c,d,r1,r2,r3,Integral)
	Double Precision :: a,b,c,d,F_0,Integral,pi,r1,r2,r3
	pi = 3.1415926535898d0

	call F0( (a+b)*(c+d)/(a+b+c+d)*r3**2.0d0,F_0 )
	Integral = 2.0d0*pi**2.5d0/((a+b)*(c+d)*(a+b+c+d)**0.5d0)*exp(-a*b/(a+b)*r1**2.0d0 - c*d/(c+d)*r2**2.0d0)*F_0
  return
  end



!====================================================
!                  END OF THE PROGRAM
!====================================================
