!C For: Pure Matrix - nonlinear behavior with crackband (isotropic damage/failure)
!C Author: Royan Dmello, Aerospace Engineering, University of Michigan
!C Updated: 06/20/2022
      SUBROUTINE VUMAT_MATRIX(STATEV,DDSDDE,
     2 CMNAME,
     3 NDI, NSHR,NTENS,NSTATV,PROPS,NPROPS,
     4 CELENT, SIG, TIMEVUMAT)
     
      IMPLICIT NONE

      CHARACTER*80 CMNAME
      INTEGER:: NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NTENS,NSTATV,NPROPS
      INTEGER:: NDI,NSHR
      REAL*8:: STRESS(6), DSTRAN(6), STRAN(6), DDSDDT(6), DRPLDE(6)
      REAL*8:: DDSDDE(6,6)
      REAL*8:: STATEV(NSTATV), PROPS(NPROPS), TIME(2)
      REAL*8:: DTIME,TEMP,DTEMP,CELENT,PNEWDT,DRPLDT
      REAL*8:: SSE, SPD, SCD,RPL
      REAL*8:: PREDEF(1), DPRED(1)
      REAL*8:: COORDS(3,3),DROT(3,3), DFGRD0(3,3),DFGRD1(3,3)
      REAL*8:: STRAIN(NTENS), SIG(NTENS)
      REAL*8:: SS(6,6), CC(6,6)
      REAL*8:: DMLIM, Em, Em0, vm, vm0, SIG_Y, k1, k2, EPS_Y
	  REAL*8:: sigcr, GIc, Em_s_prev, vm_s_prev, PSMAX, PEMAX
	  REAL*8:: eps_m(6), Gm, Ecr0, eps0m, epsfm, D, h, BZ
	  REAL*8:: STRESS_OLD(6), PSTRESS(3), PSTRAIN(3), T(3,3), TT(3,3)
	  REAL*8:: EPS_J2, SIG_J2, Em_s, vm_s, TIMEVUMAT, del_eps_max

 
 
	  Em0 = PROPS(1)	  
      vm0 = PROPS(2)	  
      SIG_Y = PROPS(3)            ! SIG_Y = 15.66d0 yield stress
      k1 = PROPS(4)               ! k1 = 3308.0d0 hardening parameter
      k2 = PROPS(5)               !k2 = 56.92d0 hardening parameter
      EPS_Y = SIG_Y/Em0           ! Yield strain in matrix
      sigcr = PROPS(6)  
      GIc = PROPS(7)  	  

              
      DMLIM = 0.02d0 
 
	  IF (STATEV(16).lt.0.5d0) THEN
         STATEV(9) = PROPS(1)            !Em: Young's modulus
         STATEV(10) = PROPS(2)           !vm: Poisson's ratio
         STATEV(33) = NOEL
         STATEV(21) = 0.0d0
         STATEV(22) = 0.0d0
         STATEV(23) = 0.0d0
         STATEV(24) = 0.0d0
         STATEV(25) = 0.0d0
         STATEV(26) = 0.0d0
         D = 1.0d0
      END IF	  

!C----------------------------------------------------------
!C	  LIST AND DESCRIPTION OF STATE VARIABLES USED
!C----------------------------------------------------------
     
      !Store Stress from the previous increment      
      STRESS(1) = STATEV(21)
      STRESS(2) = STATEV(22)
      STRESS(3) = STATEV(23)
      STRESS(4) = STATEV(24)
      STRESS(5) = STATEV(25)
      STRESS(6) = STATEV(26)
	  

      STRESS_OLD = STRESS

      !Applied strains during the loading step
      STRAIN(1) = STATEV(27)
      STRAIN(2) = STATEV(28)
      STRAIN(3) = STATEV(29)
      STRAIN(4) = STATEV(30)
      STRAIN(5) = STATEV(31)
      STRAIN(6) = STATEV(32)
	  
      eps_m(1) = STRAIN(1)
	  eps_m(2) = STRAIN(2)
	  eps_m(3) = STRAIN(3)
	  eps_m(4) = STRAIN(4)
	  eps_m(5) = STRAIN(5)
	  eps_m(6) = STRAIN(6)
	  	  
!C1--------------------------------------------------------------------
      IF (STATEV(16).gt.0.5d0) THEN  !2nd increment onwards
		  CALL mysprind(STRESS_OLD, PSTRESS, T, 1, 3, 3)
	      CALL mysprind(eps_m, PSTRAIN, TT, 1, 3, 3)
          PSMAX = dmax1((PSTRESS(1)),(PSTRESS(2)),(PSTRESS(3)))
          PEMAX = dmax1((PSTRAIN(1)),(PSTRAIN(2)),(PSTRAIN(3)))	
          STATEV(17) = PSMAX
          STATEV(18) = PEMAX
      END IF	

      IF ((STATEV(16).gt.0.5d0).and.(STATEV(1).lt.0.5d0)) THEN 	  
	     !Matrix secant modulus and secant poisson's ratio from previous inc   
		 Em_s_prev = STATEV(9)  
	     vm_s_prev = STATEV(10)
         STATEV(7) = Em_s_prev
	     STATEV(8) = vm_s_prev
      

          !Calculate equivalent matrix strain from matrix constituent strains
          CALL calc_eqstrain_mm(eps_m, EPS_J2)     
          !Compute equivalent stress
          IF (EPS_J2.LE.EPS_Y) THEN
             SIG_J2 = Em0*EPS_J2      !Secant will be equal to initial modulus
          ELSEIF (EPS_J2.GT.EPS_Y) THEN
             STATEV(13) = 1.0d0
             SIG_J2 = SIG_Y - (k1/k2)*(DEXP(-k2*EPS_J2) - 
     1                           DEXP(-k2*SIG_Y/Em0))
	      END IF
	           
          !Compute secant modulus and secant Poisson's ratio         
		  IF (dabs(EPS_J2).gt.0.0d0) THEN
		  Em_s = DABS(SIG_J2/EPS_J2)
		  ELSE
		  Em_s = STATEV(9)
		  END IF         
		  vm_s = 0.5d0 + (Em_s/Em0)*(vm0 - 0.5d0)
	    
          STATEV(9) = Em_s
          STATEV(10) = vm_s
          STATEV(11) = EPS_J2
          STATEV(12) = SIG_J2

      END IF
!C1--------------------------------------------------------------------	  


!C!C---------------------------------------------------------
!C!C  WHEN THE MATERIAL FIRST ENTERS CRACKBAND
!C!C---------------------------------------------------------
      IF (STATEV(16).gt.0.5d0) THEN  !2nd increment onwards

		IF ((PSMAX.ge.sigcr).and.(STATEV(1).lt.0.5D0)) THEN
           STATEV(3) = PEMAX    !This will be used for max eps yet in CB
           STATEV(1) = 1.0D0
           STATEV(4) = PSMAX/PEMAX     
           eps0m = sigcr/STATEV(4)
           STATEV(5) = eps0m
 
         
      ! Crackband parameters stored once at crack initiation
      h = CELENT
      STATEV(15) = CELENT 
      epsfm = (2.0D0*GIc)/(h*sigcr)
      STATEV(6) = epsfm  
               
      ! Check Bazant limit
      BZ = 2.0D0*STATEV(4)*GIc/(sigcr**2)  
         IF (h.ge.0.95*BZ) THEN
           write(*,*) 'Char. element length is too large!!'
           write(*,*) 'Reduce mesh size'
           write(*,*) 'h', h
           write(*,*) 'BZ', BZ
           write(*,*) 'GIc', GIc
           write(*,*) 'sigcr', sigcr		
           CALL XPLB_EXIT
         END IF
	    END IF
      END IF

      
      
      
!C!C-----------------------------------------------------------	  
!C!C-----------------------------------------------------------
!C!C STIFFNESS DEGRADATION FOR CRACKED/DAMAGED MATRIX
!C!C-----------------------------------------------------------
!C!C-----------------------------------------------------------

      IF (STATEV(16).gt.0.5d0) THEN    !Only 2nd increment onwards
        IF (STATEV(1).lt.0.5D0) THEN   !If not yet in crackband
            D = 1.0D0
            STATEV(2) = D
	    ELSEIF ((STATEV(1).gt.0.5D0).and.(STATEV(1).lt.1.5D0)
     1                              .and.(PEMAX.ge.STATEV(3))) THEN
		    STATEV(3) = dabs(PEMAX)
            Ecr0 = STATEV(4)
            eps0m = STATEV(5)
            epsfm = STATEV(6)
            D=(sigcr/Ecr0)*(1.0D0/(epsfm-eps0m))*((epsfm/PEMAX) -1.0d0)
			STATEV(2) = D
            IF (D.lt.DMLIM) THEN
			   STATEV(2) = DMLIM
               STATEV(1) = 2.0D0
               STATEV(20) = 0.0D0   ! Delete element
            ELSE IF (STATEV(2).gt.1.0D0) THEN
               STATEV(2) = 1.0D0
            END IF
	    ELSEIF ((STATEV(1).gt.0.5D0).and.(STATEV(1).lt.1.5D0)
     1                              .and.(PEMAX.lt.STATEV(3))) THEN
               D = STATEV(2)
        ELSEIF (STATEV(1).gt.1.5d0) THEN
               D = DMLIM
               STATEV(2) = DMLIM        
        END IF
		
		
			
		IF (STATEV(2).lt.DMLIM) THEN
           STATEV(2) = DMLIM
           STATEV(1) = 2.0D0
           STATEV(20) = 0.0d0
        END IF
        D = STATEV(2)


      END IF
      
  	  	  
!C!-------------------------------------------------------
!C!      !Compute stiffness matrix 
!C!-------------------------------------------------------
	  
      Em = D*STATEV(9)
      vm = STATEV(10)
      Gm = 0.5d0*Em*(1.0d0 + vm)
	  	  	  
	  !ASSEMBLE THE STIFFNESS MATRIX	  
      SS = 0.0D0
      SS(1,1) = 1.0D0/Em
      SS(2,2) = 1.0D0/Em
      SS(3,3) = 1.0D0/Em
      SS(4,4) = 1.0D0/Gm
      SS(5,5) = 1.0D0/Gm
      SS(6,6) = 1.0D0/Gm
      SS(2,1) = -vm/Em
      SS(1,2) = SS(2,1)
      SS(3,1) = -vm/Em
      SS(1,3) = SS(3,1)
      SS(3,2) = -vm/Em
      SS(2,3) = SS(3,2)
      
      CALL m_matrixInverse(SS, CC, 6, 6)   
	  
!C-----------------------------------------------------------
!C STRESS UPDATE
!C-----------------------------------------------------------

      STRESS = MATMUL(CC, STRAIN)
	               
      !Store stress at the end of the increment	   
      STATEV(21) = STRESS(1)
      STATEV(22) = STRESS(2)
      STATEV(23) = STRESS(3)
      STATEV(24) = STRESS(4)
      STATEV(25) = STRESS(5)
      STATEV(26) = STRESS(6)
	  
      SIG = STRESS	  
	  STATEV(16) = 1.0d0
	       
	  RETURN
      END







      
      
      




!C STATE VARIABLES list for pure matrix
!C      1		crackflag in matrix
!C      2		D
!C      3		eps_max greatest in CB yet
!C      4		Ecr0
!C      5		eps_0
!C      6		eps_f
!C      7		Em (prev inc.)
!C      8		vm (prev inc.)
!C      9		Em (end of inc.)
!C      10		vm (end of inc.)
!C      11		eps_J in matrix
!C      12		eq_stress in matrix
!C      13		nonlinear flag in matrix
!C      14		Element number
!C      15		Characteristic length
!C      16		Initialization flag for matrix
!C      17		Max principal stress
!C      18		Max principal strain
!C      19		
!C      20		delete element flag
!C      21		Stress1 matrix
!C      22		Stress 2 in matrix
!C      23		Stress 3 in matrrix
!C      24		Stress 4 in matrix
!C      25		Stress 5 in matrix
!C      26		Stress 6 in matrix
!C      27		Strain 1 in matrix
!C      28		Strain 2 in matrix
!C      29		Strain 3 in matrix
!C      30		Strain 4 in matrix
!C      31		Strain 5 in matrix
!C      32		Strain 6 in matrix
!C------------------------------------------
!C INPUTS AND STATE VARIABLES
!C------------------------------------------









 
      











!C!C---------------------------------------------------------------------
!C!C           SUBROUTINES START HERE
!C!C---------------------------------------------------------------------

!C!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_stress_mm(Q_IN, EPS_IN, SIG_OUT)
      REAL*8, INTENT(IN) :: Q_IN(6,6)
      REAL*8, INTENT(IN) :: EPS_IN(6)
      REAL*8, INTENT(OUT) :: SIG_OUT(6)

      SIG_OUT = MATMUL(Q_IN, EPS_IN)

      RETURN
      END


!C!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_eqstrain_mm(EP, j2eps)

      IMPLICIT NONE
      REAL*8, INTENT(IN):: EP(6)
      REAL*8, INTENT(OUT):: j2eps
      REAL*8:: exx, eyy, ezz, gxy, gxz, gyz
  
      exx = (2.0d0/3.0d0)*EP(1) - (1.0d0/3.0d0)*EP(2) 
     1                 - (1.0d0/3.0d0)*EP(3)
      eyy = (2.0d0/3.0d0)*EP(2) - (1.0d0/3.0d0)*EP(1) 
     1                 - (1.0d0/3.0d0)*EP(3)     
      ezz = (2.0d0/3.0d0)*EP(3) - (1.0d0/3.0d0)*EP(1) 
     1                 - (1.0d0/3.0d0)*EP(2)
      gxy = EP(4)
      gxz = EP(5)
      gyz = EP(6)


      j2eps = (2.0d0/3.0d0)*DSQRT( 1.5d0*(exx**2 + eyy**2 + ezz**2) 
     1            + 0.75d0*(gxy**2 + gxz**2 + gyz**2))
           
      RETURN
      END


!C!C-------------------------------------------------------------------
!C!C-------------------------------------------------------------------
!C!C SUBROUTINES TO INVERT THE MATRIX
!C!C-------------------------------------------------------------------
!C!C-------------------------------------------------------------------
      subroutine m_matrixInverse(a,y,n,np)

      !taken from http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
      !a goes in, y=inverse(a) goes out

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER np,indx(np),n,flag
      REAL*8 a(np,np),y(np,np),d
      integer i,j

      do i=1,n !Set up identity matrix.
      do j=1,n
      y(i,j)=0.
      enddo
      y(i,i)=1.
      enddo
      call m_ludcmp(a,n,np,indx,d,flag) !Decompose the matrix just once.
      do  j=1,n                !Find inverse by columns.
      call m_lubksb(a,n,np,indx,y(1,j))
      !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
      !address of the jth column of y.
      enddo

      return
      end

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------
!C! same as above, quad precision

      subroutine m_qmatrixInverse(a,y,n,np)

      !taken from http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
      !a goes in, y=inverse(a) goes out

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER np,indx(np),n,flag
      REAL*16 a(np,np),y(np,np),d
      integer i,j

      do i=1,n !Set up identity matrix.
      do j=1,n
      y(i,j)=0.
      enddo
      y(i,i)=1.
      enddo
      call m_ludcmp(a,n,np,indx,d,flag) !Decompose the matrix just once.
      do  j=1,n                !Find inverse by columns.
      call m_lubksb(a,n,np,indx,y(1,j))
      !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
      !address of the jth column of y.
      enddo

      return
      end

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

      subroutine m_solveLinSysLU(A,b,n,flag)
      !solves Ax=b via LU-decomposition
      !double precision
      !input: A,b
      !output: x in place of b
      !        successfull: flag=1, not s.: flag=0

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      integer::n
      integer::indx(n)
      integer::flag
      real*8::A(n,n),b(n),d

!      write(*,*) 'A='
!      write(*,*) A
!      write(*,*) 'b='
!      write(*,*) b
!      write(*,*) '---'


      call m_ludcmp(A,n,n,indx,d,flag)
      if (flag.eq.1) then
      call m_lubksb(A,n,n,indx,b)
      endif

      return
      end


!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

      subroutine m_qSolveLinSysLU(A,b,n,flag)
      !solves Ax=b via LU-decomposition
      !quad precision
      !input: A,b
      !output: x in place of b
      !        successfull: flag=1, not s.: flag=0

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      integer::n
      integer::indx(n)
      integer::flag
      real*16::A(n,n),b(n),d

      call m_qludcmp(A,n,n,indx,d,flag)
      call m_qlubksb(A,n,n,indx,b)

      return
      end

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

      SUBROUTINE m_ludcmp(a,n,np,indx,d,flag)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n),NMAX,flag
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=10,TINY=1.0e-20) !Largest expected n, and a small number.
      !Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
      !the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
      !arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
      !row permutation effected by the partial pivoting; d is output as �1 depending on whether
      !the number of row interchanges was even or odd, respectively. This routine is used in
      !combination with lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
      flag=1 !So far successfull ...
      d=1. !No row interchanges yet.
      do i=1,n !Loop over rows to get the implicit scaling information
      aamax=0.
      do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
      enddo
      if (aamax.eq.0.) then
        !write(*,*) 'singular matrix in ludcmp' !No nonzero largest element.
        !write(*,*)
        !write(*,*) a
        !CALL myExit()
        flag=0
        exit
      endif
      vv(i)=1./aamax !Save the scaling.
      enddo
      if (flag.eq.1) then
      do j=1,n !This is the loop over columns of Crout�s method.
      do i=1,j-1 !This is equation (2.3.12) except for i = j.
      sum=a(i,j)
      do k=1,i-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      enddo
      aamax=0. !Initialize for the search for largest pivot element.
      do i=j,n !This is i = j of equation (2.3.12) and i = j+1. . .N
      sum=a(i,j) !of equation (2.3.13).
      do k=1,j-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum) !Figure of merit for the pivot.
      if (dum.ge.aamax) then !Is it better than the best so far?
      imax=i
      aamax=dum
      endif
      enddo
      if (j.ne.imax)then !Do we need to interchange rows?
      do k=1,n !Yes, do so...
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
      enddo
      d=-d !...and change the parity of d.
      vv(imax)=vv(j) !Also interchange the scale factor.
      endif
      indx(j)=imax
      if(a(j,j).eq.0.)a(j,j)=TINY
      !If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
      !For some applications on singular matrices, it is desirable to substitute TINY
      !for zero.
      if(j.ne.n)then !Now, finally, divide by the pivot element.
      dum=1./a(j,j)
      do i=j+1,n
      a(i,j)=a(i,j)*dum
      enddo
      endif
      enddo !Go back for the next column in the reduction.
      endif
      return
      END

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------
!C!same as above with quad precision

      SUBROUTINE m_qludcmp(a,n,np,indx,d,flag)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n),NMAX,flag
      REAL*16 d,a(np,np),TINY
      PARAMETER (NMAX=10,TINY=1.0e-20) !Largest expected n, and a small number.
      !Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
      !the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
      !arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
      !row permutation effected by the partial pivoting; d is output as �1 depending on whether
      !the number of row interchanges was even or odd, respectively. This routine is used in
      !combination with lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
      REAL*16 aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
      flag=1 !So far successfull ...
      d=1. !No row interchanges yet.
      do i=1,n !Loop over rows to get the implicit scaling information
      aamax=0.
      do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
      enddo
      if (aamax.eq.0.) then
        !write(*,*) 'singular matrix in ludcmp' !No nonzero largest element.
        !write(*,*)
        !write(*,*) a
        !CALL myExit()
        flag=0
      endif
      vv(i)=1./aamax !Save the scaling.
      enddo
      if (flag.eq.1) then
      do j=1,n !This is the loop over columns of Crout�s method.
      do i=1,j-1 !This is equation (2.3.12) except for i = j.
      sum=a(i,j)
      do k=1,i-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      enddo
      aamax=0. !Initialize for the search for largest pivot element.
      do i=j,n !This is i = j of equation (2.3.12) and i = j+1. . .N
      sum=a(i,j) !of equation (2.3.13).
      do k=1,j-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum) !Figure of merit for the pivot.
      if (dum.ge.aamax) then !Is it better than the best so far?
      imax=i
      aamax=dum
      endif
      enddo
      if (j.ne.imax)then !Do we need to interchange rows?
      do k=1,n !Yes, do so...
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
      enddo
      d=-d !...and change the parity of d.
      vv(imax)=vv(j) !Also interchange the scale factor.
      endif
      indx(j)=imax
      if(a(j,j).eq.0.)a(j,j)=TINY
      !If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
      !For some applications on singular matrices, it is desirable to substitute TINY
      !for zero.
      if(j.ne.n)then !Now, finally, divide by the pivot element.
      dum=1./a(j,j)
      do i=j+1,n
      a(i,j)=a(i,j)*dum
      enddo
      endif
      enddo !Go back for the next column in the reduction.
      endif
      return
      END

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

      SUBROUTINE m_lubksb(a,n,np,indx,b)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      !Solves the set of n linear equations A � X = B. Here a is input, not as the matrix A but
      !rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
      !permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
      !and returns with the solution vector X. a, n, np, and indx are not modified by this routine
      !and can be left in place for successive calls with different right-hand sides b. This routine
      !takes into account the possibility that b will begin with many zero elements, so it is efficient
      !for use in matrix inversion.
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0 !When ii is set to a positive value, it will become the index
      !of the first nonvanishing element of b. We now do
      !the forward substitution, equation (2.3.6). The only new
      !wrinkle is to unscramble the permutation as we go.
      do i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii.ne.0)then
      do j=ii,i-1
      sum=sum-a(i,j)*b(j)
      enddo
      else if (sum.ne.0.) then
      ii=i !A nonzero element was encountered, so from now on we will
      endif !have to do the sums in the loop above.
      b(i)=sum
      enddo
      do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
      sum=b(i)
      do j=i+1,n
      sum=sum-a(i,j)*b(j)
      enddo
      b(i)=sum/a(i,i) !Store a component of the solution vector X.
      enddo
      return !All done!
      END


!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------
!C!same as above with quad precision

      SUBROUTINE m_qlubksb(a,n,np,indx,b)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n)
      REAL*16 a(np,np),b(n)
      !Solves the set of n linear equations A � X = B. Here a is input, not as the matrix A but
      !rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
      !permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
      !and returns with the solution vector X. a, n, np, and indx are not modified by this routine
      !and can be left in place for successive calls with different right-hand sides b. This routine
      !takes into account the possibility that b will begin with many zero elements, so it is efficient
      !for use in matrix inversion.
      INTEGER i,ii,j,ll
      REAL*16 sum
      ii=0 !When ii is set to a positive value, it will become the index
      !of the first nonvanishing element of b. We now do
      !the forward substitution, equation (2.3.6). The only new
      !wrinkle is to unscramble the permutation as we go.
      do i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii.ne.0)then
      do j=ii,i-1
      sum=sum-a(i,j)*b(j)
      enddo
      else if (sum.ne.0.) then
      ii=i !A nonzero element was encountered, so from now on we will
      endif !have to do the sums in the loop above.
      b(i)=sum
      enddo
      do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
      sum=b(i)
      do j=i+1,n
      sum=sum-a(i,j)*b(j)
      enddo
      b(i)=sum/a(i,i) !Store a component of the solution vector X.
      enddo
      return !All done!
      END
