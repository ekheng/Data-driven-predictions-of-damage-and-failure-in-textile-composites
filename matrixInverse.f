      subroutine matrixInverse(a,y,n,np)
      
      !taken from http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
      !a goes in, y=inverse(a) goes out
      
      implicit none
      
      INTEGER np,indx(np),n,flag
      REAL*8 a(np,np),y(np,np),d
      integer i,j

      do i=1,n !Set up identity matrix.
      do j=1,n
      y(i,j)=0.
      enddo 
      y(i,i)=1.
      enddo 
      call ludcmp(a,n,np,indx,d,flag) !Decompose the matrix just once.
      do  j=1,n                !Find inverse by columns.
      call lubksb(a,n,np,indx,y(1,j))
      !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
      !address of the jth column of y.
      enddo 
      
      return
      end

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------  
! same as above, quad precision      
      
      subroutine qmatrixInverse(a,y,n,np)

      !taken from http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
      !a goes in, y=inverse(a) goes out
      
      implicit none
      
      INTEGER np,indx(np),n,flag
      REAL*16 a(np,np),y(np,np),d
      integer i,j

      do i=1,n !Set up identity matrix.
      do j=1,n
      y(i,j)=0.
      enddo 
      y(i,i)=1.
      enddo 
      call ludcmp(a,n,np,indx,d,flag) !Decompose the matrix just once.
      do  j=1,n                !Find inverse by columns.
      call lubksb(a,n,np,indx,y(1,j))
      !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
      !address of the jth column of y.
      enddo 
      
      return
      end

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------      
 
      subroutine solveLinSysLU(A,b,n,flag)
      !solves Ax=b via LU-decomposition 
      !double precision
      !input: A,b
      !output: x in place of b
      !        successfull: flag=1, not s.: flag=0

      implicit none
      
      integer::n
      integer::indx(n)
      integer::flag
      real*8::A(n,n),b(n),d
      
!      write(*,*) 'A='
!      write(*,*) A
!      write(*,*) 'b='
!      write(*,*) b
!      write(*,*) '---'
      
 
      call ludcmp(A,n,n,indx,d,flag)
      if (flag.eq.1) then
      call lubksb(A,n,n,indx,b)
      endif
            
      return
      end     

      
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------      
 
      subroutine qSolveLinSysLU(A,b,n,flag)
      !solves Ax=b via LU-decomposition 
      !quad precision
      !input: A,b
      !output: x in place of b
      !        successfull: flag=1, not s.: flag=0

      implicit none
      
      integer::n
      integer::indx(n)
      integer::flag
      real*16::A(n,n),b(n),d
      
      call qludcmp(A,n,n,indx,d,flag)
      call qlubksb(A,n,n,indx,b)
            
      return
      end    
      
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------

      SUBROUTINE ludcmp(a,n,np,indx,d,flag)
      
      implicit none
      
      INTEGER n,np,indx(n),NMAX,flag
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=10,TINY=1.0e-20) !Largest expected n, and a small number.
      !Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
      !the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
      !arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
      !row permutation effected by the partial pivoting; d is output as ±1 depending on whether
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
      do j=1,n !This is the loop over columns of Crout’s method.
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
      
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------     
!same as above with quad precision

      SUBROUTINE qludcmp(a,n,np,indx,d,flag)
      
      implicit none
      
      INTEGER n,np,indx(n),NMAX,flag
      REAL*16 d,a(np,np),TINY
      PARAMETER (NMAX=10,TINY=1.0e-20) !Largest expected n, and a small number.
      !Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
      !the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
      !arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
      !row permutation effected by the partial pivoting; d is output as ±1 depending on whether
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
      do j=1,n !This is the loop over columns of Crout’s method.
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

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------

      SUBROUTINE lubksb(a,n,np,indx,b)
      
      implicit none
      
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      !Solves the set of n linear equations A · X = B. Here a is input, not as the matrix A but
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
      
      
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
!same as above with quad precision

      SUBROUTINE qlubksb(a,n,np,indx,b)
      
      implicit none
      
      INTEGER n,np,indx(n)
      REAL*16 a(np,np),b(n)
      !Solves the set of n linear equations A · X = B. Here a is input, not as the matrix A but
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
      