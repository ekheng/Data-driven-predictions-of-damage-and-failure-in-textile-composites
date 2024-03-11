!other subroutines needed in Abaqus/Standard when running SCA
!
!These are utility routines that are differently named in Standard
!and Explicit
!

! Finding principle directions 
      subroutine mysprind(S,PS,AN,LSTR,NDI,NSHR)
      implicit none
      integer::NDI  !number of direct components
      integer::NSHR !number of shear components
      integer::LSTR !flag for stress or engineering strain, not 
                    !applicable for explicit simulation
      real*8::S(NDI+NSHR) !tensor that of which to find eigenvalues
      real*8::PS(3)       !principle stresses
      real*8::AN(3,3)     !principle directions
      
      !make blocks for explicit
      real*8::stress(1,NDI+NSHR)
      real*8::eigVal(1,3)
      real*8::eigVec(1,3,3)
      
      stress(1,:)=S
      call vsprind(1,stress,eigVal,eigVec,ndi,nshr)
      PS=eigVal(1,:)
      AN(1,:)=eigVec(1,:,1)
      AN(2,:)=eigVec(1,:,2)
      AN(3,:)=eigVec(1,:,3)      
      return 
      end
      
      
      
      
! Exitting
      subroutine myExit()
      call XPLB_EXIT
      return
      end