! invert a 3x3 matrix directly

      subroutine inverseMatrix3x3(J,invJ)
      
      implicit none
      
      real*8,dimension(3,3)::J,invJ
      real*8::detJ
      
      detJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
     1    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1))
     2    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
      
      invJ(1,:)=(/ (J(2,2)*J(3,3) - J(2,3)*J(3,2)), 
     1            -(J(1,2)*J(3,3) - J(1,3)*J(3,2)),  
     2             (J(1,2)*J(2,3) - J(1,3)*J(2,2)) /)
      invJ(2,:)=(/-(J(2,1)*J(3,3) - J(2,3)*J(3,1)),  
     1             (J(1,1)*J(3,3) - J(1,3)*J(3,1)), 
     2            -(J(1,1)*J(2,3) - J(1,3)*J(2,1)) /)
      invJ(3,:)=(/ (J(2,1)*J(3,2) - J(2,2)*J(3,1)), 
     1            -(J(1,1)*J(3,2) - J(1,2)*J(3,1)),  
     2             (J(1,1)*J(2,2) - J(1,2)*J(2,1)) /)
      
      invJ=invJ/detJ
      
      return
      end
      
! invert a 3x3 matrix directly, quadruple precision
      subroutine qinverseMatrix3x3(J,invJ)
      
      implicit none
      
      real*16,dimension(3,3)::J,invJ
      real*16::detJ
      
      detJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
     1    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1))
     2    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
      
      invJ(1,:)=(/ (J(2,2)*J(3,3) - J(2,3)*J(3,2)), 
     1            -(J(1,2)*J(3,3) - J(1,3)*J(3,2)),  
     2             (J(1,2)*J(2,3) - J(1,3)*J(2,2)) /)
      invJ(2,:)=(/-(J(2,1)*J(3,3) - J(2,3)*J(3,1)),  
     1             (J(1,1)*J(3,3) - J(1,3)*J(3,1)), 
     2            -(J(1,1)*J(2,3) - J(1,3)*J(2,1)) /)
      invJ(3,:)=(/ (J(2,1)*J(3,2) - J(2,2)*J(3,1)), 
     1            -(J(1,1)*J(3,2) - J(1,2)*J(3,1)),  
     2             (J(1,1)*J(2,2) - J(1,2)*J(2,1)) /)
      
      invJ=invJ/detJ
      
      return
      end
 