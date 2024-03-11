c User subroutine VUMAT
      subroutine vumat(
c Read only -
     1     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew )
c
      include 'vaba_param.inc'
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname

      parameter (
     1     i_umt_nblock = 1,
     2     i_umt_npt    = 2,
     3     i_umt_layer  = 3,
     4     i_umt_kspt   = 4,
     5     i_umt_noel   = 5 )

      call  vumatXtrArg ( jblock(i_umt_nblock),
     1     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
     7     stressNew, stateNew, enerInternNew, enerInelasNew,
     8     jblock(i_umt_noel), jblock(i_umt_npt),
     9     jblock(i_umt_layer), jblock(i_umt_kspt))

      return
      end












! The VUMAT that has information about element number, integration point number ...
c
      subroutine vumatXtrArg (
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
     7     stressNew, stateNew, enerInternNew, enerInelasNew,
     8     nElement, nMatPoint, nLayer, nSecPoint )
C
      include 'vaba_param.inc'
      
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1     charLength(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)

      dimension nElement(nblock)
      character*80 cmname

      real*8::eps(6)
      real*8::sig(6)
      real*8::statev(nstatev)
      real*8::ddsdde(6,6)
	  integer:: i

      do i=1,nblock

        sig(1) = stressOld(i,1)
        sig(2) = stressOld(i,2)
        sig(3) = stressOld(i,3)
        sig(4) = stressOld(i,4)
        sig(5) = stressOld(i,6)
        sig(6) = stressOld(i,5)		
        
        statev=stateOld(i,:)     

             
      ! INITIALIZE VARIABLES  
      IF ((totalTime.le.dt).and.(cmname(1:3).eq.'TOW')) THEN	      		  
		  STATEV(120) = 1.0d0
		  STATEV(102) = 1.0d0   !Stiffness degradation factors - 11 direction
		  STATEV(107) = 1.0d0   !                              - 22 direction 
		  STATEV(112) = 1.0d0  !                              - 33 direction
          STATEV(100) = 0.0d0  !Initialization flag in the code		  
      END IF
	  
	  
	  IF ((totalTime.le.dt).and.(cmname(1:3).eq.'MAT')) THEN	      		  
		  STATEV(20) = 1.0d0
		  STATEV(2) = 1.0d0   !Stiffness degradation factors
          STATEV(16) = 0.0d0  !Initialization flag in the code
	  END IF


      IF (cmname(1:3).eq.'TOW') THEN
         statev(151:153) = statev(151:153) + strainInc(i,1:3)
         statev(154) = statev(154) + strainInc(i,4)*2.0d0
         statev(155) = statev(155) + strainInc(i,6)*2.0d0
         statev(156) = statev(156) + strainInc(i,5)*2.0d0

         
		 
         CALL VUMAT_TOW(statev, ddsdde, cmname, ndir, nshr, 6,
     1 nstatev, props, nprops, charLength(i), sig, totalTime)      
      END IF
      IF (cmname(1:3).eq.'MAT') THEN
         statev(27:29) = statev(27:29) + strainInc(i,1:3)
         statev(30) = statev(30) + strainInc(i,4)*2.0d0
         statev(31) = statev(31) + strainInc(i,6)*2.0d0
         statev(32) = statev(32) + strainInc(i,5)*2.0d0
 		 !write(*,*) 'SDV29', statev(29)         
         CALL VUMAT_MATRIX(statev, ddsdde, cmname, ndir, nshr, 6,
     1 nstatev, props, nprops, charLength(i), sig, totalTime)      
      END IF

         stressNew(i,1) = sig(1)
         stressNew(i,2) = sig(2)
         stressNew(i,3) = sig(3)
         stressNew(i,4) = sig(4)
         stressNew(i,5) = sig(6)
         stressNew(i,6) = sig(5)    
         stateNew(i,:) = statev
	  

      enddo

      return
      end      
      
      

      
      INCLUDE 'RJD-matrix-failure.f'
      INCLUDE 'RJD-tow-failure.f'
      INCLUDE 'matrixInverse.f'
      INCLUDE 'myaux.f'            
      INCLUDE 'myExplicitSupport.f'

