c
c User subroutine VUSDFLD for user-defined fields
c
      subroutine vusdfld(
c Read only -
     &   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     &   jElemUid, kIntPt, kLayer, kSecPt, 
     &   stepTime, totalTime, dt, cmname, 
     &   coordMp, direct, T, charLength, props, 
     &   stateOld, 
c Write only -
     &   stateNew, field )
c
      include 'vaba_param.inc'
c
      dimension props(nprops),
     &          jElemUid(nblock), coordMp(nblock, &), 
     &          direct(nblock, 3, 3), T(nblock,3,3), 
     &          stateOld(nblock, nstatev), 
     &          stateNew(nblock, nstatev),
     &          field(nblock, nfieldv)
      character*80 cmname
c
      character*3 PbcData(maxblk*7),cData(maxblk*6)
	  dimension PbjData(maxblk*7),sData(maxblk*6),
	  dimension PbrData(nblock,15),rData(nblock,15)
c	 new variables
       real(8) EpsE, dEpsE, EpsCap, dEpsCap, pressure
c Get PEQCs for that increment
	  jStatus = 1
	  call vgetvrm( 'PEQC', PbrData, PbjData, PbcData, jStatus )
c	   
       if( jStatus .ne. 0 ) then
          call xplb_abqerr(-2,'Utility routine VGETVRM '//
     &      'failed to get variable Pb.',0,zero,' ')
          call xplb_exit
       end if
	  	  
c Get invariant stresses - mises and pressure:
       jStatus = 1
       call vgetvrm( 'S', rData, sData, cData, jStatus )
c	   
       if( jStatus .ne. 0 ) then
          call xplb_abqerr(-2,'Utility routine VGETVRM '//
     &      'failed to get variable S.',0,zero,' ')
          call xplb_exit
       end if
c      
	   do k = 1, nblock			
		  EpsE = PbrData(k,1)
		  stateNew(k,1) = EpsE
		  EpsCap = PbrData(k,2)
		  stateNew(k,2) = EpsCap
c		  
		  dEpsE = EpsE - stateOld(k,1)
		  dEpsCap = EpsCap - stateOld(k,2)
		  pressure = -(rData(k,1)+rData(k,2)+rData(k,3))/3

		  IF(dEpsE .GT. 0 .OR. dEpsCap .GT. 0)THEN  
			  IF(dEpsE .GT. dEpsCap)THEN
				  IF(pressure .GT. 0) THEN
						stateNew(k,3) = 50
				  ELSE
						stateNew(k,3) = 100
				  ENDIF
			  ELSE
				  stateNew(k,3) = 150
			  ENDIF
		  ELSE
			  stateNew(k,3) = 0
		  ENDIF
		  field(k,1) = exp((-1.)*PbrData(k,4))
		  stateNew(k,4) = field(k,1)	
	   end do
c
      return 	  
      end
