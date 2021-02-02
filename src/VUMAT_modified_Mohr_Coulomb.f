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
      character&80 cmname
c
       character*3 PbcData(maxblk*1),cData(maxblk*6)
       dimension PbjData(maxblk*1),jData(maxblk*6)
       dimension PbrData(nblock,1),rData(nblock,15)
c	 new variables
       real(8) EpsE, dEpsE, pressure	 	   
c Get PEQCs for that increment
	  jStatus = 1
	  call vgetvrm( 'PEEQ', PbrData, PbjData, PbcData, jStatus )
c	   
       if( jStatus .ne. 0 ) then
          call xplb_abqerr(-2,'Utility routine VGETVRM '//
     &      'failed to get variable Pb.',0,zero,' ')
          call xplb_exit
       end if
	  	  
c Get invariant stresses - mises and pressure:
       jStatus = 1
       call vgetvrm( 'S', rData, jData, cData, jStatus )
c	   
       if( jStatus .ne. 0 ) then
          call xplb_abqerr(-2,'Utility routine VGETVRM '//
     &      'failed to get variable S.',0,zero,' ')
          call xplb_exit
       end if
c	calculate pressure. Store them as SDV.
	   do k = 1, nblock			
		  EpsE = PbrData(k,1)
		  stateNew(k,1) = EpsE		
		  dEpsE = EpsE - stateOld(k,1)
		  pressure = -(rData(k,1)+rData(k,2)+rData(k,3))/3
		  IF(dEpsE .GT. 0)THEN  
			  IF(pressure .GT. 0) THEN
					stateNew(k,2) = 50
			  ELSE
					stateNew(k,2) = 100
			  ENDIF
		  ELSE
			  stateNew(k,2) = 0
		  ENDIF			  
	   end do
      return 	  
      end
