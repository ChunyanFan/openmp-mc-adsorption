!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE ENSEMBLE_GCMC
!    ARRANGE THE PARTICLES IN THE CUBIC LATTICE CONFIGURATION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
    subroutine ENSEMBLE_GCMC
    USE Constant_M
    USE FLUID_SOLID_M
    USE PHYSICAL_M
    USE MCSETTING_M
    USE FLUCTUATION_M
    USE FLEXICLINFO_M
    USE SUBBOX_M 
    implicit none
    
    integer iP, isN, sNmove,iSub, FLAG
    integer numofinsert, successinsert, numofdelete, successdelete
    integer numberOfMove, numberOfAcceptanceMove
    integer iEqui, iMove, iCycle
    integer i, j
    integer Ncycle1, Ncycles1, numberofequilibrium
    integer NChange
    integer NDD, Dmaxbin, bin
    
        
    real*8 rdn
    real*8 secnds, time
    real*8 EnergyTotal, totAcceptanceRatio
    real*8 Energy, CLEnergy,EnergyC, TEnergy
    real*8 Sum_TEnergy, Sum_Energy, Sum_EnergyC
    real*8 Sum_Npart, Sum_TEN, Sum_EN, Sum_ECN , Sum_NNpart, Sum_Ninpore,Sum_CollectN                 
    real*8 numberOfSampling
    real*8 Av_Npart, Av_Ninpore,Ad_Ninpore, Ns_rho, Nb_rho, Av_CollectN
    real*8 Iso_Energy, Iso_EnergyC, Iso_TEnergy, Iso_Npart
    real*8 Iso_EN, Iso_ECN, Iso_TEN,Iso_NN
    real*8 Iso_HeatFF, Iso_HeatSF,Iso_Heat 
    real*8 Ns_rhoMolPerM2, Nb_rhoMolPerM3, rhoMolPerM3
    real*8 EnergyAverageJoulePerMol
    real*8 Iso_HeatKJoulePerMol, Iso_HeatFFKJoulePerMol, Iso_HeatSFKJoulePerMol                                                                                                  
    real*8 AcceptanceRatio, EnergyAverage
    real*8 Ratio_V 
    real*8 Ex_NB_rho, Fun_NN, Fun_EN, Fun_ENFF, Fun_ENSF
    real*8 CN_BULK, CFUN_EN_BULK, ratio_FUN_N
    real*8 NISO_HEAT, NISO_HEATFF, NISO_HEATSF
    real*8 isoterm(10)
    real*8 Ns_rhoUMolPerM2,  Ex_Nb_rhoMolPerM3
    real*8 NIso_HeatKJoulePerMol, NIso_HeatFFKJoulePerMol,NIso_HeatSFKJoulePerMol
    real*8 IsoTermKJoulePerMol4, IsoTermKJoulePerMol5, IsoTermKJoulePerMol7
    real*8 IsoTermKJoulePerMol1,IsoTermKJoulePerMol2, IsoTermKJoulePerMol6,IsoTermKJoulePerMol8
    real*8 LayerH,ExtEnergyTotal
    real*8 NNISO_HEAT, NNIso_HeatKJoulePerMol
    real*8 refgasV0,NNIsoTerm(3), NNIsoTermKJPerMol(3)
    real*8 TEnergyNew, EnergyNew, CLEnergyNew,EnergyCNew
    real*8 SumCarbonlengthx1, SumCarbonlengthy1, SumBoxlengthZ
    real*8 distanceZ, DdeltaZ
    real*8, Dimension(:),POINTER::sum_DNbin, ave_DNbin, DNbin
    real*8, Dimension(:),POINTER::dz,dzA,DDensity,DDensityKmolperM3  
  
      
              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     FILE FOR ISOSTERIC HEAT
! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

      open( unit=13, file='IsoHeat.txt', status='unknown')
      rewind(13)
      write(13,*)'----------------------------------------------------------------------------------------------------------------------------------------------------------'
      write(13,*)'Column  1: number of pressure point'
      write(13,*)'Column  2: Pressure (Pa)'
      write(13,*)'Column  3: N_Bulk (Particle)'
      write(13,*)'Column  4: Fun_NN_Bulk = f(NG,NG) (Particle)'
      write(13,*)'Column  5: Fun_EN_Bulk = f(UG,NG) (-)'
      write(13,*)'Column  6: Fun_NN = f(N,N)        (Particle)'
      write(13,*)'Column  7: Fun_EN = f(U,N)        (-)'
      write(13,*)'Column  8: Term1  = -f(U,N)/(f(N,N)-f(NG,NG))'
      write(13,*)'Column  9: Term2  = k*T*NG/f(NG,NG)'
      write(13,*)'Column 10: Term3  = f(U,N)/f(N,N)'
      write(13,*)'Column 11: Term4  = f(U,N)/NG'
      write(13,*)'Column 12: Term5  = k*T*(f(N,N)-f(NG,NG))/f(NG,NG)'
      write(13,*)'Column 13: Term6  = (f(U,N)-f(UG,NG))/(f(N,N)-f(NG,NG))'
      write(13,*)'Column 14: Term7  = (f(U,N)-f(UG,NG))/NG'
      write(13,*)'Column 15: Term8  = (f(U,N)-f(UG,NG))/f(N,N)'
      write(13,*)'Column 16: Term9  = f(N,N)/f(NG,NG)'
      write(13,*)'Column 17: Term10 = f(NG,NG)/NG'   
      write(13,'(1A4,16A14)') 'no.','P','CN_Bulk','CFun_NN_Bulk','CFun_EN_Bulk','Fun_NN','Fun_EN','Term1KJPerM','Term2KJPerM','Term3(-)','Term4KJPerM', 'Term5KJPerM','Term6KJPerM','Term7KJPerM','Term8KJPerM','Term9(-)','Term10(-)'   
      write(13,*)'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
      
      write(13,*)'Temperature (-):', Temperature
      
      ratio_V   = AccVolume/BBulkVolume
      do iP = 1, NumOfPre
         CFun_NN_Bulk(iP) = ratio_V*Fun_NN_Bulk(iP)
      enddo
      write(13,*)'AccVolume(-):', AccVolume
      write(13,*)'BBulkVolume(-):', BBulkVolume
      write(13,*)'ratio_V (-):', ratio_V 

       if(SDisplace)then
          NinSBox = 0
          SeqNinSBox = 0
          IndNinSBox = 0
          IF(Npart.gt.0)THEN
             call InitialNpart
          ENDIF
       endif
!       call fluctuation(0)
       IF(RadiusD)call radiusdistribution(0)
       
       
       IF(DistanceD)THEN
       DdeltaZ         = 0.1e0/scalelength
       Dmaxbin         = int(BoxLengthZ/DdeltaZ)+ 1 
       Allocate(sum_DNbin(Dmaxbin), ave_DNbin(Dmaxbin), DNbin(Dmaxbin))
       Allocate(dz(Dmaxbin),dzA(Dmaxbin),DDensity(Dmaxbin),DDensityKmolperM3(Dmaxbin))
       open(unit=31, file='distanceD',status='unknown')
       ENDIF
       
       
       IF(LocalD2D) call DDistribution_2D(0)

    do iPC = 1 , NumOfPre
       write(*,*)'-----------------------------------------------------------------'
       write(*,*)'The pressure being calculated is _ of _ - _:', iPC,NumOfPre,P(iPC)
       write(*,*)'-----------------------------------------------------------------'    
!      ------------
!      START TIMING
!      ------------
       time = secnds(0.0)
!      -------------------------
!      SET THE DISPLACEMENT STEP 
!      -------------------------
       if(SDisplace.AND.(.NOT. TradMC))then
          do i=1, SubBox
             SStepX(i)=IniStepX(i)
             SStepY(i)=IniStepY(i)
             SStepZ(i)=IniStepZ(i)
          enddo
       elseif(TradMC)then
          deltaRx   = deltaRxInitial
          deltaRy   = deltaRyInitial
          deltaRz   = deltaRzInitial
          write(*,*)'deltaRx,deltaRy,deltaRz',deltaRx,deltaRy,deltaRz
       endif
!      --------------------------------
!      INITIAL POSITION ALL N PARTICLES
!      --------------------------------
!       call Initialization   
        write(1,*)P(iPC),rho(iPC)
!      -------------------------------- 
!      EQUILIBRIUM STEPS FOR ADSORPTION 
!      --------------------------------
       numofinsert            = 0
       successinsert          = 0
       numofdelete            = 0
       successdelete          = 0
       numberOfMove           = 0
       numberOfAcceptanceMove = 0
       numberofequilibrium    = 0.0
       sum_Npart              = 0.0
       EnergyTotal            = 0.0
       totAcceptanceRatio     = 0.0 
       SumCarbonlengthx1      = 0.0 
       SumCarbonlengthy1      = 0.0 
       SumBoxlengthZ          = 0.0 
       NChange = 0
       call ConfigEnergy(Energy, CLEnergy,EnergyC)
          TEnergy = Energy + CLEnergy + EnergyC 
          if(ExtraBin)ExtTEnergy = ExtEnergy + ExtCLEnergy
          write(*,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
          write(1,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
    
       LAYERF = .FALSE.
       

          
       IF(SDisplace.AND.(.NOT. TradMC))THEN
             AtpinMove = 0.0e0
             SucinMove = 0.0e0
             SucMovein = 0.0e0
             SucMoveout= 0.0e0
             Ncycle1 = int(Ncycle*1e3/4.0/SubBox)
         do iEqui = 1, Ncycle1

             if (mod(iEqui,1000).eq.0)then
                write(*,*)'==============Mu-GCMC EQULIBRATION==================='
                write(1,*)'==============Mu-GCMC EQULIBRATION==================='
                write(*,646)iEqui, sum_Npart/numberofequilibrium, EnergyTotal/numberofequilibrium
                write(1,646)iEqui, sum_Npart/numberofequilibrium, EnergyTotal/numberofequilibrium
      646       format(1x,'CYCLES:',I8,'  Npart:',e15.6,'  AveEnergy:',E15.6)
                write(*,642)P(iPC)
                write(1,642)P(iPC)
      642       format(1x,'P-pa:',E11.4)
             endif
              
             IF(mod(iEqui,1000).eq.0)THEN 
               if(ChangeXY)then
                  call ShapeChange(TEnergy, Energy, CLEnergy,EnergyC,TEnergyNew, EnergyNew, CLEnergyNew,EnergyCNew)
                  TEnergy = TEnergyNew
                  Energy  = EnergyNew
                  CLEnergy = CLEnergyNew
                  EnergyC  = EnergyCNew
                  if(ExtraBin)ExtTEnergy = ExtEnergy + ExtCLEnergy
                  SumCarbonlengthx1      = SumCarbonlengthx1 + Carbonlengthx1
                  SumCarbonlengthy1      = SumCarbonlengthy1 + Carbonlengthy1
                  SumBoxlengthZ          = SumBoxlengthZ + BoxlengthZ
                  NChange                = NChange + 1
               endif
             ENDIF
                
             do iMove = 1, SubBox !Nmove
                call random_number(rdn)
                iSub = int(rdn*SubBox) + 1
                if(iSub.gt.SubBox) iSub = SubBox      

                 FLAG=0
                 call random_number(rdn)
                 IF(rdn.le.0.5e0)THEN   
                    FLAG=1 !Insert-Displace-Delete-Displace
                 ELSE
                    FLAG=-1 !Delete-Displace-Insert-Displace
                 ENDIF

                 IF(FLAG.EQ.1)THEN
                    call Smcexc1(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Insert
                 ELSEIF(FLAG.EQ.-1)THEN
                    call Smcexc2(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Delete
                 ENDIF
                  numberofequilibrium = numberofequilibrium + 1
                  if(mod(numberofequilibrium,1000).eq.0)then
                     call secAdjustment
                  endif  
                  if((mod(numberofequilibrium,100000).eq.0).and.(numberofequilibrium.LE.1E7))then  
                     write(23,217)P(iPC),numberofequilibrium,(SucMovein(i)-SucMoveout(i),i=1,SubBox) , EnergyTotal/numberofequilibrium    
217                  FORMAT(1x,E11.4,I12,15E14.4,E14.4) 
                     SucMovein = 0
                     SucMoveout = 0
                   endif
                  sum_Npart = sum_Npart + Npart
                  if(Npart .ne. 0)then
                      EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                  endif
                  
                  call Smcexc3(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Displace
                  numberofequilibrium = numberofequilibrium + 1.0e0
                  if(mod(numberofequilibrium,1000).eq.0)then
                     call secAdjustment
                  endif  
                  if((mod(numberofequilibrium,100000).eq.0).and.(numberofequilibrium.LE.1E7))then  
                     write(23,218)P(iPC),numberofequilibrium,(SucMovein(i)-SucMoveout(i),i=1,SubBox), EnergyTotal/numberofequilibrium  
218                  FORMAT(1x,E11.4,I12,15E14.4,E14.4) 
                     SucMovein = 0
                     SucMoveout = 0
                   endif
                  sum_Npart = sum_Npart + Npart
                  if(Npart .ne. 0)then
                      EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                  endif
                  
                 IF(FLAG.EQ.1)THEN
                    call Smcexc2(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Delete
                 ELSEIF(FLAG.EQ.-1)THEN
                    call Smcexc1(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Insert
                 ENDIF
                  numberofequilibrium = numberofequilibrium + 1.0e0
                  if(mod(numberofequilibrium,1000).eq.0)then
                     call secAdjustment
                  endif  
                  if((mod(numberofequilibrium,100000).eq.0).and.(numberofequilibrium.LE.1E7))then  
                     write(23,219)P(iPC),numberofequilibrium,(SucMovein(i)-SucMoveout(i),i=1,SubBox), EnergyTotal/numberofequilibrium      
219                  FORMAT(1x,E11.4,I12,15E14.4,E14.4) 
                     SucMovein = 0
                     SucMoveout = 0
                   endif
                  sum_Npart = sum_Npart + Npart
                  if(Npart .ne. 0)then
                      EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                  endif
                  
                  call Smcexc3(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Displace
                  numberofequilibrium = numberofequilibrium + 1.0e0
                  if(mod(numberofequilibrium,1000).eq.0)then
                     call secAdjustment
                  endif  
                  if((mod(numberofequilibrium,100000).eq.0).and.(numberofequilibrium.LE.1E7))then  
                     write(23,220)P(iPC),numberofequilibrium,(SucMovein(i)-SucMoveout(i),i=1,SubBox) , EnergyTotal/numberofequilibrium 
220                  FORMAT(1x,E11.4,I12,15E14.4,E14.4) 
                     SucMovein = 0
                     SucMoveout = 0
                   endif
                  sum_Npart = sum_Npart + Npart
                  if(Npart .ne. 0)then
                      EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                  endif
             enddo  !iMove
         enddo !iEqui
         write(16,279)P(iPC),(SStepX(i)*ScaleLength,i=1,SubBox),(SStepY(i)*ScaleLength,i=1,SubBox),(SStepZ(i)*ScaleLength,i=1,SubBox)
279      FORMAT(1x,E11.4,45E18.8)
         
          IF(ChangeXY)THEN  !AT THE END OF EQULIBRIM STAGE RESET THE BOX WITH THE AVERAGE BOX LENGTH
             CALL BOXRESET(SumCarbonlengthx1,SumCarbonlengthy1,SumBoxlengthZ,NChange)
          ENDIF

       ELSEIF(TradMC)THEN    !Traditional GCMC
             IF(SDisplace)THEN
                SucMovein = 0.0e0
                SucMoveout= 0.0e0
             ENDIF
          do iEqui = 1, Ncycle
             if (mod(iEqui,1000).eq.0)then
                write(*,*)'================GCMC EQULIBRATION==================='
                write(1,*)'================GCMC EQULIBRATION==================='
                write(*,645)iEqui, sum_Npart/numberofequilibrium, EnergyTotal/numberofequilibrium
                write(1,645)iEqui, sum_Npart/numberofequilibrium, EnergyTotal/numberofequilibrium
      645       format(1x,'CYCLES:',I8,'  Npart:',e15.6,'  AveEnergy:',E15.6)
                write(*,649)P(iPC),totAcceptanceRatio/dble(iEqui),deltaRx
                write(*,650)deltaRy,deltaRz
                write(1,649)P(iPC),totAcceptanceRatio/dble(iEqui),deltaRx
                write(1,650)deltaRy,deltaRz
      649       format(1x,'P-pa:',E11.4, ' Ratio:',e15.6, '  deltaRx:',e11.4)
      650       format(1x,'  deltaRy:',e11.4, '  deltaRz:',e11.4)

             endif
         
             if(.NOT. FlexiCL)RdmCL = 0.0E0
             do iMove = 1, Nmove
                call random_number(rdn)
                if(rdn.lt.RdmCL)then
                   call CLMove(TEnergy,Energy,CLEnergy,EnergyC)
                else if (rdn .ge. RdmCL .AND. rdn .le. (1.0+2.0*RdmCL)/3.0E0)then
                   call mcMove(TEnergy,Energy,CLEnergy,EnergyC, numberOfMove, numberOfAcceptanceMove)
                else
                   call mcexc(TEnergy,Energy, CLEnergy,EnergyC, numofinsert, &
            &              successinsert,numofdelete,successdelete)
                endif
                numberofequilibrium = numberofequilibrium + 1.0e0
                IF(SDisplace)THEN
                   if((mod(numberofequilibrium,100000).eq.0).and.(numberofequilibrium.LE.1E7))then  
                       write(23,316)P(iPC),numberofequilibrium,(SucMovein(i)-SucMoveout(i),i=1,SubBox), EnergyTotal/numberofequilibrium 
316                    FORMAT(1x,E11.4,I12,15E14.4,E14.4) 
                       SucMovein = 0
                       SucMoveout = 0
                   endif
                ENDIF
                
                sum_Npart = sum_Npart + Npart
                if(Npart .ne. 0)then
                EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                endif

             enddo

             call adjustment(iEqui,numberOfMove, numberOfAcceptanceMove,AcceptanceRatio)
             totAcceptanceRatio = totAcceptanceRatio + AcceptanceRatio
         
            
         
         enddo
          
         
          
          IF(SDisplace)THEN    
             write(16,341)P(iPC), deltaRx*ScaleLength,  deltaRy*ScaleLength,  deltaRz*ScaleLength
341          FORMAT(1x,E11.4,3E18.8)
          ENDIF

     ENDIF

     write(*,*) '1284 TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
!    -----------------------------
!    SAMPLING STEPS FOR ADSORPTION
!    -----------------------------
     call ConfigEnergy(Energy, CLEnergy,EnergyC)
     TEnergy = Energy + CLEnergy + EnergyC 
     if(ExtraBin)ExtTEnergy = ExtEnergy + ExtCLEnergy
     write(*,*) 'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
     write(1,*) 'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
    
     Ninpore = 0
     CollectN = 0
     do i = 1, Npart
        call particleposition(i,1)
     enddo 
     
     EnergyTotal            = 0.0E0
     numberOfMove           = 0
     numberOfAcceptanceMove = 0
     numberOfSampling       = 0.0E0    
     numofinsert            = 0
     successinsert          = 0
     numofdelete            = 0
     successdelete          = 0 
     totAcceptanceRatio     = 0.0
     
     Sum_TEnergy            = 0.0E0
     Sum_Energy             = 0.0E0
     Sum_EnergyC            = 0.0E0
     Sum_Npart              = 0.0E0
     Sum_TEN                = 0.0E0
     Sum_EN                 = 0.0E0
     Sum_ECN                = 0.0E0
     Sum_NNpart             = 0.0E0
     Sum_Ninpore            = 0.0E0
     Sum_CollectN           = 0.0E0
     
     
     
     
    IF (DistanceD)THEN
         do i = 1, Dmaxbin
         sum_DNbin(i) = 0.0e0
         ave_DNbin(i) = 0.0e0
         DNbin(i)     = 0.0e0
         enddo
         NDD       = 0
        write(31,*)'Pressure(pa)   distance(A)   LocalD(Kmol/m3)'
        write(31,*)'==========================================='
    ENDIF
    
         
    
     
         
     IF(SDisplace.AND.(.NOT. TradMC))THEN
           totNinSBox = 0.0e0
         
           Ncycles1 = int(Ncycles*1e3/4.0/SubBox) 
           
         do iCycle=1, NCycles1
             AtpinMove = 0.0e0
             SucinMove = 0.0e0
             SucMovein = 0.0e0
             SucMoveout= 0.0e0

             if (mod(iCycle,1000).eq.0)then
                write(*,*)'===============Mu-GCMC SAMPLING==================='
                write(1,*)'===============Mu-GCMC SAMPLING==================='
                write(*,726)iCycle, sum_Npart/numberOfSampling, EnergyTotal/numberOfSampling
                write(1,726)iCycle, sum_Npart/numberOfSampling, EnergyTotal/numberOfSampling
      726       format(1x,'CYCLES:',I8,'  Npart:',e15.6,'  AveEnergy:',E15.6)
                write(*,742)P(iPC)
                write(1,742)P(iPC)
      742       format(1x,'P-pa:',E11.4)
             endif
                       
             do iMove = 1, SubBox !Nmove
                call random_number(rdn)
                iSub = int(rdn*SubBox) + 1
                if(iSub.gt.SubBox) iSub = SubBox      
                
                    FLAG=0
                    call random_number(rdn)
                    IF(rdn.le.0.5e0)THEN   
                       FLAG=1 !Insert-Displace-Delete-Displace
                    ELSE
                       FLAG=-1 !Delete-Displace-Insert-Displace
                    ENDIF
                      
                    IF(FLAG.EQ.1)THEN
                      call Smcexc1(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Insert
                    ELSEIF(FLAG.EQ.-1)THEN
                      call Smcexc2(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Delete
                    ENDIF
                      numberOfSampling = numberOfSampling + 1.0e0
                      Sum_Energy   = Sum_Energy + Energy + CLEnergy
                      Sum_EnergyC  = Sum_EnergyC + EnergyC
                      Sum_TEnergy  = Sum_TEnergy + TEnergy
                      Sum_Npart    = Sum_Npart + dble(Npart)
                      Sum_CollectN = Sum_CollectN + dble(CollectN)
                      if(ExtraBin)then
                         if(Ninpore .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Ninpore)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Ninpore)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Ninpore)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Ninpore)
                         Sum_NNpart   = Sum_NNpart + dble(Ninpore)*dble(Ninpore)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)      
                         ExtEnergyTotal =  ExtEnergyTotal + ExtTEnergy
                      else
                         if(Npart .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Npart)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Npart)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Npart)
                         Sum_NNpart   = Sum_NNpart + dble(Npart)*dble(Npart)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)                     
                      endif
                      
                      do i=1, SubBox
                         totNinSBox(i) = totNinSBox(i)+NinSBox(i)
                      enddo 
                      
                      call Smcexc3(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Displace
                      numberOfSampling = numberOfSampling + 1.0e0
                      Sum_Energy   = Sum_Energy + Energy + CLEnergy
                      Sum_EnergyC  = Sum_EnergyC + EnergyC
                      Sum_TEnergy  = Sum_TEnergy + TEnergy
                      Sum_Npart    = Sum_Npart + dble(Npart)
                      Sum_CollectN = Sum_CollectN + dble(CollectN)
                      if(ExtraBin)then
                         if(Ninpore .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Ninpore)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Ninpore)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Ninpore)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Ninpore)
                         Sum_NNpart   = Sum_NNpart + dble(Ninpore)*dble(Ninpore)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)      
                         ExtEnergyTotal =  ExtEnergyTotal + ExtTEnergy
                      else
                         if(Npart .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Npart)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Npart)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Npart)
                         Sum_NNpart   = Sum_NNpart + dble(Npart)*dble(Npart)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)                     
                      endif
                      do i=1, SubBox
                         totNinSBox(i) = totNinSBox(i)+NinSBox(i)
                      enddo 
                      
                    IF(FLAG.EQ.1)THEN
                      call Smcexc2(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Delete
                    ELSEIF(FLAG.EQ.-1)THEN
                      call Smcexc1(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Insert
                    ENDIF
                      numberOfSampling = numberOfSampling + 1.0e0
                      Sum_Energy   = Sum_Energy + Energy + CLEnergy
                      Sum_EnergyC  = Sum_EnergyC + EnergyC
                      Sum_TEnergy  = Sum_TEnergy + TEnergy
                      Sum_Npart    = Sum_Npart + dble(Npart)
                      Sum_CollectN = Sum_CollectN + dble(CollectN)
                      if(ExtraBin)then
                         if(Ninpore .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Ninpore)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Ninpore)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Ninpore)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Ninpore)
                         Sum_NNpart   = Sum_NNpart + dble(Ninpore)*dble(Ninpore)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)      
                         ExtEnergyTotal =  ExtEnergyTotal + ExtTEnergy
                      else
                         if(Npart .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Npart)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Npart)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Npart)
                         Sum_NNpart   = Sum_NNpart + dble(Npart)*dble(Npart)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)                     
                      endif
                      do i=1, SubBox
                         totNinSBox(i) = totNinSBox(i)+NinSBox(i)
                      enddo 
                      
                      call Smcexc3(iSub,TEnergy,Energy, CLEnergy,EnergyC) !Displace
                      numberOfSampling = numberOfSampling + 1.0e0
                      Sum_Energy   = Sum_Energy + Energy + CLEnergy
                      Sum_EnergyC  = Sum_EnergyC + EnergyC
                      Sum_TEnergy  = Sum_TEnergy + TEnergy
                      Sum_Npart    = Sum_Npart + dble(Npart)
                      Sum_CollectN = Sum_CollectN + dble(CollectN)
                      if(ExtraBin)then
                         if(Ninpore .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Ninpore)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Ninpore)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Ninpore)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Ninpore)
                         Sum_NNpart   = Sum_NNpart + dble(Ninpore)*dble(Ninpore)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)      
                         ExtEnergyTotal =  ExtEnergyTotal + ExtTEnergy
                      else
                         if(Npart .ne. 0)then
                            EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                         endif  
                         Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Npart)
                         Sum_ECN      = Sum_ECN + EnergyC*dble(Npart)
                         Sum_TEN      = Sum_TEN + TEnergy*dble(Npart)
                         Sum_NNpart   = Sum_NNpart + dble(Npart)*dble(Npart)
                         Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)                     
                      endif
                      do i=1, SubBox
                         totNinSBox(i) = totNinSBox(i)+NinSBox(i)
                      enddo 
            enddo  !iMove
            IF(RadiusD)call radiusdistribution(2)
           
            IF(LocalD2D) call DDistribution_2D(2)
          enddo !iCycle    
     ELSEIF(TradMC)THEN
        do iCycle=1, NCycles
     
            if( mod(iCycle,1000).eq.0) then
               write(*,*)'==================GCMC SAMPLING====================='
               write(1,*)'==================GCMC SAMPLING====================='
               write(*,716)iCycle, sum_Npart/numberOfSampling, EnergyTotal/numberOfSampling
               write(1,716)iCycle, sum_Npart/numberOfSampling, EnergyTotal/numberOfSampling
      716      format(1x,'CYCLES:',I8,'  Npart:',e15.6,'  AveEnergy:',e15.6)
               write(*,719)P(iPC),totAcceptanceRatio/dble(icycle),deltaRx
               write(*,720)deltaRy,deltaRz
               write(1,719)P(iPC),totAcceptanceRatio/dble(icycle),deltaRx
               write(1,720)deltaRy,deltaRz
      719      format(1x,'P-pa:',E11.4, 'Ratio:',e15.6, '  deltaR:',e11.4)
      720      format(1x,'  deltaRy:',e11.4, '  deltaRz:',e11.4)
            endif
           
            do iMove=1, Nmove           
                call random_number(rdn)
                if(rdn.lt.RdmCL)then
                   call CLMove(TEnergy,Energy,CLEnergy,EnergyC)
                else if (rdn .ge. RdmCL .AND. rdn .le. (1.0+2.0*RdmCL)/3.0E0)then
                   call  mcMove(TEnergy,Energy,CLEnergy,EnergyC, numberOfMove, numberOfAcceptanceMove) 
                else
                   call  mcexc(TEnergy,Energy,CLEnergy, EnergyC,  numofinsert,&
            &            successinsert,numofdelete,successdelete)
                endif
              
                numberOfSampling = numberOfSampling + 1.0E0
                
                if(Npart .ne. 0)then
                EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
                endif
                             
                Sum_Energy   = Sum_Energy + Energy + CLEnergy
    !            Sum_CLEnergy = Sum_CLEnergy + CLEnergy
                Sum_EnergyC  = Sum_EnergyC + EnergyC
                Sum_TEnergy  = Sum_TEnergy + TEnergy
                Sum_Npart    = Sum_Npart + dble(Npart)
                Sum_EN       = Sum_EN + (Energy+CLEnergy)*dble(Npart)
                Sum_ECN      = Sum_ECN + EnergyC*dble(Npart)
                Sum_TEN      = Sum_TEN + TEnergy*dble(Npart)
                Sum_NNpart   = Sum_NNpart + dble(Npart)*dble(Npart)
                Sum_Ninpore  = Sum_Ninpore + dble(Ninpore)
                Sum_CollectN = Sum_CollectN + dble(CollectN)
             
            enddo
         
!            call adjustment(numberOfMove, numberOfAcceptanceMove,AcceptanceRatio) 
            totAcceptanceRatio = totAcceptanceRatio + AcceptanceRatio 
            successinsert =  0
            successdelete =  0         
            IF(RadiusD)call radiusdistribution(2)
            IF(LOCALF)call fluctuation(2)
            IF(LocalD2D)call DDistribution_2D(2)
            
            IF(DistanceD)THEN
                NDD = NDD + 1
                DNbin = 0.0e0
                do i = 1, Npart
                    DistanceZ = abs(mcz(i))
                    bin = int(DistanceZ/DdeltaZ) + 1
                       if(bin.le.Dmaxbin)then
                           DNbin(bin) = DNbin(bin) + 1.0e0
                       endif
                enddo
                do i = 1, Dmaxbin
                    sum_DNbin(i) = sum_DNbin(i) + DNbin(i)                    
                enddo
            ENDIF
        ENDDO  
    ENDIF
            
                       
!    -----------------------------------------
!    FINISH THE SAMPLING; NOW DO THE AVERAGING 
!    -----------------------------------------
     IF(RadiusD)call radiusdistribution(3)
     
     IF(DistanceD)THEN
        FINALLD = .TRUE.
        do i = 1, Dmaxbin
            dz(i) = DdeltaZ*(i - 0.5)
            dzA(i) = dz(i)*scalelength
            ave_DNbin(i) = sum_DNbin(i)/dble(NDD)
            DDensity(i)  = ave_DNbin(i)/(CarbonLengthx1*CarbonLengthy1*DdeltaZ)
            DDensityKmolperM3(i) = DDensity(i)/6.0230E23/(Scalelength*1.0e-10)**3/1.0e3
            IF(GCMC.AND.FINALLD)write(31,'(3E16.6)')P(iPC), dzA(i),DDensityKmolperM3(i)  
        enddo
     ENDIF
     
     IF(LOCALF)call fluctuation(3)
     IF(LocalD2D)call DDistribution_2D(3)
          
     EnergyAverage = EnergyTotal/numberOfSampling
     Av_Npart      = sum_npart / numberOfSampling
     Av_Ninpore    = Sum_Ninpore/numberOfSampling
     Ad_Ninpore    = Av_Ninpore - rho(iPC) * AccVolume         
     Ns_rho        = Ad_Ninpore / TotSurfaceArea
     Nb_rho        = Av_Ninpore / AccVolume   
     Ex_Nb_rho     = Ad_Ninpore / AccVolume
     IF(DecayESF)THEN
        Av_CollectN   = Sum_CollectN/numberOfSampling
        Ns_rho        = (Av_CollectN-rho(iPC) * AccVolume) / TotSurfaceArea
     ENDIF

     Iso_Energy    = Sum_Energy/numberOfSampling
     Iso_EnergyC   = Sum_EnergyC/numberOfSampling
     Iso_TEnergy   = Sum_TEnergy/numberOfSampling
     Iso_Npart     = Sum_Npart/numberOfSampling
     Iso_EN        = Sum_EN/numberOfSampling
     Iso_ECN       = Sum_ECN/numberOfSampling
     Iso_TEN       = Sum_TEN/numberOfSampling
     Iso_NN        = Sum_NNpart/numberOfSampling
          
      if(ExtraBin)then
        Iso_Npart = Sum_Ninpore/numberOfSampling        
      endif
      
      IF(SDisplace)then
        do i=1, SubBox
           ave_NinSBox(i) = totNinSBox(i)/numberOfSampling 
        enddo  
      ENDIF  
            
!    ------------------------
!    CALCULATE ISOSTERIC HEAT
!    ------------------------
     Fun_NN       = Iso_NN - Iso_Npart*Iso_Npart
     Fun_EN       = Iso_TEN - Iso_TEnergy*Iso_Npart
     Fun_ENFF     = Iso_EN - Iso_Energy*Iso_Npart
     Fun_ENSF     = Iso_ECN - Iso_EnergyC*Iso_Npart
     ratio_Fun_N  = Fun_NN_Bulk(iPC)/Fun_NN
     CN_Bulk      = ratio_V*N_Bulk(iPC)
     CFun_EN_Bulk = ratio_V*Fun_EN_Bulk(iPC)

     IF(Iso_NN .ne. Iso_Npart*Iso_Npart )then
        Iso_HeatFF = (Iso_Energy*Iso_Npart - Iso_EN) / Fun_NN + Temperature 
        Iso_HeatSF = (Iso_EnergyC*Iso_Npart - Iso_ECN) / Fun_NN 
        Iso_Heat   = -Fun_EN / Fun_NN + Temperature
     else 
        Iso_HeatFF = Temperature 
        Iso_HeatSF = 0.0e0
        Iso_Heat   = Temperature 
        write(*,*)'Iso_NN = Iso_Npart*Iso_Npart'   
        write(1,*)'Iso_NN = Iso_Npart*Iso_Npart' 
     endif  

!    --------------------------------------------
!    CALCULATE ISOSTERIC HEAT WITH THE NEW METHOD
!    --------------------------------------------                  
        NIso_Heat    = -Fun_EN/(Fun_NN - CFun_NN_Bulk(iPC)) + Temperature*N_Bulk(iPC)/Fun_NN_Bulk(iPC)
        NIso_HeatFF  = (Iso_Energy*Iso_Npart - Iso_EN) / (Fun_NN- CFun_NN_Bulk(iPC)) + Temperature*N_Bulk(iPC)/Fun_NN_Bulk(iPC) 
        NIso_HeatSF  = (Iso_EnergyC*Iso_Npart - Iso_ECN) / (Fun_NN- CFun_NN_Bulk(iPC))
        
        refgasV0     =  BBulkVolume*Iso_Npart/N_Bulk(iPC)
        NNIsoTerm(1) = Temperature*N_Bulk(iPC)/Fun_NN_Bulk(iPC)
        NNIsoTerm(2) = refgasV0/(refgasV0-AccVolume)*Fun_EN_Bulk(iPC)/Fun_NN_Bulk(iPC)
        NNIsoTerm(3) = -Fun_EN/(Fun_NN - CFun_NN_Bulk(iPC))
        
        NNIso_Heat   =  NNIsoTerm(1) + NNIsoTerm(2)+ NNIsoTerm(3)
         
        IsoTerm(1)   = -Fun_EN/(Fun_NN-CFun_NN_Bulk(iPC))
        IsoTerm(2)   = N_Bulk(iPC)*Temperature/Fun_NN_Bulk(iPC)
        IsoTerm(3)   = Fun_EN/Fun_NN
        IsoTerm(4)   = Fun_EN/CN_Bulk 
        IsoTerm(5)   = Temperature*(Fun_NN-CFun_NN_Bulk(iPC))/(CFun_NN_Bulk(iPC))
        IsoTerm(6)   = (Fun_EN-CFun_EN_Bulk)/(Fun_NN-CFun_NN_Bulk(iPC))
        IsoTerm(7)   = (Fun_EN-CFun_EN_Bulk)/CN_Bulk 
        IsoTerm(8)   = (Fun_EN-CFun_EN_Bulk)/Fun_NN
        IsoTerm(9)   = Fun_NN/CFun_NN_Bulk(iPC)
        IsoTerm(10)  = Fun_NN_Bulk(iPC)/N_Bulk(iPC)
            
!    ----------------------
!    DIMENSIONAL QUANTITIES
!    ----------------------
      Ns_rhoUMolPerM2 = Ns_rho/(SCALELENGTH*1.0E-10)**2/AvogadroNumber*1.0E6 
      Nb_rhoMolPerM3 = Nb_rho/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
      Ex_Nb_rhoMolPerM3  = Ex_Nb_rho/(SCALELENGTH*1.0E-10)**3/AvogadroNumber       

      EnergyAverageJoulePerMol = EnergyAverage*SCALEENERGY*Rg
      
      Iso_HeatKJoulePerMol      =  Iso_Heat*SCALEENERGY*Rg/1.0E3
      Iso_HeatFFKJoulePerMol    =  Iso_HeatFF*SCALEENERGY*Rg/1.0E3
      Iso_HeatSFKJoulePerMol    =  Iso_HeatSF*SCALEENERGY*Rg/1.0E3
      NIso_HeatKJoulePerMol     =  NIso_Heat*SCALEENERGY*Rg/1.0E3
      NIso_HeatFFKJoulePerMol   =  NIso_HeatFF*SCALEENERGY*Rg/1.0E3
      NIso_HeatSFKJoulePerMol   =  NIso_HeatSF*SCALEENERGY*Rg/1.0E3
      
      NNIsoTermKJPerMol(1)      =  NNIsoTerm(1)*SCALEENERGY*Rg/1.0E3
      NNIsoTermKJPerMol(2)      =  NNIsoTerm(2)*SCALEENERGY*Rg/1.0E3
      NNIsoTermKJPerMol(3)      =  NNIsoTerm(3)*SCALEENERGY*Rg/1.0E3
      NNIso_HeatKJoulePerMol    =  NNIso_Heat*SCALEENERGY*Rg/1.0E3
      
      
      IsoTermKJoulePerMol1   = IsoTerm(1)*SCALEENERGY*Rg/1.0E3
      IsoTermKJoulePerMol2   = IsoTerm(2)*SCALEENERGY*Rg/1.0E3
      IsoTermKJoulePerMol6   = IsoTerm(6)*SCALEENERGY*Rg/1.0E3
      IsoTermKJoulePerMol8   = IsoTerm(8)*SCALEENERGY*Rg/1.0E3
      IsoTermKJoulePerMol4   = IsoTerm(4)*SCALEENERGY*Rg/1.0E3 
      IsoTermKJoulePerMol5   = IsoTerm(5)*SCALEENERGY*Rg/1.0E3 
      IsoTermKJoulePerMol7   = IsoTerm(7)*SCALEENERGY*Rg/1.0E3 

!    ------------------
!    DISPLAY THE OUTPUT
!    ------------------
     time = secnds(0.0) - time
        
     write(*,30)  rho(iPC),Ns_rho,Nb_rho,Ex_Nb_rho, AccVolume, Av_Npart,Av_Ninpore,Ad_Ninpore, &
       &          Pressure(iPC), EnergyAverage,Iso_Heat, &
       &          Iso_HeatFF, Iso_HeatSF, &
       &          NIso_Heat, NIso_HeatFF, NIso_HeatSF, &
       &          Ns_rhoUMolPerM2, Nb_rhoMolPerM3, Ex_Nb_rhoMolPerM3,P(iPC), AccVolumeM3, &
       &          EnergyAverageJoulePerMol,Iso_HeatKJoulePerMol, &
       &          Iso_HeatFFKJoulePerMol, Iso_HeatSFKJoulePerMol, &
       &          NIso_HeatKJoulePerMol, NIso_HeatFFKJoulePerMol, &
       &          NIso_HeatSFKJoulePerMol, &
       &          time
!                 Nb_rho, Nb_rhoMolPerM3               
     write(1,30)  rho(iPC),Ns_rho,Nb_rho,Ex_Nb_rho, AccVolume,Av_Npart,Av_Ninpore,Ad_Ninpore, &
       &          Pressure(iPC),  EnergyAverage,Iso_Heat, &
       &          Iso_HeatFF, Iso_HeatSF, &
       &          NIso_Heat, NIso_HeatFF, NIso_HeatSF, &
       &          Ns_rhoUMolPerM2,Nb_rhoMolPerM3,Ex_Nb_rhoMolPerM3,P(iPC), AccVolumeM3, &
       &          EnergyAverageJoulePerMol,Iso_HeatKJoulePerMol, &
       &          Iso_HeatFFKJoulePerMol, Iso_HeatSFKJoulePerMol, &
       &          NIso_HeatKJoulePerMol, NIso_HeatFFKJoulePerMol, &
       &          NIso_HeatSFKJoulePerMol, &
       &          time
!                 Nb_rho, Nb_rhoMolPerM3, 
30           format(1x, 'DIMENSIONLESS OUTPUT',/, &
       &      1x, '---------------------',/, &
       &      1x, 'rho_eos (-) -------------------------------- : ', e20.10,/, &
       &      1x, 'Ns_rho (-) --------------------------------- : ', e20.10,/, &
       &      1x, 'Nb_rho (-) --------------------------------- : ', e20.10,/, &
       &      1x, 'Ex_Nb_rho (-) ------------------------------ : ', e20.10,/, &
!      &      1x, 'LJ-CarbonLength (-) ------------------------ : ', f20.10,/, &
!      &      1x, 'LJ-Porewidth-H (-) ------------------------- : ', f20.10,/, &
       &      1x, 'AccessVolume (-) --------------------------- : ', e20.10,/, &
       &      1x, 'Number of Particles in total---------------- : ', e20.10,/, &
       &      1x, 'Number of Particles in pore----------------- : ', e20.10,/, &
       &      1x, 'Number of Particles be adsorbed------------- : ', e20.10,/, &      
       &      1x, 'LJ-Pressure-set (-) ------------------------ : ', e20.10,/, &
       &      1x, 'Average Energy per Particle (-) ------------ : ', e20.10,/, &
       &      1x, 'Isostetic Heat (-) ------------------------- : ', e20.10,/, &
       &      1x, 'Isostetic Heat FF(-) ----------------------- : ', e20.10,/, &
       &      1x, 'Isostetic Heat SF(-) ----------------------- : ', e20.10,/, &
       &      1x, 'New Isostetic Heat (-) --------------------- : ', e20.10,/, &
       &      1x, 'New Isostetic Heat FF(-) ------------------- : ', e20.10,/, &
       &      1x, 'New Isostetic Heat SF(-) ------------------- : ', e20.10,/, &
       &      1x, 'DIMENSIONAL OUTPUT',/, &
       &      1x, '---------------------',/, &
       &      1x, 'Ns_Density (uMol/m2) ----------------------- : ', e20.10,/, &
       &      1x, 'Nb_Density (KMol/m3) ----------------------- : ', e20.10,/, &
       &      1x, 'Ex_Nb_Density (KMol/m3) -------------------- : ', e20.10,/, &
       &      1x, 'Pressure-set (Pa) -------------------------- : ', e20.10,/, &
       &      1x, 'AccessVolume (m3) -------------------------- : ', e20.10,/, &
       &      1x, 'Average Energy per Particle (J/mol) -------- : ', e20.10,/, &
       &      1x, 'Isostetic Heat (KJ/mol) -------------------- : ', e20.10,/, &
       &      1x, 'Isostetic Heat FF(KJ/mol) ------------------ : ', e20.10,/, &
       &      1x, 'Isostetic Heat SF(KJ/mol) ------------------ : ', e20.10,/, &
       &      1x, 'New Isostetic Heat (KJ/mol) ---------------- : ', e20.10,/, &
       &      1x, 'New Isostetic Heat FF(KJ/mol) -------------- : ', e20.10,/, &
       &      1x, 'New Isostetic Heat SF(KJ/mol) -------------- : ', e20.10,/, &
       &      1x, 'Time (sec) --------------------------------- : ', f20.10/)
       write(1,*)'average of Total Energy(-)',Iso_TEnergy
!              -------------------------------
!              STORE THE POSITION OF PARTICLES
!              -------------------------------
               write(2,*) 'Npart=', Npart
               do i = 1, Npart
                  write(2,435) P(iPC),mcx(i), mcy(i), mcz(i) !, Lequi(i)*scalelength
435               format(5f16.7)
               enddo  
               write(2,*) 'LJSITES' 
               do i = 1, Npart
                  do j = 1, NLJ
                    write(2,435)P(iPC), LJX(i,j),LJY(i,j),LJZ(i,j)
                  enddo
               enddo
               write(2,*) 'COULOMBSITES' 
               do i = 1, Npart
                  do j = 1, NCL
                    write(2,435)P(iPC), CLX(i,j),CLY(i,j),CLZ(i,j)
                  enddo
               enddo

  
!    ----------------------------------
!    STORE THE INFORMATION FOR ISOTHERM
!    ----------------------------------
        write(7,712)iPC,rCutoff,P(iPC),P(iPC)/Psaturate,Ns_rhoUMolPerM2,Nb_rhoMolPerM3,Ex_Nb_rhoMolPerM3,Iso_HeatKJoulePerMol,Iso_HeatFFKJoulePerMol, Iso_HeatSFKJoulePerMol, &
  &                NIso_HeatKJoulePerMol,NIso_HeatFFKJoulePerMol, NIso_HeatSFKJoulePerMol, NNIso_HeatKJoulePerMol,NNIsoTermKJPerMol(1),NNIsoTermKJPerMol(2),NNIsoTermKJPerMol(3), &
  &                Av_Npart,Av_Ninpore,Ad_Ninpore,CarbonLengthx1*scalelength,CarbonLengthy1*scalelength,BoxLengthZ*scalelength  !, Accu_mcBin_ML(MLBin),Accu_mcBin_SL(MLBin)
712     FORMAT(1x,I4,6E18.8,20e18.8)  
        write(13,*)iPC,P(iPC),CN_Bulk,CFun_NN_Bulk,CFun_EN_Bulk,Fun_NN,Fun_EN,IsoTermKJoulePerMol1,IsoTermKJoulePerMol2,IsoTerm(3),IsoTermKJoulePerMol4, IsoTermKJoulePerMol5, IsoTermKJoulePerMol6, IsoTermKJoulePerMol7,IsoTermKJoulePerMol8, &
!        write(13,713)iPC,P(iPC),CN_Bulk,CFun_NN_Bulk,CFun_EN_Bulk,Fun_NN,Fun_EN,IsoTermKJoulePerMol1,IsoTermKJoulePerMol2,IsoTerm(3),IsoTermKJoulePerMol4, IsoTermKJoulePerMol5, IsoTermKJoulePerMol6, IsoTermKJoulePerMol7,IsoTermKJoulePerMol8, &
  &     IsoTerm(9),  IsoTerm(10)  
713     FORMAT(1x,I4,16E16.6)  
 ENDDO
     
    

!        call fluctuation(4)
        IF(DistanceD)THEN
           DEALLOCATE(sum_DNbin, ave_DNbin, DNbin)
           DEALLOCATE(dz,dzA,DDensity,DDensityKmolperM3) 
        ENDIF
        
        DEAllocate(rho, P, Pressure, zactivity,InrefChemP)
        DEAllocate(idealChemicalPotential, ChemicalPotential)
        DEAllocate(ResidualChemicalPotential, ChemicalPotentialJoulepermol) 
        DEAllocate(Fun_NN_Bulk, N_Bulk, Fun_EN_Bulk,CFun_NN_Bulk) 
    
    if(SDisplace)then
       IF(HorizontalBin)THEN
           DEALLOCATE(LayerZL,LayerZH, SubVol)
           DEALLOCATE(LayerYL,LayerYH)
           DEALLOCATE(IniStepX,IniStepY,IniStepZ)
           DEALLOCATE(SStepX,SStepY,SStepZ)
           DEALLOCATE(NinSBox,SeqNinSBox,totNinSBox,ave_NinSBox)
           DEALLOCATE(SubSE, SubSEuMolPerM2)
           DEALLOCATE(AtpinMove, SucinMove,SucMovein,SucMoveout) 
        ELSEIF(VerticalBin)THEN 
           DEALLOCATE(SBoxinBinSec,BinSecYL,BinSecYH)
           DEALLOCATE(BinSecZL,BinSecZH)
           DEALLOCATE(LayerYL,LayerYH, SubVol)
           DEALLOCATE(LayerZL,LayerZH)
           DEALLOCATE(IniStepX,IniStepY,IniStepZ)
           DEALLOCATE(SStepX,SStepY,SStepZ)
           DEALLOCATE(NinSBox,SeqNinSBox,totNinSBox,ave_NinSBox)
           DEALLOCATE(AtpinMove, SucinMove,SucMovein,SucMoveout)
        ENDIF
    endif    


    RETURN
    END SUBROUTINE ENSEMBLE_GCMC


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE ENSEMBLE_NVT
!    ARRANGE THE PARTICLES IN THE CUBIC LATTICE CONFIGURATION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
    subroutine ENSEMBLE_NVT
    USE Constant_M
    USE FLUID_SOLID_M
    USE PHYSICAL_M
    USE MCSETTING_M
    USE FLUCTUATION_M
    USE SUBBOX_M
    implicit none
    
!    integer iden
    integer NumberOfMove, NumberOfAcceptanceMove
    integer iEqui, iMove, iCycle
    integer i, j, Ncount, iSub1,iSub2,fSub,tSub
    integer iSwap, iDisp, dSub, flag
    integer Ncycle1, Ncycles1
    integer numberofequilibrium
    integer Nbranch, ibranch
    integer NChange
    
    real*8 secnds, time,rdn
    real*8 EnergyTotal,totAcceptanceRatio
    real*8 Energy, EnergyC, CLEnergy, TEnergy, EnergyAverage
    real*8 numberOfSampling, Sum_CollectN, Av_CollectN
    real*8 rhoMolPerM3, PressurePa, EnergyAverageJoule
    real*8 AcceptanceRatio, EnergyLRC, VirialLRC, PLRC  
    real*8 VirialTotal, SEinputMolPerM2   
    real*8 FFEnergyTotal, SFEnergyTotal, FFEnergyAverage,SFEnergyAverage
    real*8 FFEnergyAverageJoule, SFEnergyAverageJoule
    real*8 totNcount,ave_Ncount
    real*8 SEnvt,SEnvtuMolPerM2, aveEinBulk  
    real*8 CNBulk, CEBulk, CEBulkJoule  
    real*8 EnergyTotalE, EnergyTotalC
    real*8 ExtEnergyTotal, ExtEnergyAverage
    real*8 Nb_rho, Ex_Nb_rho, Nb_rhoMolPerM3,Ex_Nb_rhoMolPerM3
    real*8 TEnergyNew, EnergyNew, CLEnergyNew,EnergyCNew
    real*8 SumCarbonlengthx1, SumCarbonlengthy1, SumBoxlengthZ, BulkVirialAve 
    real*8 BulkVirialTotal 
    
    open(unit=15, file='NVTisotherm.txt',status='unknown')
    rewind(15)
    write(15,*)'------------------------------------------------------------------------------------------------------'
    write(15,'(4x,A4,9A14)') 'iden','rhogas','P', 'Npart', 'SEinput','Nb_rho','Ex_Nb_rho', 'TotEnergy','FFEnergy', 'SFEnergy'     
    write(15,'(8x,9A14)') '(mol/m3)','(Pa)', '(-)', '(umol/m2)','(mol/m3)','(mol/m3)', '(J)','(J)', '(J)'     
    write(15,*)'------------------------------------------------------------------------------------------------------'
    if(SDisplace)then
       open(unit=22,file='SubBox.txt',status='unknown')
       rewind(22)
       write(22,*)'------------------------------------------------------------------------------------------------------'
       write(22,'(4x,4A8,8A14)') 'iden','Npart','rhogas','P', 'SEoutput', 'TotEnergy','FFEnergy', 'SFEnergy','CorectedNg','CorectedUg','SubSE(i)','NinSub'     
       write(22,'(8x,11A14)') '(-)','(mol/m3)','(Pa)', '(umol/m2)', '(J)','(J)', '(J)','(-)','(-)', '(umol/m2)','(-)'    
       write(22,*)'------------------------------------------------------------------------------------------------------'
    endif
       
!    -------------------------------------------
!    INITIAL POSITION FOR STARTING CONFIGURATION
!    -------------------------------------------
     IF(SDisplace .AND. Npart.gt.0)THEN
         call InitialNpart
     ELSEIF(Npart.EQ.0)THEN
        Npart = int( SEinput(1)*TotSurfaceArea)+1
        call Initialization 
     ENDIF  

!    finalLD = .TRUE.
     Nbranch = 1
     if(NVTDes)Nbranch = 2
      
     IF(RadiusD)call radiusdistribution(0)
     
      
   DO ibranch =1, Nbranch
       do iden = 1 , Numofrho
        write(*,*)'----------------------------NVT ENSEMBLE-------------------------------'
        write(*,*)'The point being calculated is _ of _ ,Npart:', iden,Numofrho,Npart+addN 
        write(*,*)'-----------------------------------------------------------------------'    
!       ------------
!       START TIMING
!       ------------
        time = secnds(0.0)
!       -------------------------
!       SET THE DISPLACEMENT STEP 
!       -------------------------
        if(SDisplace .AND. (.NOT. TradMC))then
          do i = 1, SubBox
             SStepX(i)=IniStepX(i)
             SStepY(i)=IniStepY(i)
             SStepZ(i)=IniStepZ(i)
          enddo
        elseif(TradMC)then
           deltaRx   = deltaRxInitial
           deltaRy   = deltaRyInitial
           deltaRz   = deltaRzInitial
           write(*,*)'deltaRx',deltaRx
        endif
!       --------------------------------
!       INSERT NaddNVT PARTICLES
!       --------------------------------
        if(.NOT. NaddNVT) addN = int(rho(iden)*AccVolume + SEinput(iden)*TotSurfaceArea)+1 - Npart
        if(ibranch.eq.1)then
           call InsertNVT   
           Npart = Npart + addN
        elseif(ibranch.eq.2)then
           call DelNVT
        endif
!       -----------------
!       EQUILIBRIUM STEPS
!       -----------------
        numberOfMove           = 0
        numberOfAcceptanceMove = 0
        numberofequilibrium    = 0
        VirialTotal            = 0.0
        totAcceptanceRatio     = 0.0 
        EnergyTotal            = 0.0
        EnergyTotalE           = 0.0
        EnergyTotalC           = 0.0
        SumCarbonlengthx1      = 0.0
        SumCarbonlengthy1      = 0.0
        SumBoxlengthZ          = 0.0
        NChange                = 0
        call ConfigEnergy(Energy, CLEnergy,EnergyC)
             TEnergy = Energy + CLEnergy + EnergyC 
             if(ExtraBin)ExtTEnergy =ExtEnergy + ExtCLEnergy
             write(*,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
             write(1,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
         
        IF(SDisplace .AND. (.NOT. TradMC))THEN
     
             Nmove=SubBox*(SubBox-1)/2
             Ncycle1=int(Ncycle*1e3/6.0/Nmove)
             AtpinMove = 0
             SucinMove = 0
             SucMovein = 0
             SucMoveout = 0

            do iEqui = 1, Ncycle1

               if (mod(iEqui,1000).eq.0)then
                  write(*,*)'=============Mu-NVT EQULIBRATION==================='
                  write(1,*)'=============Mu-NVT EQULIBRATION==================='
                  write(*,646)iEqui, EnergyTotal/numberofequilibrium/dble(Npart)
                  write(1,646)iEqui, EnergyTotal/numberofequilibrium/dble(Npart)
      646         format(1x,'CYCLES:',I8, '  AveEnergy:',3E15.6)
                  write(*,648)Npart,totAcceptanceRatio/dble(iEqui),deltaRx
                  write(1,648)Npart,totAcceptanceRatio/dble(iEqui),deltaRx
      648         format(1x,'Npart(-):',I5, ' Ratio:',f15.5, '  deltaR:',f10.6)
               endif
               
               if (mod(iEqui,1000).eq.0)then
                  IF(ChangeXY)THEN
                    CALL ShapeChange(TEnergy, Energy, CLEnergy,EnergyC,TEnergyNew, EnergyNew, CLEnergyNew,EnergyCNew)
                    TEnergy = TEnergyNew
                    Energy  = EnergyNew
                    CLEnergy = CLEnergyNew
                    EnergyC  = EnergyCNew
                    if(ExtraBin)ExtTEnergy =ExtEnergy + ExtCLEnergy
                    SumCarbonlengthx1      = SumCarbonlengthx1 + Carbonlengthx1
                    SumCarbonlengthy1      = SumCarbonlengthy1 + Carbonlengthy1
                    SumBoxlengthZ          = SumBoxlengthZ + BoxlengthZ
                    NChange                = NChange + 1
                  ENDIF
               endif

               do iMove = 1, Nmove 
!                 -------------------------
!                  PICK 2 SUBBOXES RANDOMLY
!                 -------------------------
                  call random_number(rdn)
                  iSub1 = int(rdn*SubBox) + 1
                  if(iSub1.gt.SubBox) iSub1 = SubBox   
                  flag=0
                  do while(flag.eq.0)
                     call random_number(rdn)
                     iSub2 = int(rdn*SubBox) + 1
                     if(iSub2.gt.SubBox) iSub2 = SubBox    
                     if(iSub2.NE.iSub1)then
                       flag=1
                       EXIT
                     endif
                  enddo 

                  DO iSwap=1, 2
                     if(iSwap .eq.1)then
                        fSub=iSub1
                        tSub=iSub2
                     elseif(iSwap.eq.2)then
                        fSub=iSub2
                        tSub=iSub1
                     endif
                     call SecMove1(fSub,tSub,TEnergy,Energy,CLEnergy,EnergyC)

                     numberofequilibrium = numberofequilibrium + 1
                     EnergyTotal = EnergyTotal + TEnergy
                     if(ExtraBin)EnergyTotal = EnergyTotal + ExtTEnergy
                         if(mod(numberofequilibrium,1000).eq.0)then
                           call secAdjustment
                         endif  
                         if((mod(numberofequilibrium,100).eq.0).and.(numberofequilibrium.LE.10000))then  
                             write(23,554)Npart,numberofequilibrium,(SucMovein(i),i=1,SubBox), (SucMoveout(i),i=1,SubBox)     
554                          FORMAT(1x,I5,I12,30E14.4) 
                             SucMovein = 0
                             SucMoveout = 0
                         endif

                     do iDisp=1,2
                        if(iDisp.eq.1)then
                           dSub=fSub
                        else
                           dSub=tSub
                        endif
                     
                        if(NinSBox(dSub).eq.0)CYCLE
                        call SecMove2(dSub,TEnergy,Energy,CLEnergy,EnergyC)             
                             numberofequilibrium = numberofequilibrium + 1
                             EnergyTotal = EnergyTotal + TEnergy
                             if(ExtraBin)EnergyTotal = EnergyTotal + ExtTEnergy
                         if(mod(numberofequilibrium,1000).eq.0)then
                           call secAdjustment
                         endif  
                         if((mod(numberofequilibrium,100).eq.0).and.(numberofequilibrium.LE.10000))then  
                             write(23,552)Npart,numberofequilibrium,(SucMovein(i),i=1,SubBox), (SucMoveout(i),i=1,SubBox)     
552                          FORMAT(1x,I5,I12,30E14.4) 
                             SucMovein = 0
                             SucMoveout = 0
                         endif

                     enddo !iDisp
                  ENDDO !iSwap
               enddo !iMove
            enddo !iEqui
            IF(ChangeXY)THEN  !AT THE END OF EQULIBRIM STAGE RESET THE BOX WITH THE AVERAGE BOX LENGTH
               CALL BOXRESET(SumCarbonlengthx1,SumCarbonlengthy1,SumBoxlengthZ,NChange)
            ENDIF
            write(*,*) '1805 TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
            write(16,553)Npart,(SStepX(i)*ScaleLength,i=1,SubBox),(SStepY(i)*ScaleLength,i=1,SubBox),(SStepZ(i)*ScaleLength,i=1,SubBox)
553         FORMAT(1x,I4,45E18.8)
          
        ELSEIF(TradMC)THEN
            IF(SDisplace)THEN
               SucMovein = 0
               SucMoveout = 0
            ENDIF    
            do iEqui = 1, Ncycle

               if (mod(iEqui,1000).eq.0)then
                  write(*,*)'================NVT EQULIBRATION==================='
                  write(1,*)'================NVT EQULIBRATION==================='
                  write(*,642)iEqui, EnergyTotal/numberofequilibrium/dble(Npart)
                  write(1,642)iEqui, EnergyTotal/numberofequilibrium/dble(Npart)
      642         format(1x,'CYCLES:',I8, '  AveEnergy:',E15.6)
                  write(*,644)Npart,totAcceptanceRatio/dble(iEqui),deltaRx
                  write(1,644)Npart,totAcceptanceRatio/dble(iEqui),deltaRx
      644         format(1x,'Npart(-):',I5, ' Ratio:',f15.5, '  deltaR:',f10.6)
               endif
         
         
               do iMove = 1, Nmove     
                  call mcMove(TEnergy,Energy,CLEnergy,EnergyC, numberOfMove, numberOfAcceptanceMove)
                  numberofequilibrium = numberofequilibrium + 1.0e0
                  IF(SDisplace)THEN
                     if((mod(numberofequilibrium,100).eq.0).and.(numberofequilibrium.LE.10000))then  
                        write(23,645)Npart,numberofequilibrium,(SucMovein(i),i=1,SubBox), (SucMoveout(i),i=1,SubBox)     
645                     FORMAT(1x,I5,I12,30E14.4) 
                        SucMovein = 0
                        SucMoveout = 0
                     endif
                  ENDIF

                  EnergyTotal = EnergyTotal + TEnergy
                  VirialTotal = VirialTotal + Virial 
               enddo
               call adjustment(iEqui,numberOfMove, numberOfAcceptanceMove,AcceptanceRatio)
!              totAcceptanceRatio = totAcceptanceRatio + AcceptanceRatio      
             enddo
             IF(SDisplace)THEN
                write(16,557)Npart,deltaRx*ScaleLength,deltaRy*ScaleLength, deltaRz*ScaleLength
557             FORMAT(1x,I4,3E18.8)
             ENDIF

        ENDIF
!       --------------
!       SAMPLING STEPS
!       --------------
        call ConfigEnergy(Energy, CLEnergy,EnergyC)
        TEnergy = Energy + CLEnergy + EnergyC 
        if(ExtraBin)ExtTEnergy = ExtEnergy+ExtCLEnergy
        write(*,*) 'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
        write(1,*) 'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
     
        CollectN = 0
        do i = 1, Npart
           call particleposition(i,1)
        enddo
        Sum_CollectN = 0.0
        
        EnergyTotal            = 0.0E0
        FFEnergyTotal          = 0.0E0
        SFEnergyTotal          = 0.0E0
        numberOfMove           = 0
        numberOfAcceptanceMove = 0
        numberOfSampling       = 0.0E0    
        totAcceptanceRatio     = 0.0
        totNcount              = 0.0
        BulkVirialTotal        = 0.0E0
        if(ExtraBin)then
          ExtEnergyTotal       = 0.0
        endif 
          
        
        IF(RadiusD)call radiusdistribution(1) 
          
        IF(SDisplace .AND. (.NOT. TradMC))THEN

            totNinSBox           = 0.0
            totEinBulk           = 0.0 
         
             Nmove=SubBox*(SubBox-1)/2
             Ncycles1=int(Ncycles*1e3/6.0/Nmove)
            do iCycle = 1, Ncycles1
                AtpinMove = 0
                SucinMove = 0
                SucMovein = 0
                SucMoveout = 0

               if (mod(iCycle,1000).eq.0)then
                  write(*,*)'=============Mu-NVT SAMPLING==================='
                  write(1,*)'=============Mu-NVT SAMPLING==================='
                  if(ExtraBin)then
                     write(*,711)iCycle, (EnergyTotal+ExtEnergyTotal)/numberOfSampling/dble(Npart)
                     write(1,711)iCycle, (EnergyTotal+ExtEnergyTotal)/numberOfSampling/dble(Npart)
      711            format(1x,'CYCLES:',I8, '  AveEnergy:',E15.6)
                  else
                      write(*,711)iCycle, EnergyTotal/numberOfSampling/dble(Npart)
                      write(1,711)iCycle, EnergyTotal/numberOfSampling/dble(Npart)
                  endif
                  write(*,712)Npart,totAcceptanceRatio/dble(iCycle),deltaRx
                  write(1,712)Npart,totAcceptanceRatio/dble(iCycle),deltaRx
      712         format(1x,'Npart(-):',I5, ' Ratio:',f15.5, '  deltaR:',f10.6)
               endif
            
            
               do iMove = 1, Nmove  
                  call random_number(rdn)
                  iSub1 = int(rdn*SubBox) + 1
                  if(iSub1.gt.SubBox) iSub1 = SubBox   
                  flag=0
                  do while(flag.eq.0)
                     call random_number(rdn)
                     iSub2 = int(rdn*SubBox) + 1
                     if(iSub2.gt.SubBox) iSub2 = SubBox    
                     if(iSub2.NE.iSub1)then
                       flag=1
                       EXIT
                     endif
                  enddo 
   
                  DO iSwap=1, 2
                     if(iSwap .eq. 1)then
                        fSub=iSub1
                        tSub=iSub2
                     elseif(iSwap .eq. 2)then
                        fSub=iSub2
                        tSub=iSub1
                     endif
                     call SecMove1(fSub,tSub,TEnergy,Energy,CLEnergy,EnergyC)
                     numberOfSampling  = numberOfSampling  + 1.0e0
                     EnergyTotal = EnergyTotal + TEnergy  !/dble(Npart)
                     if(ExtraBin)ExtEnergyTotal = ExtEnergyTotal + ExtTEnergy
                     FFEnergyTotal = FFEnergyTotal + (Energy+CLEnergy)   ! /dble(Npart)
                     SFEnergyTotal = SFEnergyTotal + EnergyC   !/dble(Npart) 
                     BulkVirialTotal = BulkVirialTotal + BulkVirial
 !                   VirialTotal = VirialTotal + Virial 
                     Sum_CollectN = Sum_CollectN + dble(CollectN)
                     do i=1, SubBox
                        totNinSBox(i) = totNinSBox(i) + dble(NinSBox(i))
                     enddo
                        totEinBulk = totEinBulk + EinBulk
                  
                     do iDisp=1,2
                        if(iDisp.eq.1)then
                           dSub=fSub
                        else
                           dSub=tSub
                        endif
                     
                        if(NinSBox(dSub).eq.0)CYCLE
                        call SecMove2(dSub,TEnergy,Energy,CLEnergy,EnergyC)
                        numberOfSampling  = numberOfSampling  + 1.0e0
                        EnergyTotal = EnergyTotal + TEnergy  
                        if(ExtraBin)ExtEnergyTotal = ExtEnergyTotal + ExtTEnergy
                        FFEnergyTotal = FFEnergyTotal + (Energy+CLEnergy)   
                        SFEnergyTotal = SFEnergyTotal + EnergyC  
                        BulkVirialTotal = BulkVirialTotal + BulkVirial
                        Sum_CollectN = Sum_CollectN + dble(CollectN)
                        do i=1, SubBox
                           totNinSBox(i) = totNinSBox(i) + dble(NinSBox(i))
                        enddo
                        totEinBulk = totEinBulk + EinBulk
                     enddo !iDisp
                  ENDDO !iSwap
               enddo   !iMove
            IF(RadiusD)call radiusdistribution(2)
           
            enddo  !iCycle
        ELSEIF(TradMC)THEN
            IF(SDisplace)THEN
               totNinSBox           = 0.0
               totEinBulk           = 0.0
            ENDIF
          do iCycle=1, NCycles
     
              if( mod(iCycle,1000).eq.0) then
         
                 write(*,*)'==================NVT SAMPLING====================='
                 write(1,*)'==================NVT SAMPLING====================='
                 write(*,715)iCycle, EnergyTotal/numberOfSampling
                 write(1,715)iCycle, EnergyTotal/numberOfSampling
      715        format(1x,'CYCLES:',I8,'  AveEnergy:',e15.6)
                 write(*,718)Npart,totAcceptanceRatio/dble(icycle),deltaRx
                 write(1,718)Npart,totAcceptanceRatio/dble(icycle),deltaRx
      718        format(1x,'Npart(-):',I5, 'Ratio:',f15.5, '  deltaR:',f10.6)
              endif
           
              do iMove=1, Nmove           
                     call  mcMove(TEnergy,Energy,CLEnergy,EnergyC, numberOfMove, numberOfAcceptanceMove)                        
                           numberOfSampling = numberOfSampling + 1.0E0 
                           EnergyTotal = EnergyTotal + TEnergy   !/dble(Npart)
                           FFEnergyTotal = FFEnergyTotal + (Energy+CLEnergy)  !/dble(Npart)
                           SFEnergyTotal = SFEnergyTotal + EnergyC  !/dble(Npart)
                           BulkVirialTotal = BulkVirialTotal + BulkVirial 
                           IF(SDisplace)THEN
                             do i=1, SubBox
                                 totNinSBox(i) = totNinSBox(i) + dble(NinSBox(i))
                             enddo
                             totEinBulk = totEinBulk + EinBulk
                           ELSE
                              call NBulkCheck(Ncount) 
                              totNcount = totNcount + dble(Ncount)
                           ENDIF
              enddo
         
!            call adjustment(iCycle,numberOfMove, numberOfAcceptanceMove,AcceptanceRatio) 
!                 totAcceptanceRatio = totAcceptanceRatio + AcceptanceRatio 
             IF(RadiusD)call radiusdistribution(2)
           
            enddo
        ENDIF
!       -----------------------------------------
!       FINISH THE SAMPLING; NOW DO THE AVERAGING 
!       -----------------------------------------
         
         EnergyAverage = EnergyTotal/numberOfSampling
         FFEnergyAverage = FFEnergyTotal/numberOfSampling
         SFEnergyAverage = SFEnergyTotal/numberOfSampling
         BulkVirialAve = BulkVirialTotal/numberOfSampling
         if(ExtraBin)ExtEnergyAverage = ExtEnergyTotal/numberOfSampling       

         if(SDisplace)then
            do i=1, SubBox
               ave_NinSBox(i) = totNinSBox(i)/numberOfSampling
            enddo
            aveEinBulk = totEinBulk/numberOfSampling
            rho(iden) =  ave_NinSBox(SubBox)/SubVol(SubBox)  
            P(iden) = rho(iden)*Temperature
            SEnvt = (dble(Npart)-rho(iden)*AccVolume)/totSurfaceArea  
            if(ExtraBin)then
               aveEinBulk = ExtEnergyAverage
               SEnvt = (dble(Npart)-rho(iden)*AccVolume-ave_NinSBox(SubBox))/totSurfaceArea
               Nb_rho    = (dble(Npart)-ave_NinSBox(SubBox)) / AccVolume   
               Ex_Nb_rho = (dble(Npart)-ave_NinSBox(SubBox)-rho(iden)*AccVolume) / AccVolume   
            endif
            CNBulk = rho(iden)*AccVolume
            CEBulk = aveEinBulk*AccVolume/SubVol(SubBox)
            CEBulkJoule = CEBulk*SCALEENERGY*kB
            do i=1, SubBox
               SubSE(i) = (ave_NinSBox(i)-rho(iden)*SubVol(i))/totSurfaceArea
               SubSEuMolPerM2(i) = SubSE(i)/(SCALELENGTH*1.0E-10)**2/AvogadroNumber*1.0e6
            enddo
         else 
            ave_Ncount = totNcount/numberOfSampling
            SEnvt = (ave_Ncount-rho(iden)*AccVolume)/totSurfaceArea
            P(iden) = rho(iden)*Temperature + BulkVirialAve/SubVol(SubBox)
         endif
         IF(DecayESF)THEN
            Av_CollectN = Sum_CollectN/numberOfSampling
            SEnvt = (Av_CollectN-rho(iden)*AccVolume)/totSurfaceArea
         ENDIF
         
         IF(RadiusD)call radiusdistribution(3)
         
         call ZLRC(rho(iden), EnergyLRC, VirialLRC)
         PLRC = rho(iden)*VirialLRC
         P(iden) = P(iden) + PLRC
      
!       ----------------------
!       DIMENSIONAL QUANTITIES
!       ----------------------
         rhoMolPerM3              = rho(iden)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
         SEinputMolPerM2          = SEinput(iden)/(SCALELENGTH*1.0E-10)**2/AvogadroNumber
         SEnvtuMolPerM2           = SEnvt/(SCALELENGTH*1.0E-10)**2/AvogadroNumber*1.0e6
         PressurePa               = P(iden)*SCALEENERGY*kB/(SCALELENGTH*1E-10)**3
         EnergyAverageJoule       = EnergyAverage*SCALEENERGY*kB
         FFEnergyAverageJoule     = FFEnergyAverage*SCALEENERGY*kB
         SFEnergyAverageJoule     = SFEnergyAverage*SCALEENERGY*kB
         Nb_rhoMolPerM3    = Nb_rho/(SCALELENGTH*1.0E-10)**3/AvogadroNumber  
         Ex_Nb_rhoMolPerM3 = Ex_Nb_rho/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
      
!       ------------------
!       DISPLAY THE OUTPUT
!       ------------------
        time = secnds(0.0) - time
        
        write(*,50)  AccVolume,Npart,P(iden),&
          &          PLRC, EnergyAverage,&
          &          SEnvtuMolPerM2,PressurePa, AccVolumeM3,EnergyAverageJoule,&
          &          time            
        write(1,50)  AccVolume,Npart,P(iden),&
          &          PLRC, EnergyAverage,&
          &          SEnvtuMolPerM2,PressurePa, AccVolumeM3,EnergyAverageJoule,&
          &          time 
50              format(1x, 'DIMENSIONLESS OUTPUT',/, &
          &      1x, '---------------------',/, &
          &      1x, 'AccessVolume (-) --------------------------- : ', f20.10,/, &
          &      1x, 'Number of Particles in total---------------- : ', I10,/, &
          &      1x, 'LJ-Pressure (-) ---------------------------- : ', f20.10,/, &
          &      1x, 'PLRC (-) ----------------------------------- : ', f20.10,/, &      
          &      1x, 'Average Energy per Particle (-) ------------ : ', f20.10,/, &
          &      1x, 'DIMENSIONAL OUTPUT',/, &
          &      1x, '---------------------',/, &
          &      1x, 'Ns_Density (Mol/m2) ------------------------ : ', e20.10,/, &
          &      1x, 'Pressure (Pa) ------------------------------ : ', e20.10,/, &
          &      1x, 'AccessVolume (m3) -------------------------- : ', e20.10,/, &
          &      1x, 'Average Energy per Particle (J) ------------ : ', e20.10,/, &
          &      1x, 'Time (sec) --------------------------------- : ', f20.10/)
!                 -------------------------------
!                 STORE THE POSITION OF PARTICLES
!                 -------------------------------
                  write(2,*) 'Npart=', Npart
                  do i = 1, Npart
                     write(2,436) mcx(i), mcy(i), mcz(i)
436                  format(3f16.7)
                  enddo  
                  write(2,*) 'LJSITES' 
                  do i = 1, Npart
                     do j = 1, NLJ
                       write(2,436) LJX(i,j),LJY(i,j),LJZ(i,j)
                     enddo
                  enddo
                  write(2,*) 'COULOMBSITES' 
                  do i = 1, Npart
                     do j = 1, NCL
                       write(2,436) CLX(i,j),CLY(i,j),CLZ(i,j)
                     enddo
                  enddo

  
!       ----------------------------------
!       STORE THE INFORMATION FOR ISOTHERM
!       ----------------------------------
           write(15,716)iden,rhoMolPerM3,PressurePa, Npart, SEnvtuMolPerM2,Nb_rhoMolPerM3, Ex_Nb_rhoMolPerM3,EnergyAverageJoule,FFEnergyAverageJoule, SFEnergyAverageJoule     
716        FORMAT(1x,I4,2E18.8,I8,6E18.8)  
           if(SDisplace)then
             write(22,726)iden,Npart,rhoMolPerM3,PressurePa, SEnvtuMolPerM2, EnergyAverageJoule,FFEnergyAverageJoule, SFEnergyAverageJoule,&
&                         CNBulk,CEBulkJoule,(SubSEuMolPerM2(i),i=1,SubBox), (ave_NinSBox(i),i=1,SubBox),CarbonLengthx1*scalelength,CarbonLengthy1*scalelength,BoxLengthZ*scalelength  
726          FORMAT(1x,I4,I8,8E18.8,33E18.8)  
           endif               
         
       enddo
    ENDDO
   
    DEAllocate(rho, SEinput,P) 
    if(SDisplace)then
       IF(HorizontalBin)THEN
           DEALLOCATE(LayerZL,LayerZH, SubVol)
           DEALLOCATE(LayerYL,LayerYH)
           DEALLOCATE(IniStepX,IniStepY,IniStepZ)
           DEALLOCATE(SStepX,SStepY,SStepZ)
           DEALLOCATE(NinSBox,SeqNinSBox,totNinSBox,ave_NinSBox)
           DEALLOCATE(SubSE, SubSEuMolPerM2)
           DEALLOCATE(AtpinMove, SucinMove,SucMovein,SucMoveout)  
       ELSEIF(VerticalBin)THEN 
           DEALLOCATE(SBoxinBinSec,BinSecYL,BinSecYH)
           DEALLOCATE(BinSecZL,BinSecZH)
           DEALLOCATE(LayerYL,LayerYH, SubVol)
           DEALLOCATE(LayerZL,LayerZH)
           DEALLOCATE(IniStepX,IniStepY,IniStepZ)
           DEALLOCATE(SStepX,SStepY,SStepZ)
           DEALLOCATE(NinSBox,SeqNinSBox,totNinSBox,ave_NinSBox)
           DEALLOCATE(AtpinMove, SucinMove,SucMovein,SucMoveout)
        ENDIF
    endif    
    
    RETURN
    END SUBROUTINE ENSEMBLE_NVT


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE INITIALIZATION
!    ARRANGE THE PARTICLES IN THE CUBIC LATTICE CONFIGURATION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
     subroutine Initialization  
     USE FLUID_SOLID_M
     USE MCSETTING_M
     USE PHYSICAL_M
     USE SUBBOX_M
     implicit none
     
     
     integer N, N1, N3
     integer i, j, k, iCount, iSelected 
          

     real*8 rdn
     real*8 BoxLength1, BoxLength2, BoxLength3
     real*8 iMCx0(10000), iMCy0(10000), iMCz0(10000)
     
    
     N = Npart
     if(Npart .eq. 0)return     
     if(SDisplace)then
        NinSBox = 0
        SeqNinSBox= 0
        indNinSBox = 0
     endif

!    --------------------------------------
!    NUMBER OF POSITION IN LINEAR DIMENSION
!    --------------------------------------
     N1         = int(N**(1.E0/3.E0)) + 1
     BoxLength1 = CarbonLengthx1/N1
     BoxLength2 = CarbonLengthy1/N1
     BoxLength3 = BoxLengthZ/N1
!    -------------------------------------------
!    NUMBER OF AVAILABLE POSITIONS IN CUBE IS N3
!    -------------------------------------------
     N3 = N1*N1*N1

!    ------------------------------------------
!    ASSIGN THE COORDINATES OF THE N3 POSITIONS
!    ------------------------------------------
     do k=1, N1
        do j=1, N1
              do i=1, N1
                 iCount = ( (k-1)*N1+ (j-1))*N1 + i
                 iMCx0(iCount) = BoxLength1/2.0 + (i-1)*Boxlength1
                 iMCy0(iCount) = BoxLength2/2.0 + (j-1)*Boxlength2
                 iMCz0(iCount) = BoxLength3/2.0 + (k-1)*Boxlength3
              enddo
        enddo
     enddo
!    ------------------------------------------------------------
!    NOW RANDOMLY DISTRIBUTE N PARTICLES IN N3 POSSIBLE LOCATIONS
!    ------------------------------------------------------------
     do i=1, N
        call random_number(rdn)
        iSelected       = int( rdn*N3) + 1
        if(iSelected.gt.N3) iSelected = N3
          MCx(i)        = iMCx0(iSelected)
          MCy(i)        = iMCy0(iSelected)
          MCz(i)        = iMCz0(iSelected)
          if(SDisplace)then
             do k = 1, SubBox
                if(MCz(i).GE.LayerZL(k) .AND. MCz(i).LT.LayerZH(k))then
                  NinSBox(k) = NinSBox(k)+1
                  SeqNinSBox(k,int(NinSBox(k))) = i
                  indNinSBox(i)=k
                  EXIT
                endif
             enddo
          endif
          do j = 1, NLJ
             LJx(i,j) = LJx0(j) + MCx(i) - MCx0
             LJy(i,j) = LJy0(j) + MCy(i) - MCy0
             LJz(i,j) = LJz0(j) + MCz(i) - MCz0
          enddo
         
          do j =1, NCL
             CLx(i,j) = CLx0(j) + MCx(i) - MCx0
             CLy(i,j) = CLy0(j) + MCy(i) - MCy0
             CLz(i,j) = CLz0(j) + MCz(i) - MCz0
          enddo
!    ----------------------------
!    REMOVE THE SELECTED POSITION
!    ----------------------------
        do j = iSelected, N3-1
            iMCx0(j)    = iMCx0(j+1)
            iMCy0(j)    = iMCy0(j+1)
            iMCz0(j)    = iMCz0(j+1)
        enddo
        N3     = N3 - 1
     enddo
     return
     ENDSUBROUTINE Initialization 


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE CARBONCONFIGURATION
!    ARRANGE THE GRAPHIT LAYER IN THE CUBIC, THIS IS FOR GRAPHIT SURFACE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
     subroutine CarbonConfiguration
     implicit none
      
      
     integer N, Npart, NX, NY, NZ, NC3
     integer NZ1, NC1
     integer i, j, k, iCount, offset, iiCount 
     
     real*8 x(100000), y(100000), z(100000)
     real*8 xc(100000), yc(100000), zc(100000)
     real*8 H, GrapheneLayer, F, H1,BoxLength
     real*8 CarbonLength, BulkVolume,AccVolume, CarbonLengthx, CarbonLengthy
     real*8 PoreWidth, BoxLengthZ, CarbonLength1,CarbonLengthx1, CarbonLengthy1
     real*8 SCALELENGTH, SCALEENERGY  
   
     common/positions/ x, y, z, xc, yc, zc
     common/number/ Npart, NZ, NC3
     common/Box/PoreWidth,BoxLengthZ,CarbonLengthx1,CarbonLengthy1,BoxLength
     common/Volume/BulkVolume, AccVolume
     common/other/ SCALELENGTH, SCALEENERGY, GrapheneLayer
     
     F(i,j,k) = NX*NY*(k-1)+NX*(j-1)+i 
!    -----------------------------------------------
!    NUMBER OF POSITION IN LINEAR DIMENSION (X AXIS)
!    -----------------------------------------------
     CarbonLengthx1 = 4.26E0/SCALELENGTH* int(CarbonLength1/4.26E0)    
     NX             = int(CarbonLengthx1*SCALELENGTH/4.26E0)*2   
!    -----------------------------------------------
!    NUMBER OF POSITION IN LINEAR DIMENSION (Y AXIS)
!    -----------------------------------------------
     CarbonLengthy1 = 2.46E0/SCALELENGTH* int(CarbonLength1/2.46E0)
     NY             = int(CarbonLengthy1*SCALELENGTH/2.46E0)*2
!    ---------------------------------------------------------
!    NUMBER OF AVALIABLE POSITIONS FOR CARBON IN CUBE IS NC3
!    ---------------------------------------------------------
     NC1 = NX*NY*NZ1
       
!    ---------------------------------------------------------
!    ASSIGN THE COORDINATE OF NC3 POSITION IN CUBE
!    ---------------------------------------------------------
    
     do 100 k = 1, NZ1, 2
            do 100 j = 1, NY, 2
                   do 100 i = 1, NX
                      if(mod(i,2).eq.0)then
                         offset = 1
                      else
                         offset = 0
                      endif
     iCount = F(i,j,k)
     xc(iCount) = ((i-1-offset)*2.13E0 + offset*2.84E0 )/SCALELENGTH
     yc(iCount) = ((j-1)*1.23E0)/SCALELENGTH
     zc(iCount) = (k-1)*GrapheneLayer
       
100  continue     
       
     do 200 k = 1, NZ1, 2
            do 200 j = 2, NY, 2
                   do 200 i = 1, NX
                      if(mod(i,2).eq.0)then
                         offset = 0
                      else
                         offset = 1
                      endif
     iCount = F(i,j,k)
     xc(iCount) = ((i-1)*2.13E0 + offset*0.71E0 )/SCALELENGTH
     yc(iCount) = ((j-1)*1.23E0)/SCALELENGTH
     zc(iCount) = (k-1)*GrapheneLayer
       
200  continue    
       
     do 300 k = 2, NZ1, 2
            do 300 j = 1, NY
                   do 300 i = 1, NX
                     iCount = F(i,j,k)
                     iiCount = F(i,j,k-1)
     xc(iCount) =  xc(iiCount) + 0.71E0/SCALELENGTH  
     yc(iCount) =  yc(iiCount) + 1.23E0/SCALELENGTH
     zc(iCount) = (k-1)*GrapheneLayer
       
300  continue    
  
    return    
    end
 
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE MCEXC
!      THIS SUBROUTINE EXCHANGE A PARTICLE WITH THE RESERVOIR
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine  mcexc(TEnergy,Energy,CLEnergy, EnergyC, numofinsert, &
            &          successinsert,numofdelete,successdelete)
       USE FLUID_SOLID_M  
       USE PHYSICAL_M
       USE MCSETTING_M 
       USE FLUCTUATION_M  
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer old,new, i,j
       integer numofinsert,successinsert
       integer numofdelete,successdelete  
       integer bin, bin1
       integer flag, iSub, bflag
      
       real*8 xn,yn,zn
       real*8 TEnergy, Energy, CLEnergy,EnergyC
       real*8 EnergyOld, CLEnergyOld, EnergyCOld, TEnergyOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew, TEnergyNew
       real*8 arg, rdn 
       real*8 arg1
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
                   
       
       call random_number(rdn)
              
       if( rdn .lt. 0.5 )then           ! Choose an old particle to remove !delete if-1
       
            numofdelete = numofdelete + 1
            if(npart .eq. 0) return
            
            call random_number(rdn)
                old = int(Npart*rdn) + 1
                if(old .gt. Npart) old = Npart   
            call energysingleparticle(0,old, EnergyOld, CLEnergyOld, EnergyCOld )
                TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld 
                arg = dble(Npart)*exp(TEnergyOld/Temperature)/(zactivity(iPC)*BulkVolume)
  !          arg1 = dble(Npart)*debro3*exp((TEnergyOld-Chemicalpotential(iPC))/Temperature)/BulkVolume
            
            
        
            call random_number(rdn)
            if(rdn .lt. arg)then     ! DELETE   !delete if-2
               IF(SDisplace)THEN 
                  flag=0
                  DO WHILE(flag.eq.0)
                     do j=1,SubBox
                        if( HorizontalBin .AND.(MCz(old).GE.LayerZL(j)) .AND.  (MCz(old).LT.LayerZH(j)) )then
                           iSub = j
                           SucMoveout(iSub)=SucMoveout(iSub)+1
                           flag=1
                           EXIT
                        elseif( VerticalBin .AND.(MCy(old).GE.LayerYL(j)) .AND.  (MCy(old).LT.LayerYH(j)) )then
                           iSub = j
                           SucMoveout(iSub)=SucMoveout(iSub)+1
                           flag=1
                           EXIT
                        endif
                      enddo
                  ENDDO                  
                ENDIF
               call particleposition(old,-1)         
               MCx(old) = MCx(npart)
               MCy(old) = MCy(npart)
               MCz(old) = MCz(npart)
               SFenergy(old) = SFenergy(npart)
               do j = 1, NLJ
                  LJx(old,j) = LJx(npart,j) 
                  LJy(old,j) = LJy(npart,j) 
                  LJz(old,j) = LJz(npart,j) 
               enddo
               do j = 1, NCL
                 CLx(old,j)  = CLx(npart,j)  
                 CLy(old,j)  = CLy(npart,j) 
                 CLz(old,j)  = CLz(npart,j) 
                 ECLx(old,j) = ECLx(npart,j)  
                 ECLy(old,j) = ECLy(npart,j) 
                 ECLz(old,j) = ECLz(npart,j) 
               enddo
                     
               successdelete = successdelete + 1
               Energy   = Energy - EnergyOld
               CLEnergy = CLEnergy - CLEnergyOld
               EnergyC  = EnergyC - EnergyCOld
               TEnergy  = TEnergy - TEnergyOld

               IF(ADSORPTION.AND. LOCALF)THEN    !delete if-3
                  bin = I_POS(old)

                  do j = 1, NinBin(bin)
                     if(SeqNinBin(bin,j).NE.old)then
                        UFFinBin(bin)=UFFinBin(bin)-EnergyMatrix(old,SeqNinBin(bin,j))
                     endif
                  enddo
                  USFinBin(bin) = USFinBin(bin) - EnergyMatrix(old,old)
                  do j = 1, NinBin(bin)
                     if(SeqNinBin(bin,j) .EQ. old)then
                     SeqNinBin(bin,j) = SeqNinBin(bin,NinBin(bin))
                     EXIT
                     endif
                  enddo
                  NinBin(bin) = NinBin(bin) -1
                  if(NinBin(bin).lt.0)then
                  write(*,*)'line 1633 negative NinBin,bin',NinBin(bin),bin
                  read(*,*)
                  endif
                  do j = 1,Npart-1
                     if(j.NE.old)then
                       EnergyMatrix(old,j) =  EnergyMatrix(Npart,j)
                       EnergyMatrix(j,old) =  EnergyMatrix(Npart,j)
                     endif
                  enddo   
                  EnergyMatrix(old,old) = EnergyMatrix(Npart,Npart)     
                  bin1 = I_POS(npart)
                  do j = 1, NinBin(bin1)
                     if(SeqNinBin(bin1,j) .EQ. npart)then
                     SeqNinBin(bin1,j) = old
                     EXIT
                     endif
                  enddo 
                   I_POS(old) = bin1 
                   Npart      = Npart - 1

               ELSE               !delete if-3
                  Npart = Npart - 1
               ENDIF          !delete if-3
            endif    !delete if-2
               
            return               
            
!               do j = 1, Npart-1
!                  if(j.NE.old)then
!                       EnergyMatrix(old,j) =  EnergyMatrix(Npart,j)
!                       EnergyMatrix(j,old) =  EnergyMatrix(Npart,j)
!                     endif
!                enddo   
!                EnergyMatrix(old,old) = EnergyMatrix(Npart,Npart)     
!                Npart    = Npart - 1
!                  WRITE(*,*)'1595 after DELETE Energy,EnergyOld,Npart',Energy,EnergyOld,Npart
!                   write(*,*)(Energymatrix(1:npart, 1:npart))
!                  read(*,*)
    
        else                           ! Insert  !delete if-1

            numofinsert = numofinsert + 1 
            call random_number(rdn)
            xn = rdn * CarbonLengthx1
            call random_number(rdn)
            yn = rdn * CarbonLengthy1
            call random_number(rdn)
            zn = rdn * BoxLengthZ   
                            
            if( InclinePore .and.  &
               &   zn .gt. IncCornerZ  .and. &
               &   yn .lt. ((zn-IncCornerZ)/tan(IncAngle))   )then
                   return 
            endif     
               
            if( SurfPore .and.  &
               &   zn .gt. CriLengthZ  .and. &
               &   yn .lt. CarbonLengthy0 .and. &
               &   yn .lt. ((zn-CriLengthZ)/tan(surftheta))   )then
                return 
            endif 
            
             if(BojanSlit)then      
                 bflag = 0
                 do j = 1, NS
                     if( yn.gt.StripLy(j) .AND. yn.le.StripHy(j)  )then
                         if(   StripAngleY(j) .GT. 0.0 .and. &
             &                 zn .lt. InLzLyStripZ(j) .and. &
             &                (InLzLyStripZ(j)-zn)   .GE. (yn-StripLy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT
                        elseif( StripAngleY(j) .GT. 0.0 .and. &
             &                 zn .gt. InHzLyStripZ(j) .and. &
             &                 (zn-InHzLyStripZ(j))   .GE. (yn-StripLy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT  
                        elseif(  StripAngleY(j) .LT. 0.0 .and. &
             &                   zn .lt. InLzLyStripZ(j) .and. &
             &                 (InLzLyStripZ(j)-zn)   .GE. (yn-StripHy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT
                        elseif( StripAngleY(j) .LT. 0.0 .and. &
             &                 zn .gt. InHzLyStripZ(j) .and. &
             &                 (zn-InHzLyStripZ(j))   .GE. (yn-StripHy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT 
                        elseif(   StripAngleY(j) .EQ. 0.0 .and. &
             &                  ( zn .lt. InLzLyStripZ(j) .or. zn .gt. InHzLyStripZ(j) ) ) then 
                               bflag = 1
                               EXIT
                        endif
                   endif
                enddo
                if(bflag .eq. 1)return
             endif                       
            
            call positioncheck(xn,yn,zn)

!---------------------------------------------------------------
!  
!     FOR MULTI-SITE MOLECULAR SIMULATION
!
!---------------------------------------------------------------
!-------------------------------- 
!     FOR ORIENTATIONAL MOVES    
!--------------------------------
!--------------------
!     Rotate
!--------------------
          CALL ROTATE
          
          do j = 1, NLJ
             LJxnew(j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j)
             LJynew(j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j)
             LJznew(j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j)
          enddo

          do j = 1, NCL
            CLxnew(j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) 
            CLynew(j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) 
            CLznew(j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) 
          enddo  
                 
          new  = npart + 1
          MCx(new) = xn
          MCy(new) = yn
          MCz(new) = zn
          do i = 1, NLJ
             LJx(new,i) = LJxnew(i) + MCx(new) - MCx0
             LJy(new,i) = LJynew(i) + MCy(new) - MCy0
             LJz(new,i) = LJznew(i) + MCz(new) - MCz0
          enddo
         
          do i =1, NCL
             CLx(new,i) = CLxnew(i) + MCx(new) - MCx0
             CLy(new,i) = CLynew(i) + MCy(new) - MCy0
             CLz(new,i) = CLznew(i) + MCz(new) - MCz0
             ECLx(new,i) = CLx(new,i)
             ECLy(new,i) = CLy(new,i)
             ECLz(new,i) = CLz(new,i)
          enddo
          
          call energysingleparticle(0,new,EnergyNew, CLEnergyNew, EnergyCNew)
          TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
          
                          
          arg = zactivity(iPC)*BulkVolume*exp(-TEnergyNew/Temperature)/dble(new)     
!         arg1 = BulkVolume*exp((ChemicalPotential(iPC)-TEnergyNew)/Temperature)/debro3/dble(new)
            call random_number(rdn)
            if(rdn .lt. arg)then
              IF(SDisplace)THEN 
                 flag=0
                 DO WHILE(flag.eq.0)
                    do j=1,SubBox
                       if( HorizontalBin .AND.(MCz(new).GE.LayerZL(j)) .AND.  (MCz(new).LT.LayerZH(j)) )then
                          iSub = j
                          SucMovein(iSub)=SucMovein(iSub)+1
                          flag=1
                          EXIT
                       elseif( VerticalBin .AND.(MCy(new).GE.LayerYL(j)) .AND.  (MCy(new).LT.LayerYH(j)) )then
                          iSub = j
                          SucMovein(iSub)=SucMovein(iSub)+1
                          flag=1
                          EXIT
                       endif
                    enddo
                 ENDDO
               ENDIF

               SFenergy(new)   = EnergyCNew
               successinsert = successinsert + 1
               Npart         = Npart + 1
               Energy        = Energy + EnergyNew
               CLEnergy      = CLEnergy + CLEnergyNew
               EnergyC       = EnergyC + EnergyCNew
               TEnergy       = TEnergy + TEnergyNew
               call particleposition(new,1)         

               IF(ADSORPTION.AND. LOCALF)THEN 
                  call BinJudge(MCx(new),MCy(new),MCz(new),bin) 
                  EnergyMatrix(new,new) = EnergyCNew
                 do j = 1, NinBin(bin)
                      UFFinBin(bin) = UFFinBin(bin) + EnergyMatrix(new,SeqNinBin(bin,j))
                 enddo
                 USFinBin(bin) = USFinBin(bin) + EnergyMatrix(new,new)
                 NinBin(bin) = NinBin(bin) + 1
                 SeqNinBin(bin,NinBin(bin)) = new
                 I_POS(new) = bin
              ENDIF  
              return
            endif
                     
        endif        
         
       return
       endsubroutine  mcexc
          
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ENERGYSINGLEPARTICLE
!      THIS SUBROUTINE CALCULATES THE ENERGY OF ONE PARTICLE WITH THE OTHER 
!      PARTICLES AND THE SOLID
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine energysingleparticle(FLAG,i, Energy, CLEnergy, EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE BUCKINGHAM_M
       USE FLUCTUATION_M 
       USE SUBBOX_M 
       implicit none
       
       integer i, j,k,m, FLAG

       real*8  Energy, EnergyC, CLEnergy
       real*8  d2, dc2,cld, distanceZ
       real*8  U, V,UC, CLU, CLV
       real*8  xlj, ylj, zlj, xclj, yclj, zclj,xmclj, ymclj, zmclj
       real*8  xcmcij,ycmcij,zcmcij,xcl, ycl, zcl
       real*8  x1, y1, z1
       real*8  sigma,welldepth,sigmaSF,welldepthSF
       real*8  xmcn,ymcn,zmcn,mcd2,smcy
       real*8  xcmcn,ycmcn,zcmcn,cmcd2
       real*8  BU, gSM
       real*8  deltay,deltay1,cy1
       real*8  LJy1(20), sLJy(20)
       real*8  pEnergy, EnergyCi, EnergyCij   
       real*8  MCVectorX, MCVectorY, MCVectorZ, MCVectorX0, MCVectorY0, MCVectorZ0
       real*8  LJVectorX, LJVectorY, LJVectorZ, LJVectorX0, LJVectorY0, LJVectorZ0
       real*8  CLVectorX, CLVectorY, CLVectorZ, CLVectorX0, CLVectorY0, CLVectorZ0
       real*8  VectorProduct
        
!      --------------------------------
!      INITIALISE THE ENERGY 
!      --------------------------------
       Energy   = 0.0E0
       EnergyC  = 0.0E0
       ExtEnergy1   = 0.0E0
       ExtCLEnergy1 = 0.0E0
       CLEnergy = 0.0E0
       Virial   = 0.0E0
       VirialCL = 0.0E0
       VirialF  = 0.0E0
       VirialC  = 0.0E0
!      ---------------------------------------------------
!      INTERACTION ENERGY BETWEEN PARTICLE AND CARBON ATOM
!      ---------------------------------------------------
       IF((.NOT. ADSORPTION).OR.(.NOT. LOCALF).OR.(i .GT. Npart).OR. (flag.EQ.1))THEN
          if(steele)then
                if(ExtraBin)then
                   if(indNinSBox(i).lt.SubBox)then
                      call steelepotential(i,UC)
                      EnergyC = EnergyC + UC 
                   endif
                else
                   call steelepotential(i,UC)
                   EnergyC = EnergyC + UC 
!                  VirialC = VirialC + VC
                endif                   
          endif
          if(CrowellChang)then
             call CrowellChangPotential(i,UC)
             EnergyC = EnergyC + UC 
          endif


          if(Bojan)then
!$OMP  parallel do default(none) &
!$OMP   private (j,k,deltay,deltay1,cy1,BU,smcy,LJy1,sLJy) &
!$OMP   shared (i,NS,extrabin,indninsbox,carbonlengthy1,BoPBCy,LJy,mcy) &
!$OMP   shared (NLJ,StripC,StripW,SubBox,VC) &
!$OMP   reduction (+:EnergyC,VirialC)
              do j = 1, NS
                  IF(ExtraBin)THEN
                     if(indNinSBox(i).lt.SubBox)then
                        call BojanPotential(j, i,BU) 
                        EnergyC = EnergyC + BU
                      endif
                  ELSE
                      call BojanPotential(j, i,BU) 
                      EnergyC = EnergyC + BU
                  ENDIF
                   if(BoPBCy)then
                      deltay  = abs(mcy(i) - (StripC(j)+ StripW(j)/2.0e0)) 
                      deltay1 = abs(mcy(i) - (StripC(j)- StripW(j)/2.0e0))  
                      if(deltay .lt. deltay1)  deltay = deltay1
                      if(deltay .gt. (carbonlengthy1/2.0e0))then
                          if((mcy(i).ge.0.0e0) .and. (mcy(i).le.CarbonLengthy1/2.0e0))then
                             cy1 = mcy(i) + carbonlengthy1
                             do k = 1, NLJ
                                LJy1(k) = LJy(i,k) + carbonlengthy1
                             enddo
                          else
                             cy1 = mcy(i) - carbonlengthy1
                             do k = 1, NLJ
                                LJy1(k) = LJy(i,k) - carbonlengthy1
                             enddo
                          endif
                          smcy = mcy(i)
                          mcy(i) = cy1
                          do k = 1, NLJ
                             sLJy(k) = LJy(i,k)
                             LJy(i,k) = LJy1(k)
                          enddo
                          call BojanPotential(j, i,BU) 
                          EnergyC = EnergyC + BU 
                          VirialC = VirialC + VC
                          mcy(i) = smcy 
                          do k = 1, NLJ
                             LJy(i,k) = sLJy(k)
                          enddo
                     endif
                  endif 
             enddo
             IF(inclinepore)then
               call BojanPotential(NS+1, i,BU) 
                   EnergyC = EnergyC + BU 
             ENDIF
             IF(Close1End.OR.Close2Ends)then
               call BojanPotential(NS+2, i,BU) 
                   EnergyC = EnergyC + BU 
             ENDIF  
         endif
       
!$OMP  parallel do default(none)                                         &
!$OMP   private (k,j,xcmcij,ycmcij,zcmcij,xcmcn,ycmcn,zcmcn,cmcd2,xclj,yclj,zclj,sigmaSF,welldepthSF,UC,dc2) &
!$OMP   shared (nlj,nc3,pbcx,pbcy,pbcz,xc,yc,zc,mcx,mcy,mcz,r2Cutoff,VC,ljx,ljy,ljz,i) &
!$OMP   shared (CarbonLengthx1,CarbonLengthy1,BoxLengthZ,sigmaSS,sigmaFF,welldepthSS,welldepthFF) &
!$OMP   reduction (+:EnergyC,VirialC)
          do k = 1, NLJ
             do j= 1, NC3
!               --------------------
!               DISTANCES IN X, Y, Z
!               --------------------
                xcmcij    = abs( xc(j)-mcx(i) )
                ycmcij    = abs( yc(j)-mcy(i) )
                zcmcij    = abs( zc(j)-mcz(i) )
                xcmcn     = xcmcij
                ycmcn     = ycmcij
                zcmcn     = zcmcij
                if(pbcx)then
                  if(xcmcij > CarbonLengthx1/2.0D0) xcmcn = xcmcn - CarbonLengthx1
                endif
                if(pbcy)then
                  if(ycmcij > CarbonLengthy1/2.0D0) ycmcn = ycmcn - CarbonLengthy1
                endif
                if(pbcz)then
                    if(zcmcij > BoxLengthZ/2.0D0) zcmcn = zcmcn - BoxLengthZ
                endif
                cmcd2 = xcmcn*xcmcn + ycmcn*ycmcn + zcmcn*zcmcn
                if(cmcd2.gt.r2Cutoff)then
                   CYCLE
                else
                   sigmaSF  = (sigmaSS + sigmaFF(k))/2.0e0
                   welldepthSF = sqrt(welldepthSS*welldepthFF(k))
                   xclj      = abs( xc(j)-ljx(i,k))
                   yclj      = abs( yc(j)-ljy(i,k))
                   zclj      = abs( zc(j)-ljz(i,k))
!               ----------------------------------------------------
!               ENSURE THAT WE DEAL WITH THE NEAREST PERIODIC IMAGES
!               ----------------------------------------------------
                  if(pbcx)then
                     if(xcmcij > CarbonLengthx1/2.0D0) xclj = xclj - CarbonLengthx1
                  endif
                  if(pbcy)then
                     if(ycmcij > CarbonLengthy1/2.0D0) yclj = yclj - CarbonLengthy1
                  endif
                  if(pbcz)then
                     if(zcmcij > BoxLengthZ/2.0D0) zclj = zclj - BoxLengthZ
                  endif
!               -----------------------
!               INTER-PARTICLE DISTANCE
!               -----------------------
                dc2 = xclj*xclj + yclj*yclj + zclj*zclj
!               ----------------
!               POTENTIAL ENERGY
!               ----------------
                call potentialEnergyC( dc2,sigmaSF,welldepthSF,UC)
                EnergyC = EnergyC + UC 
                VirialC = VirialC + VC
              endif
            enddo
          enddo 
        ELSE
          EnergyC = EnergyMatrix(i,i)
        ENDIF   
          EnergyCi = EnergyC    
!      ------------------------------------
!      INTERACTION ENERGY BETWEEN PARTICLES
!      ------------------------------------
      if(i.gt.0) then   !! 1-if
        if((.NOT. ADSORPTION).OR.(.NOT. LOCALF).OR.(i .GT. Npart) .OR. (FLAG .eq. 1))then  !! 2-if
!$OMP  parallel do default(none)                                         &
!$OMP    private (j,pEnergy,d2,mcd2,VectorProduct,k,m,xlj,ylj,zlj,xmclj,ymclj,zmclj) &
!$OMP    private (MCVectorX0,MCVectorY0,MCVectorZ0,MCVectorX,MCVectorY,MCVectorZ) &
!$OMP    private (LJVectorX0,LJVectorY0,LJVectorZ0,LJVectorX,LJVectorY,LJVectorZ) &
!$OMP    private (CLVectorX0,CLVectorY0,CLVectorZ0,CLVectorX,CLVectorY,CLVectorZ) &
!$OMP    private (cld,xcl,ycl,zcl,xmcn,ymcn,zmcn,U,V,CLU,CLV,EnergyCij,welldepth) &
!$OMP    private (sigma,gSM) &
!$OMP    shared (Npart,indNinSBox,SubBox,i,ExtraBin,CarbonLengthx1,CarbonLengthy1,BoxLengthZ) &
!$OMP    shared (mcx,mcy,mcz,LOCALF,ADSORPTION,flag,ljx,ljy,ljz,clx,cly,clz) &
!$OMP    shared (BUCKINGHAM,ksi,ksif,EnergyMatrix,NCL,EnergyCi,SFEnergy,smediation1,smediation2) &
!$OMP    shared (R2Cutoff,welldepthFF,NLJ,pbcx,pbcy,pbcz,sigmaFF,SMLimitZ,Temperature,charges) &
!$OMP    reduction (+:CLEnergy,VirialCL,ExtCLEnergy1,ExtEnergy1,Energy,Virial,VirialF)
           do j = 1, Npart !! 1-do
              if(ExtraBin)then
                if( (indNinSBox(i).eq.SubBox) .and. (indNinSBox(j).lt.SubBox) )CYCLE
                if( (indNinSBox(i).lt.SubBox) .and. (indNinSBox(j).eq.SubBox) )CYCLE
              endif
              pEnergy = 0.0e0
              if(i .ne. j)then !! 5-if
              !!---------------------------------------------
              !! Lennard-Jones interaction
              !!---------------------------------------------
                xmclj = abs(mcx(j)-mcx(i))
                ymclj = abs(mcy(j)-mcy(i))
                zmclj = abs(mcz(j)-mcz(i))
                xmcn = xmclj
                ymcn = ymclj
                zmcn = zmclj
                MCVectorX0 = mcx(j)-mcx(i)
                MCVectorY0 = mcy(j)-mcy(i)
                MCVectorZ0 = mcz(j)-mcz(i)
                MCVectorX = MCVectorX0
                MCVectorY = MCVectorY0
                MCVectorZ = MCVectorZ0
                if(pbcx)then                
                     if(xmclj > CarbonLengthx1/2.0D0)then
                        xmcn = xmclj - CarbonLengthx1
                        IF(MCVectorX0 .GT. 0.0)MCVectorX = MCVectorX0 - CarbonLengthx1
                        IF(MCVectorX0 .LT. 0.0)MCVectorX = MCVectorX0 + CarbonLengthx1
                     endif
                endif
                if(pbcy)then
                     if(ymclj > CarbonLengthy1/2.0D0)then
                        ymcn = ymclj - CarbonLengthy1
                        IF(MCVectorY0 .GT. 0.0)MCVectorY = MCVectorY0 - CarbonLengthy1
                        IF(MCVectorY0 .LT. 0.0)MCVectorY = MCVectorY0 + CarbonLengthy1
                     endif
                endif
                if(pbcz)then
                     if(zmclj > BoxLengthZ/2.0D0)then
                        zmcn = zmclj - BoxLengthZ
                        IF(MCVectorZ0 .GT. 0.0)MCVectorZ = MCVectorZ0 - BoxLengthZ
                        IF(MCVectorZ0 .LT. 0.0)MCVectorZ = MCVectorZ0 + BoxLengthZ
                     endif
                endif
                mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
                if(mcd2.gt.r2Cutoff)then !! 4-if
                       IF(ADSORPTION .AND. LOCALF)THEN
                           if(flag.eq.0)then
                             EnergyMatrix(Npart+1,j) = 0.0e0
                             EnergyMatrix(j,Npart+1) = 0.0e0
                           else
                             EnergyMatrix(i,j) = 0.0e0
                             EnergyMatrix(j,i) = 0.0e0
                           endif
                        ENDIF
!                    EnergyMatrix(i,j) = 0.0e0
!                    EnergyMatrix(j,i) = 0.0e0
                   CYCLE
                else
                  do k=1,NLJ
                     do m = 1, NLJ
                        sigma = (sigmaFF(k) + sigmaFF(m))/2.0e0
                        welldepth = sqrt(welldepthFF(k)*welldepthFF(m))
!                       --------------------
!                       DISTANCES IN X, Y, Z
!                       --------------------
                        xlj    = abs(ljx(j,m)-ljx(i,k))
                        ylj    = abs(ljy(j,m)-ljy(i,k))  
                        zlj    = abs(ljz(j,m)-ljz(i,k))
                        LJVectorX0 = ljx(j,m)-ljx(i,k)
                        LJVectorY0 = ljy(j,m)-ljy(i,k)
                        LJVectorZ0 = ljz(j,m)-ljz(i,k)
                        LJVectorX = LJVectorX0
                        LJVectorY = LJVectorY0
                        LJVectorZ = LJVectorZ0
   
!                  ----------------------------------------------------
!                  ENSURE THAT WE DEAL WITH THE NEAREST PERIODIC IMAGES
!                  ----------------------------------------------------
                      if(pbcx)then                
                           if(xmclj > CarbonLengthx1/2.0D0)then
                              xlj = xlj - CarbonLengthx1
                              IF(LJVectorX0.GT.0.0)LJVectorX = LJVectorX0 - CarbonLengthx1
                              IF(LJVectorX0.LT.0.0)LJVectorX = LJVectorX0 + CarbonLengthx1
                           endif
                      endif
                      if(pbcy)then
                          if(ymclj > CarbonLengthy1/2.0D0)then
                             ylj = ylj - CarbonLengthy1
                             IF(LJVectorY0.GT.0.0)LJVectorY = LJVectorY0 - CarbonLengthy1
                             IF(LJVectorY0.LT.0.0)LJVectorY = LJVectorY0 + CarbonLengthy1
                          endif
                      endif
                      if(pbcz)then
                           if(zmclj > BoxLengthZ/2.0D0)then
                              zlj = zlj - BoxLengthZ
                              IF(LJVectorZ0.GT.0.0)LJVectorZ = LJVectorZ0 - BoxLengthZ
                              IF(LJVectorZ0.LT.0.0)LJVectorZ = LJVectorZ0 + BoxLengthZ
                           endif    
                      endif
                
!                  -----------------------
!                  INTER-PARTICLE DISTANCE
!                  -----------------------
                      d2 = xlj*xlj + ylj*ylj + zlj*zlj
!                  ----------------
!                  POTENTIAL ENERGY
!                  ----------------
                     call potentialEnergy(d2, sigma, welldepth,U, V)
                        VectorProduct = MCVectorX*LJVectorX + MCVectorY*LJVectorY + MCVectorZ*LJVectorZ
                        V = V*VectorProduct/d2
                        if(SMediation1)then
                          EnergyCij = sqrt(EnergyCi*SFEnergy(j))
                          gSM = exp(-ksi*EnergyCij/Temperature)
                          U = U *gSM
                        endif
                        if(SMediation2)then
                          if(mcz(i).lt.SMLimitZ)U = U *ksiF
                          if(mcz(j).lt.SMLimitZ)U = U *ksiF
                        endif
                        
                        pEnergy = pEnergy + U
                        if(ExtraBin)then
                           if( (indNinSBox(i).lt.SubBox) .and. (indNinSBox(j).lt.SubBox) )Energy = Energy + U
                           if( (indNinSBox(i).eq.SubBox) .and. (indNinSBox(j).eq.SubBox) )ExtEnergy1 = ExtEnergy1 + U
                        else
                            Energy = Energy + U
                        endif  
                        Virial = Virial + V
                        IF(indNinSBox(j).eq.SubBox)VirialF = VirialF + V
                   enddo
                enddo 
                
           !!---------------------------------------------
           !! Coulombs force
           !!---------------------------------------------
                if(BUCKINGHAM)then  !! 3-if
                   CLEnergy = 0.0e0
                else
                   do k = 1, NCL
                      do m = 1, NCL
                         xcl =  abs(CLx(j,m) - CLx(i,k))
                         ycl =  abs(CLy(j,m) - CLy(i,k))
                         zcl =  abs(CLz(j,m) - CLz(i,k))
                         CLVectorX0 = CLx(j,m) - CLx(i,k)
                         CLVectorY0 = CLy(j,m) - CLy(i,k)
                         CLVectorZ0 = CLz(j,m) - CLz(i,k)
                         CLVectorX = CLVectorX0
                         CLVectorY = CLVectorY0
                         CLVectorZ = CLVectorZ0
                         if(pbcx)then                
                              if(xmclj > CarbonLengthx1/2.0D0)then
                                 xcl = xcl - CarbonLengthx1
                                 IF(CLVectorX0 .GT. 0.0)CLVectorX = CLVectorX0 - CarbonLengthx1
                                 IF(CLVectorX0 .LT. 0.0)CLVectorX = CLVectorX0 + CarbonLengthx1 
                              endif
                         endif
                         if(pbcy)then
                              if(ymclj > CarbonLengthy1/2.0D0)then
                                 ycl = ycl - CarbonLengthy1
                                 IF(CLVectorY0.GT.0.0)CLVectorY = CLVectorY0 - CarbonLengthy1
                                 IF(CLVectorY0.LT.0.0)CLVectorY = CLVectorY0 + CarbonLengthy1
                              endif
                         endif
                         if(pbcz)then
                              if(zmclj > BoxLengthZ/2.0D0)then
                                 zcl = zcl - BoxLengthZ
                                 IF(CLVectorZ0.GT.0.0)CLVectorZ = CLVectorZ0 - BoxLengthZ
                                 IF(CLVectorZ0.LT.0.0)CLVectorZ = CLVectorZ0 + BoxLengthZ
                              endif
                         endif
                         cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                         call coulombforce(cld, charges(k),charges(m),CLU,CLV)
                         VectorProduct = MCVectorX*CLVectorX + MCVectorY*CLVectorY + MCVectorZ*CLVectorZ
                         CLV = CLV*VectorProduct/cld**2.0
                         if(SMediation1)then
                           EnergyCij = sqrt(EnergyCi*SFEnergy(j))
                           gSM = exp(-ksi*EnergyCij/Temperature)
                           CLU = CLU *gSM
                         endif
                         if(SMediation2)then
                           if(mcz(i).lt.SMLimitZ)CLU = CLU *ksiF
                           if(mcz(j).lt.SMLimitZ)CLU = CLU *ksiF
                         endif
                         pEnergy  = pEnergy + CLU
                         if(ExtraBin)then
                           if( (indNinSBox(i).lt.SubBox) .and. (indNinSBox(j).lt.SubBox) )CLEnergy = CLEnergy + CLU
                           if( (indNinSBox(i).eq.SubBox) .and. (indNinSBox(j).eq.SubBox) )ExtCLEnergy1 = ExtCLEnergy1 + CLU
                         else
                            CLEnergy = CLEnergy + CLU
                            VirialCL = VirialCL + CLV
                            IF(indNinSBox(j).eq.SubBox)VirialF = VirialF + CLV
                         endif                        
                       enddo
                    enddo   
                 endif   !! 3-if   
              endif !! 4-if
                        IF(ADSORPTION.AND. LOCALF)THEN
                           if(flag.eq.0)then
                             EnergyMatrix(Npart+1,j) = pEnergy
                             EnergyMatrix(j,Npart+1) = pEnergy
                           else
                             EnergyMatrix(i,j) = pEnergy
                             EnergyMatrix(j,i) = pEnergy
                           endif
                        ENDIF
             endif  !! 5-if
            enddo  !! 1-do
          else  !! 2-if
            do j = 1, Npart
               if(j .NE. i)then
               Energy = Energy + EnergyMatrix(i,j)
               CLEnergy = 0.0e0
               endif
             enddo     
          endif !! 2-if
        endif   !! 1-if


     return
     end  subroutine energysingleparticle              

   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE CONFIGENERGY
!      THIS SUBROUTINE CALCULATES THE CONFIGURATION ENERGY OF N PARTICLES 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ConfigEnergy(Energy,CLEnergy,EnergyC)
      USE FLUID_SOLID_M
      USE PHYSICAL_M
      USE MCSETTING_M
      USE BUCKINGHAM_M
      USE SUBBOX_M
      implicit none
 
!     ------------------------
!     DECLARATION OF VARIABLES
!     ------------------------
      integer i, j, N
      integer k, m, flag, FLAG1, FLAG2
      
      real*8 sigma, welldepth, sigmaSF,welldepthSF                 
      real*8 Energy,CLEnergy,EnergyC 
      real*8 UC,U,CLU,BU,CLV
      real*8 xlj, ylj, zlj, d2
      real*8 xcl,ycl,zcl,cld
      real*8 xclj, yclj, zclj, dc2
      real*8 xmclj,ymclj,zmclj
      real*8 xcmcij,ycmcij,zcmcij
      real*8 distanceZ
      real*8 xmcn,ymcn,zmcn,mcd2
      real*8 xcmcn,ycmcn,zcmcn,cmcd2
      real*8 deltay,deltay1,cyi
      real*8 smcy,LJy1(20), sLJy(20)
      real*8 EnergyCij, gSM
      real*8 PEnergy, PCLEnergy, PEnergyC
      real*8 MCVectorX, MCVectorY, MCVectorZ, MCVectorX0, MCVectorY0, MCVectorZ0
      real*8 LJVectorX, LJVectorY, LJVectorZ, LJVectorX0, LJVectorY0, LJVectorZ0
      real*8 CLVectorX, CLVectorY, CLVectorZ, CLVectorX0, CLVectorY0, CLVectorZ0
      real*8 VectorProduct, PairV, V

!     -----------------------------------
!     INITIALIZATION OF ENERGY 
!     -----------------------------------
       Energy      = 0.0E0
       CLEnergy    = 0.0E0
       ExtEnergy   = 0.0E0
       ExtCLEnergy = 0.0E0
       EnergyC     = 0.0E0
       Virial      = 0.0E0
       VirialC     = 0.0E0
       BulkVirial  = 0.0E0
       N        = Npart  
       EinBulk  = 0.0e0 
       SFEnergy = 0.0
       if(N.eq.0) return

!     --------------------------------------------------------
!     SUMMING THE PAIRWISE INTERACTION ENERGIES FOR ADSORPTION
!     --------------------------------------------------------
        do i=1, N-1
           FLAG1=indNinSBox(i)
           do j=i+1, N
              FLAG2=indNinSBox(j)
              PairV = 0.0
              if(ExtraBin)then
                if( (indNinSBox(i).eq.SubBox) .AND. (indNinSBox(j).lt.SubBox) )CYCLE
                if( (indNinSBox(i).lt.SubBox) .AND. (indNinSBox(j).eq.SubBox) )CYCLE
              endif
              xmclj = abs(mcx(i) - mcx(j))
              ymclj = abs(mcy(i) - mcy(j))
              zmclj = abs(mcz(i) - mcz(j))
              xmcn = xmclj
              ymcn = ymclj
              zmcn = zmclj
              MCVectorX0 = mcx(j) - mcx(i)
              MCVectorY0 = mcy(j) - mcy(i)
              MCVectorZ0 = mcz(j) - mcz(i)
              MCVectorX = MCVectorX0
              MCVectorY = MCVectorY0
              MCVectorZ = MCVectorZ0
                 if(pbcx)then
                      if(xmclj > CarbonLengthx1/2.0D0)then
                         xmcn = xmclj - CarbonLengthx1
                         IF(MCVectorX0.GT.0.0)MCVectorX = MCVectorX0 - CarbonLengthx1
                         IF(MCVectorX0.LT.0.0)MCVectorX = MCVectorX0 + CarbonLengthx1
                      endif
                 endif
                 if(pbcy)then
                      if(ymclj > CarbonLengthy1/2.0D0)then
                         ymcn = ymclj - CarbonLengthy1
                         IF(MCVectorY0.GT.0.0)MCVectorY = MCVectorY0 - CarbonLengthy1
                         IF(MCVectorY0.LT.0.0)MCVectorY = MCVectorY0 + CarbonLengthy1
                      endif
                 endif
                 if(pbcz)then
                      if(zmclj > BoxLengthZ/2.0D0)then
                         zmcn = zmclj - BoxLengthZ
                         IF(MCVectorZ0.GT.0.0)MCVectorZ = MCVectorZ0 - BoxlengthZ
                         IF(MCVectorZ0.LT.0.0)MCVectorZ = MCVectorZ0 + BoxlengthZ 
                      endif
                 endif 
                 mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
              if(mcd2.gt.r2cutoff)then
                 CYCLE
              else
                 PEnergy = 0.0E0
                 PCLEnergy = 0.0E0
                 do k=1, NLJ                      ! the LJ sites of particle i
                    do m = 1, NLJ                 ! the LJ sites of particle j 
                       sigma = (sigmaFF(k)+sigmaFF(m))/2.0e0
                       welldepth = sqrt(welldepthFF(k)*welldepthFF(m))  
!                 -------------------
!                 DISTANCE IN X, Y, Z
!                 -------------------
                       xlj = abs(ljx(j,m) - ljx(i,k))
                       ylj = abs(ljy(j,m) - ljy(i,k))
                       zlj = abs(ljz(j,m) - ljz(i,k))
                       LJVectorX0 = ljx(j,m) - ljx(i,k)
                       LJVectorY0 = ljy(j,m) - ljy(i,k)
                       LJVectorZ0 = ljz(j,m) - ljz(i,k)
                       LJVectorX = LJVectorX0
                       LJVectorY = LJVectorY0
                       LJVectorZ = LJVectorZ0
!                ----------------------------------------------------
!                ENSURE THAT WE DEAL WITH THE NEAREST PERIODIC IMAGES
!                ----------------------------------------------------
                       if(pbcx)then
                            if(xmclj > CarbonLengthx1/2.0D0)then
                               xlj = xlj - CarbonLengthx1
                               IF(LJVectorX0 .GT. 0.0)LJVectorX = LJVectorX0 - CarbonLengthx1
                               IF(LJVectorX0 .LT. 0.0)LJVectorX = LJVectorX0 + CarbonLengthx1
                            endif
                       endif
                       if(pbcy)then
                            if(ymclj > CarbonLengthy1/2.0D0)then 
                               ylj = ylj - CarbonLengthy1
                               IF(LJVectorY0 .GT. 0.0)LJVectorY = LJVectorY0 - CarbonLengthy1
                               IF(LJVectorY0 .LT. 0.0)LJVectorY = LJVectorY0 + CarbonLengthy1
                            endif
                       endif
                       if(pbcz)then
                            if(zmclj > BoxLengthZ/2.0D0)then
                               zlj = zlj - BoxLengthZ
                               IF(LJVectorZ0.GT.0.0)LJVectorZ = LJVectorZ0 - BoxlengthZ
                               IF(LJVectorZ0.LT.0.0)LJVectorZ = LJVectorZ0 + BoxlengthZ
                            endif
                       endif 
!                -----------------------
!                INTER-PARTICLE DISTANCE
!                -----------------------
                       d2 = xlj*xlj + ylj*ylj + zlj*zlj
!                ----------------
!                POTENTIAL ENERGY
!                ----------------
                      call potentialEnergy(d2,sigma,welldepth,U,V)
                        VectorProduct = MCVectorX*LJVectorX + MCVectorY*LJVectorY + MCVectorZ*LJVectorZ
                        V = V*VectorProduct/d2
                        PairV = PairV + V
                        if(SMediation1)then
                           EnergyCij = sqrt(SFenergy(i)*SFenergy(j))
                           gSM = exp(-ksi*EnergyCij/Temperature)
                           U = U * gSM
                         endif
                         if(SMediation2)then
                          if(mcz(i).lt.SMLimitZ)U = U *ksiF
                          if(mcz(j).lt.SMLimitZ)U = U *ksiF
                         endif
                       PEnergy = PEnergy + U  
                       if(ExtraBin)then
                         if( (indNinSBox(i).lt.SubBox) .and. (indNinSBox(j).lt.SubBox) )Energy = Energy + U
                         if( (indNinSBox(i).eq.SubBox) .and. (indNinSBox(j).eq.SubBox) )ExtEnergy = ExtEnergy + U
                       else
                         Energy = Energy + U
                       endif
                   enddo
                enddo
             !!---------------------------
             !! COULOMB FORCE
             !!---------------------------                    
                do k = 1, NCL
                   do m = 1, NCL
                    xcl =  abs(CLx(j,m) - CLx(i,k))
                    ycl =  abs(CLy(j,m) - CLy(i,k))
                    zcl =  abs(CLz(j,m) - CLz(i,k))
                    CLVectorX0 = CLx(j,m) - CLx(i,k)
                    CLVectorY0 = CLy(j,m) - CLy(i,k)
                    CLVectorZ0 = CLz(j,m) - CLz(i,k)
                    CLVectorX = CLVectorX0
                    CLVectorY = CLVectorY0
                    CLVectorZ = CLVectorZ0
                    if(pbcx)then                
                         if(xmclj > CarbonLengthx1/2.0D0)then
                            xcl = xcl - CarbonLengthx1
                            IF(CLVectorX0.GT.0.0)CLVectorX = CLVectorX0 - CarbonLengthx1
                            IF(CLVectorX0.LT.0.0)CLVectorX = CLVectorX0 + CarbonLengthx1
                         endif
                    endif
                    if(pbcy)then
                         if(ymclj > CarbonLengthy1/2.0D0)then
                            ycl = ycl - CarbonLengthy1
                            IF(CLVectorY0.GT.0.0)CLVectorY = CLVectorY0 - CarbonLengthy1
                            IF(CLVectorY0.LT.0.0)CLVectorY = CLVectorY0 + CarbonLengthy1
                         endif
                    endif
                    if(pbcz)then
                         if(zmclj > BoxLengthZ/2.0D0)then
                            zcl = zcl - BoxLengthZ
                            IF(CLVectorZ0.GT.0.0)CLVectorZ = CLVectorZ0 - BoxlengthZ
                            IF(CLVectorZ0.LT.0.0)CLVectorZ = CLVectorZ0 + BoxlengthZ
                         endif
                    endif
                    cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                    call coulombforce(cld, charges(k),charges(m),CLU,CLV)
                        VectorProduct = MCVectorX*CLVectorX + MCVectorY*CLVectorY + MCVectorZ*CLVectorZ
                        CLV = CLV*VectorProduct/cld**2.0
                        PairV = PairV + CLV
                        if(SMediation1)then
                           EnergyCij = sqrt(SFenergy(i)*SFenergy(j))
                           gSM = exp(-ksi*EnergyCij/Temperature)
                           CLU = CLU * gSM
                         endif
                         if(SMediation2)then
                          if(mcz(i).lt.SMLimitZ)CLU = CLU *ksiF
                          if(mcz(j).lt.SMLimitZ)CLU = CLU *ksiF
                         endif
                       if(ExtraBin)then
                         if( (indNinSBox(i).lt.SubBox) .and. (indNinSBox(j).lt.SubBox) )CLEnergy = CLEnergy + CLU
                         if( (indNinSBox(i).eq.SubBox) .and. (indNinSBox(j).eq.SubBox) )ExtCLEnergy = ExtCLEnergy + CLU
                       else
                         CLEnergy = CLEnergy + CLU
                       endif

                      PCLEnergy = PCLEnergy + CLU  
                   enddo
                enddo 
                Virial = Virial + PairV
              endif 
              
              if(NVT.AND.SDisplace.AND. (.NOT. ExtraBin))then
                  if( (FLAG2.eq.SubBox) .AND. (FLAG1.eq.SubBox) )then
                      EinBulk = EinBulk + PEnergy + PCLEnergy
                      BulkVirial = BulkVirial + PairV
                  elseif((FLAG2.eq.SubBox) .AND. (FLAG1.ne.SubBox))then
                      EinBulk = EinBulk + (PEnergy + PCLEnergy)/2.0
                      BulkVirial = BulkVirial + PairV/2.0
                  elseif((FLAG2.ne.SubBox) .AND. (FLAG1.eq.SubBox))then
                      EinBulk = EinBulk + (PEnergy + PCLEnergy)/2.0
                      BulkVirial = BulkVirial + PairV/2.0
                  endif
              endif
                
         enddo    
       enddo
!     ------------------------------------------------------
!     THE PAIRWISE INTERACTION BETWEEN FLUID AND SOLID
!     ------------------------------------------------------
       if(steele)then
          do i = 1, N
             if(ExtraBin .AND. (indNinSBox(i).eq.SubBox))CYCLE
                call steelepotential(i,UC)
                EnergyC = EnergyC + UC 
                SFEnergy(i) = UC
!                VirialC = VirialC + VC
                if(NVT.AND.SDisplace .AND. (.NOT. ExtraBin).AND. (UC.NE.0.0))then
                  if(indNinSBox(i).eq.SubBox)then
                     EinBulk = EinBulk + UC
                  endif
                endif
          enddo
       endif
       
       if(CrowellChang)then
          do i = 1, N
             call CrowellChangPotential(i,UC)
             EnergyC = EnergyC + UC 
             SFEnergy(i) = UC
          enddo
       endif

       IF(Bojan)then
         do i = 1, N
            if(ExtraBin .AND. (indNinSBox(i).eq.SubBox))CYCLE
            PEnergyC = 0.0
           do j = 1, NS 
              call BojanPotential(j, i,BU) 
              EnergyC = EnergyC + BU 
              PEnergyC = PEnergyC + BU
              if(BoPBCy)then
                   deltay  = abs(mcy(i)  - (StripC(j)+ StripW(j)/2.0e0)) 
                   deltay1 = abs(mcy(i)  - (StripC(j)- StripW(j)/2.0e0))  
                   if(deltay .lt. deltay1)  deltay = deltay1
                   if(deltay .gt. (carbonlengthy1/2.0e0))then
                       if((mcy(i).ge.0.0e0) .and. (mcy(i).le.CarbonLengthy1/2.0e0))then
                          cyi = mcy(i)  + carbonlengthy1
                          do k = 1, NLJ
                             LJy1(k) = LJy(i,k) + carbonlengthy1
                          enddo
                       else
                          cyi = mcy(i)  - carbonlengthy1
                          do k = 1, NLJ
                             LJy1(k) = LJy(i,k) - carbonlengthy1
                          enddo  
                       endif
                       smcy = mcy(i)
                       mcy(i) = cyi
                       do k = 1, NLJ
                          sLJy(k) = LJy(i,k)
                          LJy(i,k) = LJy1(k)
                       enddo
                       call BojanPotential(j, i,BU) 
                       EnergyC = EnergyC + BU
                       PEnergyC = PEnergyC + BU 
                       VirialC = VirialC + VC
                       mcy(i) = smcy 
                       do k = 1, NLJ
                          LJy(i,k) = sLJy(k)
                       enddo 
                   endif
              endif               
           enddo ! do j
           IF(inclinePore)then
              call BojanPotential(NS+1, i,BU) 
              EnergyC = EnergyC + BU 
              PEnergyC = PEnergyC + BU
           ENDIF
           IF(Close1End.OR.Close2Ends)then
              call BojanPotential(NS+2, i,BU) 
              EnergyC = EnergyC + BU 
              PEnergyC = PEnergyC + BU
           ENDIF
           SFEnergy(i)=PEnergyC
         enddo ! do i
       ENDIF   
      
      DO i=1, N
            PEnergyC = 0.0
         do j=1, NC3
            xcmcij = abs(xc(j)-mcx(i))
            ycmcij = abs(yc(j)-mcy(i))
            zcmcij = abs(zc(j)-mcz(i))
            xcmcn  = xcmcij
            ycmcn  = ycmcij
            zcmcn  = zcmcij 
            if(pbcx)then
               if(xcmcij > CarbonLengthx1/2.0D0) xcmcn = xcmcn - CarbonLengthx1
            endif
            if(pbcy)then
               if(ycmcij > CarbonLengthy1/2.0D0) ycmcn = ycmcn - CarbonLengthy1
            endif
            if(pbcz)then
               if(zcmcij > BoxLengthZ/2.0D0) zcmcn = zcmcn - BoxLengthZ
            endif
            cmcd2 = xcmcn*xcmcn + ycmcn*ycmcn + zcmcn*zcmcn
            if(cmcd2.gt.r2cutoff)then
              CYCLE
            else
              do k = 1, NLJ 
                 sigmaSF = (sigmaFF(k)+sigmaSS)/2.0e0
                 welldepthSF = sqrt(welldepthFF(k)*welldepthSS) 
!               -------------------
!               DISTANCE IN X, Y, Z
!               -------------------
                 xclj = abs(xc(j) - ljx(i,k))
                 yclj = abs(yc(j) - ljy(i,k))
                 zclj = abs(zc(j) - ljz(i,k))
!               ----------------------------------------------------
!               ENSURE THAT WE DEAL WITH THE NEAREST PERIODIC IMAGES
!               ----------------------------------------------------
                if(pbcx)then
                  if(xcmcij > CarbonLengthx1/2.0D0) xclj = xclj - CarbonLengthx1
                endif
                if(pbcy)then
                  if(ycmcij > CarbonLengthy1/2.0D0) yclj = yclj - CarbonLengthy1
                endif
                if(pbcz)then
                   if(zcmcij > BoxLengthZ/2.0D0) zclj = zclj - BoxLengthZ
                endif
!               -----------------------
!               INTER-PARTICLE DISTANCE
!               -----------------------
                dc2 = xclj*xclj + yclj*yclj + zclj*zclj
                
                call potentialEnergyC(dc2,sigmaSF,welldepthSF,UC)
                EnergyC = EnergyC + UC
                PEnergyC = PEnergyC + UC
                VirialC = VirialC + VC
              enddo  ! do k
            endif
          enddo ! do j  
          SFEnergy(i)= PEnergyC
        ENDDO ! do i
    
    return
    end  
    
    
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE POTENTIALENERGY
!      THIS SUBROUTINE CALCULATE THE POTENTIALENERGY AMONG THE N PARTICLES
!      USE LENNARD-JONES 12-6 POTENTIAL MODEL 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine potentialEnergy(d2,sigma1, welldepth1, U,V)
      USE Constant_M
      USE MCSETTING_M
      USE  BUCKINGHAM_M
      implicit none
      
      real*8 d1, d2, d6, U,V
      real*8 sigma1, welldepth1
      real*8 term1, term2, term3
      
      if(.NOT. BUCKINGHAM)then      
!          ----------------------------------
!          LENNARD-JONES 12-6 POTENTIAL MODEL
!          ----------------------------------
           d6  = (sigma1**2/d2)**3
           U   =  4.0E0*welldepth1*d6*(d6 - 1.0E0)
           V   =  16.0E0*welldepth1*d6*(d6 - 0.5E0)
       else
!          --------------------------------
!          BUCKINGHAM EXP-6 POTENTIAL MODEL
!          --------------------------------
           d1 = sqrt(d2)
           if(d1 .le. rMin)d1 = rMin
           term1 = exp(alpha_B*(1.0E0-d1/rM))
           term2 = (rM/d1)**6
           term3 = 1.0E0/(1.0E0 - 6.0E0/alpha_B)
           
           U = term3*(6.0e0*term1/alpha_B - term2)
           V = 2.0e0*term3*(d1*term1/rM - term2)
       endif       
     
     return
     end
      
  
  
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE POTENTIALENERGYC
!      THIS SUBROUTINE CALCULATES THE POTENTIALENERGY BETWEEN THE PARTICLES AND
!      THE CARBON ATOMS, USE LENNARD-JONES 12-6 POTENTIAL MODEL 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine potentialEnergyC(dc2,sigmaSF1,welldepthSF1,UC)
      USE Constant_M
      USE MCSETTING_M
      implicit none
      
      real*8 dc2, dc6, UC
      real*8 sigmaSF1, WellDepthSF1
      
      
!          ----------------------------------
!          LENNARD-JONES 12-6 POTENTIAL MODEL
!          ----------------------------------
!           if(dc2.gt.r2CutOff)then
!                 UC   = 0.E0
!                 VC   = 0.E0
!           else
                 dc6  = ((sigmaSF1**2.0E0)/dc2)**3
                 UC   =  4.0E0*WellDepthSF1*dc6*(dc6 - 1.0E0)
                 VC   =  16.0E0*WellDepthSF1*dc6*(dc6 - 0.5E0)
!          endif
     
     return
     end   
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE MCMOVE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine mcMove(TEnergy,Energy,CLEnergy,EnergyC, numberOfMove, numberOfAcceptanceMove)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer N, i,J,new
       integer numberOfMove, numberOfAcceptanceMove
       integer bin, bin0
       integer flag, fSub, tSub
       
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       real*8 TEnergy, TEnergyOld, TEnergyNew
       real*8 Energy,EnergyC
       real*8 deltaEnergy, rdn
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 ECLxnew(20),ECLynew(20),ECLznew(20) 
       real*8 CLEnergy    
             

!      -----------------------------------
!      INCREMENT THE NUMBER OF MOVE BY ONE
!      -----------------------------------
       if(Npart.eq.0)return
       numberOfMove = numberOfMove + 1
       N    = Npart
!      ----------------------------------------------------------------------------
!      RANDOMLY SELECT A PARTICLE & OBTAIN THE ENERGY OF THAT PARTICLES WITH OTHERS
!      ----------------------------------------------------------------------------
       call random_number(rdn)
       i  = int(rdn*N) + 1
!      ------------------------------------------------------------------------------
!      THIS IS TO ENSURE THAT i SHOULD NOT EXCEED N, AS THIS CAN OCCURS WHEN rand()=1
!      ------------------------------------------------------------------------------
       if(i.gt.N) i = N      
       
!          --------------------------------------
!          OLD ENERGY OF PARTICLE I IN ADSORPTION
!          --------------------------------------
           call energySingleParticle(0,i, EnergyOld, CLEnergyOld, EnergyCOld )

!          --------------------------------
!          SUM OF EnergyOld and EnergyCOld 
!          --------------------------------
           TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld
!          ----------------------------------------------------
!          DISPLACE THE POSITION OF THE SELECTED PARTICLE i 
!          AND MAKE SURE THAT NEW POSITIONS ARE INSIDE THE BOX
!          ----------------------------------------------------
           new = N + 1
           mcx(new) = mcx(i)
           mcy(new) = mcy(i) 
           mcz(new) = mcz(i)
           do j = 1, NLJ
              LJx(new,j) = LJx(i,j)
              LJy(new,j) = LJy(i,j)
              LJz(new,j) = LJz(i,j)
           enddo
       
           do j = 1, NCL
              CLx(new,j) = CLx(i,j)
              CLy(new,j) = CLy(i,j)
              CLz(new,j) = CLz(i,j)
              ECLx(new,j) = ECLx(i,j)
              ECLy(new,j) = ECLy(i,j)
              ECLz(new,j) = ECLz(i,j)
           enddo
           
           do j = 1, Npart
              if(j.ne.i)then
                 if(ADSORPTION.AND. LOCALF)THEN
                 EnergyMatrix(new,j) = EnergyMatrix(i,j)
                 EnergyMatrix(j,new) = EnergyMatrix(j,i)
                 ENDIF
              endif
           enddo
           
           if(ADSORPTION.AND. LOCALF)THEN
                 EnergyMatrix(new,new) = EnergyMatrix(i,i)
            ENDIF
       
           call random_number(rdn)
           mcx(i) = mcx(i) + (rdn - 0.5E0)*deltaRx
           if(.NOT. pbcx)then
             if((mcx(i) .lt. 0.E0) .or. (mcx(i) .gt. CarbonLengthx1))then
                mcx(i)=mcx(new)
                return
             endif 
           endif
       
           call random_number(rdn)
           mcy(i) = mcy(i) + (rdn - 0.5E0)*deltaRy  
           if(.NOT. pbcy)then    
             if((mcy(i) .lt. 0.E0) .or. (mcy(i) .gt. CarbonLengthy1))then
                mcx(i) = mcx(new)
                mcy(i) = mcy(new)
                return
             endif
           endif
        
           call random_number(rdn)
           mcz(i) = mcz(i) + (rdn - 0.5E0)*deltaRz
           if(.NOT. pbcz)then
              if((mcz(i) .le. 0.E0) .or. (mcz(i) .ge. BoxLengthZ))then
                 mcx(i) = mcx(new)
                 mcy(i) = mcy(new)
                 mcz(i) = mcz(new)
               return
              endif
           endif
            
            
           if(Clo1Hard .and. (mcz(i) .le. 0.E0))then
                 mcx(i) = mcx(new)
                 mcy(i) = mcy(new)
                 mcz(i) = mcz(new)
                return
           endif
     
             
           if(   InclinePore .and.   &
           &     mcz(i) .gt. IncCornerZ  .and.   &
           &     mcy(i) .lt. ((mcz(i)-IncCornerZ)/tan(IncAngle))  )then
                 mcx(i) = mcx(new)
                 mcy(i) = mcy(new)
                 mcz(i) = mcz(new)
                 return       
           endif
       
           if(   SurfPore .and.   &
           &     mcz(i) .gt. CriLengthZ  .and.   &
           &     mcy(i) .lt. CarbonLengthy0  .and. & 
           &     mcy(i) .lt. ((mcz(new)-CriLengthZ)/tan(surftheta))  )then
                 mcx(i) = mcx(new)
                 mcy(i) = mcy(new)
                 mcz(i) = mcz(new)
              return       
           endif
       
 !          if(Bojan)then
 !              do j = 1, NS
 !                  if( StripZ0(j) .gt. 0.0e0 )then
 !                      if(   mcz(i) .lt. StripZ0(j) .and. &
 !                          &   mcy(i) .gt. (StripC(j)- StripW(j)/2.0e0)  .and. &
 !                          &   mcy(i) .lt. (StripC(j)+ StripW(j)/2.0e0) ) then
 !                          mcx(i) = mcx(new)
 !                          mcy(i) = mcy(new)
 !                          mcz(i) = mcz(new)
 !                          return  
 !                      endif
 !                  endif
 !              enddo
 !          endif 
           
           if(BojanSlit)then      
                 do j = 1, NS
                     if( mcy(i).gt.StripLy(j) .AND. mcy(i).le.StripHy(j)  )then
                         if(   StripAngleY(j) .GT. 0.0 .and. &
             &                 mcz(i) .lt. InLzLyStripZ(j) .and. &
             &                (InLzLyStripZ(j)-mcz(i))   .GE. (mcy(i)-StripLy(j))*tan(StripAngleY(j))   ) then 
                               mcx(i) = mcx(new)
                               mcy(i) = mcy(new)
                               mcz(i) = mcz(new)
                               return  
                       elseif(   StripAngleY(j) .GT. 0.0 .and. &
             &                 mcz(i) .gt. InHzLyStripZ(j) .and. &
             &                 (mcz(i)-InHzLyStripZ(j))   .GE. (mcy(i)-StripLy(j))*tan(StripAngleY(j))   ) then 
                               mcx(i) = mcx(new)
                               mcy(i) = mcy(new)
                               mcz(i) = mcz(new)
                               return 
                      elseif(   StripAngleY(j) .LT. 0.0 .and. &
             &                 mcz(i) .lt. InLzLyStripZ(j) .and. &
             &                (InLzLyStripZ(j)-mcz(i))   .GE. (mcy(i)-StripHy(j))*tan(StripAngleY(j))   ) then 
                               mcx(i) = mcx(new)
                               mcy(i) = mcy(new)
                               mcz(i) = mcz(new)
                               return  
                      elseif(   StripAngleY(j) .LT. 0.0 .and. &
             &                 mcz(i) .gt. InHzLyStripZ(j) .and. &
             &                 (mcz(i)-InHzLyStripZ(j))   .GE. (mcy(i)-StripHy(j))*tan(StripAngleY(j))   ) then 
                               mcx(i) = mcx(new)
                               mcy(i) = mcy(new)
                               mcz(i) = mcz(new)
                               return 
                      elseif(   StripAngleY(j) .EQ. 0.0 .and. &
             &                  ( mcz(i) .lt. InLzLyStripZ(j) .or. mcz(i) .gt. InHzLyStripZ(j) ) ) then 
                               mcx(i) = mcx(new)
                               mcy(i) = mcy(new)
                               mcz(i) = mcz(new)
                               return
                         endif
                   endif
                enddo
             endif    
                       
       
           call positioncheck(mcx(i),mcy(i),mcz(i))  
      
!    ---------------------------------------
!     FOR ORIENTATIONAL MOVES    
!    ---------------------------------------
!     Put the mass center on original point
!    ---------------------------------------
          do j = 1, NLJ
             LJxnew(j) = LJx(i,j) - MCx(new)
             LJynew(j) = LJy(i,j) - MCy(new)
             LJznew(j) = LJz(i,j) - MCz(new)
          enddo
          do j = 1, NCL
             CLxnew(j) = CLx(i,j) - MCx(new)
             CLynew(j) = CLy(i,j) - MCy(new)
             CLznew(j) = CLz(i,j) - MCz(new)
             ECLxnew(j) = ECLx(i,j) - MCx(new)
             ECLynew(j) = ECLy(i,j) - MCy(new)
             ECLznew(j) = ECLz(i,j) - MCz(new)
          enddo
!    ----------
!     Rotate
!    ----------
          call ROTATE
      
          do j = 1, NLJ
             LJx(i,j) = Rot(1)*LJxnew(j)+ Rot(4)*LJynew(j) + Rot(7)*LJznew(j) + MCx(i)
             LJy(i,j) = Rot(2)*LJxnew(j)+ Rot(5)*LJynew(j) + Rot(8)*LJznew(j) + MCy(i)
             LJz(i,j) = Rot(3)*LJxnew(j)+ Rot(6)*LJynew(j) + Rot(9)*LJznew(j) + MCz(i)
          enddo

          do j = 1, NCL
             CLx(i,j) = Rot(1)*CLxnew(j)+ Rot(4)*CLynew(j) + Rot(7)*CLznew(j) + MCx(i)
             CLy(i,j) = Rot(2)*CLxnew(j)+ Rot(5)*CLynew(j) + Rot(8)*CLznew(j) + MCy(i)
             CLz(i,j) = Rot(3)*CLxnew(j)+ Rot(6)*CLynew(j) + Rot(9)*CLznew(j) + MCz(i)
             ECLx(i,j) = Rot(1)*ECLxnew(j)+ Rot(4)*ECLynew(j) + Rot(7)*ECLznew(j) + MCx(i)
             ECLy(i,j) = Rot(2)*ECLxnew(j)+ Rot(5)*ECLynew(j) + Rot(8)*ECLznew(j) + MCy(i)
             ECLz(i,j) = Rot(3)*ECLxnew(j)+ Rot(6)*ECLynew(j) + Rot(9)*ECLznew(j) + MCz(i)
          enddo     

!          -------------------------
!          NEW ENERGY OF PARTICLE I
!          -------------------------
           call  energySingleParticle(1,i, EnergyNew, CLEnergyNew, EnergyCNew)
!          -----------------------------
!          SUM OF ENERGY AT NEW POSITION
!          -----------------------------
           TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
       
!          ------------------------------------------------
!          CHECK TO SEE WHETHER THE MOVE SHOULD BE ACCEPTED
!          ------------------------------------------------
           deltaEnergy =  (TEnergyNew - TEnergyOld)/Temperature
   
!           ------------------------------------------
!           ACCEPT THE MOVE IF deltaEnergy IS NEGATIVE
!           ------------------------------------------
            call random_number(rdn)
            if(deltaEnergy.gt.75.D0)then             
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  mcz(i) = mcz(new)
                  do j = 1, NLJ
                     ljx(i,j) = ljx(new,j)
                     ljy(i,j) = ljy(new,j)
                     ljz(i,j) = ljz(new,j)
                  enddo
                  do j = 1, NCL
                     CLx(i,j) = CLx(new,j)
                     CLy(i,j) = CLy(new,j)
                     CLz(i,j) = CLz(new,j)
                     ECLx(i,j) = ECLx(new,j)
                     ECLy(i,j) = ECLy(new,j)
                     ECLz(i,j) = ECLz(new,j)

                  enddo
                   do j = 1, Npart
                      if(j.ne.i)then
                        IF(Adsorption.AND. LOCALF)THEN
                        EnergyMatrix(i,j) = EnergyMatrix(new,j)
                        EnergyMatrix(j,i) = EnergyMatrix(new,j)
                        ENDIF
                       endif
                    enddo
                    IF(Adsorption.AND. LOCALF)EnergyMatrix(i,i) = EnergyMatrix(new,new)
            elseif(deltaEnergy .lt. 0.D0 .OR. rdn.lt.exp(-deltaEnergy))then
                   numberOfAcceptanceMove = numberOfAcceptanceMove + 1
                   SFenergy(i) = EnergyCNew
                   Energy   = Energy  + EnergyNew - EnergyOld 
                   CLEnergy = CLEnergy + CLEnergyNew - CLEnergyOld
                   EnergyC  = EnergyC + EnergyCNew - EnergyCOld               
                   TEnergy  = TEnergy + TEnergyNew - TEnergyOld
                   call particleposition(new,-1)
                   call particleposition(i,1)
                   flag = 0
                   IF(SDisplace)THEN
                      DO WHILE(flag.eq.0)
                        do j=1, SubBox
                           if( HorizontalBin .AND. (MCz(new).GE.LayerZL(j)) .AND. ((MCz(new).LT.LayerZH(j))))then
                              fSub = j
                              flag = 1
                              EXIT
                           elseif(VerticalBin .AND. (MCy(new).GE.LayerYL(j)) .AND. ((MCy(new).LT.LayerYH(j))))then
                              fSub = j
                              flag = 1
                              EXIT
                           endif
                        enddo
                      ENDDO
                   
                      flag = 0
                      DO WHILE(flag.eq.0)
                        do j=1, SubBox
                           if( HorizontalBin .AND. (MCz(i).GE.LayerZL(j)) .AND. ((MCz(i).LT.LayerZH(j))))then
                              tSub = j
                              flag = 1
                              EXIT 
                           elseif(VerticalBin .AND. (MCy(i).GE.LayerYL(j)) .AND. ((MCy(i).LT.LayerYH(j))))then 
                              tSub = j
                              flag = 1   
                              EXIT
                           endif
                        enddo
                      ENDDO
                   
                      if(fSub.NE.tSub)then
                          IF(NVT)THEN
                             SeqNinSBox(fSub,i)=SeqNinSBox(fSub,int(NinSBox(fSub)))
                             NinSBox(fSub) =  NinSBox(fSub) - 1
                             NinSBox(tSub) =  NinSBox(tSub) + 1
                             SeqNinSBox(tSub,int(NinSBox(tSub))) = i
                             indNinSBox(i) = tSub
                             if(fSub.eq.SubBox)then
                                call EBulkChange(1,new,i)
                             elseif(tSub.eq.SubBox)then
                                call EBulkChange(2,new,i)
                             endif
                          ENDIF
                         SucMoveout(fSub) = SucMoveout(fSub)+1
                         SucMovein(tSub)  = SucMovein(tSub)+1
                      endif  
                      if(NVT.AND.(fSub.eq.tSub).AND. (fSub.eq.SubBox))then
                            call EBulkChange(3,new,i)
                      endif
                      if(NVT.AND.(fSub.eq.tSub).AND. (fSub.ne.SubBox))then
                            call EBulkChange(4,new,i)
                      endif


                   ENDIF 
!                   EnergyMatrix(i,i)=EnergyCNew
!                  WRITE(*,*)'2585 after move Energy,npart',Energy,npart
                   IF(Adsorption.AND. LOCALF)THEN 
                       bin0 = I_POS(i) 
                       call BinJudge(MCx(i),MCy(i),MCz(i),bin)      
                       do j = 1,NinBin(bin0)
                          if(SeqNinBin(bin0,j).eq.i)then
                            SeqNinBin(bin0,j)=SeqNinBin(bin0,NinBin(bin0))
                            EXIT
                          endif
                       enddo    
                       NinBin(bin0)=NinBin(bin0)-1
                       do j = 1, NinBin(bin0)
                            UFFinBin(bin0)=UFFinBin(bin0)-EnergyMatrix(new,SeqNinBin(bin0,j))
                       enddo
                       USFinBin(bin0)=USFinBin(bin0)-EnergyMatrix(new,new)
!                       do j = 1, Npart
!                          if(j.ne.i)then
!                            EnergyMatrix(i,j)= EnergyMatrix(Npart+1,j)
!                            EnergyMatrix(j,i)= EnergyMatrix(Npart+1,j)
!                          endif
!                       enddo
                           EnergyMatrix(i,i)=EnergyCNew
                       do j = 1, NinBin(bin)
                          UFFinBin(bin) = UFFinBin(bin) + EnergyMatrix(i,SeqNinBin(bin,j))
                       enddo
                          USFinBin(bin) = USFinBin(bin) + EnergyMatrix(i,i)
                          NinBin(bin) = NinBin(bin)+1
                          SeqNinBin(bin,NinBin(bin)) = i 
                          I_POS(i)  = bin
                   ENDIF
             else             
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  mcz(i) = mcz(new)
                  do j = 1, NLJ
                     ljx(i,j) = ljx(new,j)
                     ljy(i,j) = ljy(new,j)
                     ljz(i,j) = ljz(new,j)
                  enddo
                  do j = 1, NCL
                     CLx(i,j) = CLx(new,j)
                     CLy(i,j) = CLy(new,j)
                     CLz(i,j) = CLz(new,j)
                     ECLx(i,j) = ECLx(new,j)
                     ECLy(i,j) = ECLy(new,j)
                     ECLz(i,j) = ECLz(new,j)
                  enddo
                   do j = 1, Npart
                      if(j.ne.i)then
                        IF(Adsorption.AND. LOCALF)THEN
                        EnergyMatrix(i,j) = EnergyMatrix(new,j)
                        EnergyMatrix(j,i) = EnergyMatrix(new,j)
                        ENDIF
                       endif
                    enddo
                   IF(Adsorption.AND. LOCALF)EnergyMatrix(i,i) = EnergyMatrix(new,new)
             endif
      return
      endsubroutine mcMove 
      
      
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       SUBROUTINE ADJUSTMENT
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine adjustment(iEqui,numberOfMove, numberOfAcceptanceMove,AcceptanceRatio)
        USE PHYSICAL_M
        USE MCSETTING_M
        USE SUBBOX_M
        implicit none
        
        integer numberOfMove, numberOfAcceptanceMove
        integer iEqui
        
        real*8 AcceptanceRatio
        
!       --------------------------------------------------------------
!       CHECK THE ACCEPTANCE RATIO & THEN ADJUST THE DISPLACEMENT STEP
!       --------------------------------------------------------------
        if(numberofmove.gt.0)then
        AcceptanceRatio    = real(numberOfAcceptanceMove)/real(numberOfMove)
        else
        AcceptanceRatio    = 1.0e0
        endif
        
 !       if(iEqui.lt. Ncycle/2 )then
 !         if(AcceptanceRatio.lt.0.25E0)then
           if(AcceptanceRatio.lt.AccRatio)then
            deltaRx = deltaRx*0.95E0
            deltaRy = deltaRy*0.95E0
            deltaRz = deltaRz*0.95E0
           endif
           if(AcceptanceRatio.gt.AccRatio)then
            deltaRx = deltaRx*1.05E0
            deltaRy = deltaRy*1.05E0
            deltaRz = deltaRz*1.05E0
           endif
!          ----------------------------------------
!          INSURE THAT THE deltaR HAS A UPPER BOUND
!          ----------------------------------------
           if((deltaRx .gt. deltaRxInitial) .OR. (deltaRx .lt. 0.001E0))deltaRx = deltaRxInitial
           if((deltaRy .gt. deltaRyInitial) .OR. (deltaRy .lt. 0.001E0))deltaRy = deltaRyInitial
           if((deltaRz .gt. deltaRzInitial) .OR. (deltaRz .lt. 0.001E0))deltaRz = deltaRzInitial
!       ------------------
!       RESET THE COUNTERS
!       ------------------
        numberOfMove           = 0
        numberOfAcceptanceMove = 0
!        write(*,*)'AcceptanceRatio,deltar',AcceptanceRatio,deltar
        return
        end 
        
        
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       SUBROUTINE EquationOfState
!       CALCULATE DENSITY WITH JOHNSON EQUATION OF STATE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
        subroutine EquationOfState
        USE Constant_M
        USE PHYSICAL_M
        USE MCSETTING_M
        implicit none

        integer i, j, iP

        real*8 rhomin, rhomax, rhoj, rhos(10)
        real*8 P0
        real*8 delta_rho, F1, F2, P1, P2, rho1, rho2
        real*8 rhomolperm3(500)
        
        write(*,*)'The Calculation of the density is in processing'
            
       do iP = 1, NumOfPre
          IF( (Temperature .LT. CriticalTemp) .AND. (P(iP).LT.1.0E0))THEN
                 rho(iP) = Pressure(iP)/Temperature
                 rhoMolPerM3(iP) = rho(iP)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
          ELSEIF(WAGNER .AND. (T.lt.90.6941E0 .OR. P(iP).LT.11696.E0))THEN
                rho(iP) = Pressure(iP)/Temperature
                rhoMolPerM3(iP) = rho(iP)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
          ELSE
             rhomin = 1.0E-15
             rho2 = 0.0E0
             delta_rho = 0.001e0
             j = 0
             do while(rho2.le.1.0e0)
                rho1= rhomin
                IF(JOHNSON)call JohnsonPCalculation( rho1, P1 )
                IF(WAGNER)call WagnerPCalculation( rho1, P1 )
 !               WRITE(*,*)'rho1,P1,Pressure(iP)',rho1,P1,Pressure(iP)
 !               read(*,*)
                F1 = P1 - Pressure(iP)
                rho2= rho1+delta_rho
                IF(JOHNSON)call JohnsonPCalculation( rho2, P2 )
                IF(WAGNER)call WagnerPCalculation( rho2, P2 )
 !               WRITE(*,*)'rho2,P2,Pressure(iP)',rho2,P2,Pressure(iP)
 !               read(*,*)
                F2 = P2 - Pressure(iP)
                if(sign(1.0e0,F1).ne.sign(1.0e0,F2))then
                   rhomin = rho1
                   rhomax = rho2
                   call bisection(Pressure(iP),rhomin,F1,rhomax,F2,rhoj)
                       j = j+1
                       rhos(j) = rhoj
                       rhomin = rho2
                else
                      rhomin = rho2
                endif
              enddo

!            write(*,*)j
!            do i=1, j
!            write(*,*)'i, rhos(i)',i, rhos(i)
!            write(1,*)'i, rhos(i)',i, rhos(i)
!            enddo

            rho(iP) = rhos(1)
            rhoMolPerM3(iP) = rho(iP)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber   
          ENDIF
         write(1,*)'iP, rho(-&mol/m3)',iP, rho(iP), rhoMolPerM3(iP)
      enddo

    
      return
      end




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  SUBROUTINE BISETCION
!  CALCULAT THE DENSITY WITH BISECTION METHOD
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine bisection(Pressure0,rhomin,F1,rhomax,F2,rhoj)
       USE PHYSICAL_M
       USE MCSETTING_M
       implicit none
   
       integer flag
   
       real*8 rhomin, rhomax, rhomid, rhoj
       real*8 F0, F1, F2, P0
       real*8 error
       real*8 Pressure0
      
       flag = 0  
   
       do while(flag .eq. 0)
         error = (rhomax - rhomin)/2.0E0
         rhomid = rhomin + error
         IF(JOHNSON)call JohnsonPCalculation( rhomid, P0 )
         IF(WAGNER)call WagnerPCalculation( rhomid, P0 )
         F0 = P0-Pressure0
            
         if(abs(F0) .lt. 1.0e-10)then
           rhoj  = rhomid 
           flag  = 1
!        write(*,*)'flag, P0', flag, P0
          exit
          return
        endif
      
        if( sign(1.0e0,F1) .eq. sign(1.0e0,F0))then
           rhomin = rhomid
           F1 = F0
        else
           rhomax = rhomid
           F2 = F0
       endif
   
     enddo

     return
     end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE PCALCULATION
!      CALCULAT THE PRESSURE FOR DIFFERENT DENSITY WITH JOHNSON 
!      EQUATION OF STATE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine JohnsonPCalculation( rhom, P0)
      USE PHYSICAL_M
      USE MCSETTING_M
      implicit none

      integer i

      real*8 rhom, P0
      real*8 F, xe(32)
      real*8 sumofA, sumofB
      
      rr = 3.0E0
      F = exp(-rr*(rhom**2))
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     FOR X
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

      xe(1)  = 0.8623085097507421

      xe(2)  = 2.976218765822098

      xe(3)  = -8.402230115796038

      xe(4)  = 0.1054136629203555

      xe(5)  = -0.8564583828174598

      xe(6)  = 1.582759470107601

      xe(7)  = 0.7639421948305453

      xe(8)  = 1.753173414312048

      xe(9)  = 2.798291772190376E+03

      xe(10) = -4.8394220260857657E-02

      xe(11) = 0.9963265197721935

      xe(12) = -3.698000291272493E+01

      xe(13) = 2.084012299434647E+01

      xe(14) = 8.305402124717285E+01

      xe(15) = -9.574799715203068E+02

      xe(16) = -1.477746229234994E+02

      xe(17) = 6.398607852471505E+01

      xe(18) = 1.603993673294834E+01

      xe(19) = 6.805916615864377E+01

      xe(20) = -2.791293578795945E+03

      xe(21) = -6.245128304568454E+00

      xe(22) = -8.116836104958410E+03

      xe(23) = 1.488735559561229E+01

      xe(24) = -1.059346754655084E+04

      xe(25) = -1.131607632802822E+02
     
      xe(26) = -8.867771540418822E+03

      xe(27) = -3.986982844450543E+01

      xe(28) = -4.689270299917261E+03

      xe(29) = 2.593535277438717E+02

      xe(30) = -2.694523589434903E+03

      xe(31) = -7.218487631550215E2
     
      xe(32) = 1.721802063863269E2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     FOR a
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      a(1) = xe(1)*Temperature + xe(2)*sqrt(Temperature) + xe(3) + xe(4)/Temperature + xe(5)/(Temperature**2)
      a(2) = xe(6)*Temperature + xe(7) + xe(8)/Temperature + xe(9)/(Temperature**2)
      a(3) = xe(10)*Temperature + xe(11) + xe(12)/Temperature 
      a(4) = xe(13)
      a(5) = xe(14)/Temperature + xe(15)/(Temperature**2)
      a(6) = xe(16)/Temperature 
      a(7) = xe(17)/Temperature + xe(18)/(Temperature**2)
      a(8) = xe(19)/(Temperature**2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!    FOR b
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      b(1) = xe(20)/(Temperature**2) + xe(21)/(Temperature**3)
      b(2) = xe(22)/(Temperature**2) + xe(23)/(Temperature**4)
      b(3) = xe(24)/(Temperature**2) + xe(25)/(Temperature**3)
      b(4) = xe(26)/(Temperature**2) + xe(27)/(Temperature**4)
      b(5) = xe(28)/(Temperature**2) + xe(29)/(Temperature**3)
      b(6) = xe(30)/(Temperature**2) + xe(31)/(Temperature**3) + xe(32)/(Temperature**4)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!    SUM OF A & B
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
       sumofA = 0.0E0
       sumofB = 0.0E0
       do i = 1, 8
        sumofA = sumofA + a(i)*(rhom**(i+1))
       enddo

       do i = 1, 6
        sumofB = sumofB + b(i)*(rhom**(2*i+1))
       enddo
    
       P0 = rhom*Temperature + sumofA + F*sumofB

      return
      end Subroutine JohnsonPCalculation

 !! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !!
 !!    Subroutine    Pcalculation_Wagner EOS
 !!    FOR METHANE   90.6941K =< T <= 625K
 !!                  0.011696MPa=< P <= 1000MPa
 !!    Reference: U.Setzmann, W.Wagner, 'A New Equation of State and Tables of Thermodynamic
 !!    Properties for Methane Covering the Range from the Melting Line to 625K at Pressures 
 !!    Up to 1000 MPa', J. Phys. Chem. Ref. Data, Vol. 20, NO.6, 1991   
 !!
 !! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine WagnerPcalculation (rhow, Pw)
         USE Constant_M
         USE PHYSICAL_M
         USE FLUID_SOLID_M
         USE MCSETTING_M
         implicit none
         
         integer i
         
         real*8 rhow,Pw,Rw  
         real*8 CriticalTw, CriticalRhow, reducedRho, ireducedT 
         real*8 Wn(50),  Wd(50), Wt(50), Wc(50),WAlpha(50), WBeta(50),WGamma(50), Wdelta(50)
         real*8 residualforP, residualP1, residualP2, residualP3
         
         
         CriticalTw   = 190.564E0 !K
         CriticalRhow = 162.66E0   !kg/m3
         Rw           = 0.5182705E3  ! J/kg/K = Rg/MWofMethane = (8.314 /0.0160428 )=(J/K/mol)/(kg/mol)
         
         CriticalRhow = CriticalRhow/MW*AvogadroNumber/(ScaleLength*1.0e-10)**3.0e0 !!Reduce the unit of the density to dimensionless
         
         reducedRho = rhow/CriticalRhow
  
         ireducedT  = CriticalTw/T   !T is the dimentioanl temperature in (K).
         
         residualP1 = 0.0e0
         residualP2 = 0.0e0
         residualP3 = 0.0e0
         
! ------------------------------------------------     
         Wn(1)  =  0.4367901028E-1
         Wn(2)  =  0.6709236199E0
         Wn(3)  = -0.1765577859E1
         Wn(4)  =  0.8582330241E0
         Wn(5)  = -0.1206513052E1
         Wn(6)  =  0.5120467220E0
         Wn(7)  = -0.4000010791E-3
         Wn(8)  = -0.1247842423E-1
         Wn(9)  =  0.3100269701E-1
         Wn(10) =  0.1754748522E-2
         Wn(11) = -0.3171921605E-5
         Wn(12) = -0.2240346840E-5
         Wn(13) =  0.2947056156E-6
         Wn(14) =  0.1830487909E0
         Wn(15) =  0.1511883679E0
         Wn(16) = -0.4289363877E0
         Wn(17) =  0.6894002446E-1
         Wn(18) = -0.1408313996E-1
         Wn(19) = -0.3063054830E-1
         Wn(20) = -0.2969906708E-1
         Wn(21) = -0.1932040831E-1
         Wn(22) = -0.1105739959E0
         Wn(23) =  0.9952548995E-1
         Wn(24) =  0.8548437825E-2
         Wn(25) = -0.6150555662E-1
         Wn(26) = -0.4291792423E-1
         Wn(27) = -0.1813207290E-1
         Wn(28) =  0.3445904760E-1
         Wn(29) = -0.2385919450E-2
         Wn(30) = -0.1159094939E-1
         Wn(31) =  0.6641693602E-1
         Wn(32) = -0.2371549590E-1
         Wn(33) = -0.3961624905E-1
         Wn(34) = -0.1387292044E-1
         Wn(35) =  0.3389489599E-1
         Wn(36) = -0.2927378753E-2
         Wn(37) =  0.9324799946E-4
         Wn(38) = -0.6287171518E1
         Wn(39) =  0.1271069467E2
         Wn(40) = -0.6423953466E1
 ! ----------------------------------------------------     
         Wd(1)  =  1.0E0
         Wd(2)  =  1.0E0
         Wd(3)  =  1.0E0
         Wd(4)  =  2.0E0
         Wd(5)  =  2.0E0
         Wd(6)  =  2.0E0
         Wd(7)  =  2.0E0
         Wd(8)  =  3.0E0
         Wd(9)  =  4.0E0
         Wd(10) =  4.0E0
         Wd(11) =  8.0E0
         Wd(12) =  9.0E0
         Wd(13) =  10.0E0
         Wd(14) =  1.0E0
         Wd(15) =  1.0E0
         Wd(16) =  1.0E0
         Wd(17) =  2.0E0
         Wd(18) =  4.0E0
         Wd(19) =  5.0E0
         Wd(20) =  6.0E0
         Wd(21) =  1.0E0
         Wd(22) =  2.0E0
         Wd(23) =  3.0E0
         Wd(24) =  4.0E0
         Wd(25) =  4.0E0
         Wd(26) =  3.0E0
         Wd(27) =  5.0E0
         Wd(28) =  5.0E0
         Wd(29) =  8.0E0
         Wd(30) =  2.0E0
         Wd(31) =  3.0E0
         Wd(32) =  4.0E0
         Wd(33) =  4.0E0
         Wd(34) =  4.0E0
         Wd(35) =  5.0E0
         Wd(36) =  6.0E0
         Wd(37) =  2.0E0
         Wd(38) =  0.0E0
         Wd(39) =  0.0E0
         Wd(40) =  0.0E0
 ! ----------------------------------------------
         Wt(1)  =  -0.50E0
         Wt(2)  =   0.50E0
         Wt(3)  =   1.00E0
         Wt(4)  =   0.50E0
         Wt(5)  =   1.00E0
         Wt(6)  =   1.50E0
         Wt(7)  =   4.50E0
         Wt(8)  =   0.00E0
         Wt(9)  =   1.00E0
         Wt(10) =   3.00E0
         Wt(11) =   1.00E0
         Wt(12) =   3.00E0
         Wt(13) =   3.00E0
         Wt(14) =   0.00E0
         Wt(15) =   1.00E0
         Wt(16) =   2.00E0
         Wt(17) =   0.00E0
         Wt(18) =   0.00E0
         Wt(19) =   2.00E0
         Wt(20) =   2.00E0
         Wt(21) =   5.00E0
         Wt(22) =   5.00E0
         Wt(23) =   5.00E0
         Wt(24) =   2.00E0
         Wt(25) =   4.00E0
         Wt(26) =  12.00E0
         Wt(27) =   8.00E0
         Wt(28) =  10.00E0
         Wt(29) =  10.00E0
         Wt(30) =  10.00E0
         Wt(31) =  14.00E0
         Wt(32) =  12.00E0
         Wt(33) =  18.00E0
         Wt(34) =  22.00E0
         Wt(35) =  18.00E0
         Wt(36) =  14.00E0
         Wt(37) =   2.00E0
         Wt(38) =   0.00E0
         Wt(39) =   1.00E0
         Wt(40) =   2.00E0
! ------------------------------------------------
         Wc(14) =  1.0E0
         Wc(15) =  1.0E0
         Wc(16) =  1.0E0
         Wc(17) =  1.0E0
         Wc(18) =  1.0E0
         Wc(19) =  1.0E0
         Wc(20) =  1.0E0
         Wc(21) =  2.0E0
         Wc(22) =  2.0E0
         Wc(23) =  2.0E0
         Wc(24) =  2.0E0
         Wc(25) =  2.0E0
         Wc(26) =  3.0E0
         Wc(27) =  3.0E0
         Wc(28) =  3.0E0
         Wc(29) =  3.0E0
         Wc(30) =  4.0E0
         Wc(31) =  4.0E0
         Wc(32) =  4.0E0
         Wc(33) =  4.0E0
         Wc(34) =  4.0E0
         Wc(35) =  4.0E0
         Wc(36) =  4.0E0
! ---------------------------------------
         WAlpha(37) =   20.0E0
         WAlpha(38) =   40.0E0
         WAlpha(39) =   40.0E0
         WAlpha(40) =   40.0E0
! ----------------------------------------
         WBeta(37) =   200.0E0
         WBeta(38) =   250.0E0
         WBeta(39) =   250.0E0
         WBeta(40) =   250.0E0
! ----------------------------------------
         Wdelta(37) =   1.0E0
         Wdelta(38) =   1.0E0
         Wdelta(39) =   1.0E0
         Wdelta(40) =   1.0E0    
! ----------------------------------------
         Wgamma(37) =   1.07E0
         Wgamma(38) =   1.11E0
         Wgamma(39) =   1.11E0
         Wgamma(40) =   1.11E0    

! --------------------------------------------------------------
!
!     Calculte Pressure
!
! --------------------------------------------------------------        
        do i = 1, 13
        
             residualP1 = residualP1 + Wn(i)*Wd(i)*(reducedRho**(Wd(i)-1.0e0))*(ireducedT**Wt(i))
         
        enddo
        
        do i = 14, 36
            
            residualP2 = residualP2 + exp(-reducedRho**Wc(i))*Wn(i)*(ireducedT**Wt(i))*(reducedRho**(Wd(i)-1.0e0))*(Wd(i)-Wc(i)*(reducedRho**Wc(i)))
            
        enddo
        
        do i = 37, 40
           
           residualP3 = residualP3 + Wn(i)*(reducedRho**Wd(i))*(ireducedT**Wt(i))*exp(-WAlpha(i)*(reducedRho-Wdelta(i))**2.0e0 - WBeta(i)*(ireducedT-WGamma(i))**2.0e0) &
      &                  *(Wd(i)/reducedRho - 2.0e0*WAlpha(i)*(reducedRho-Wdelta(i)))
      
        enddo
        
        residualforP = residualP1 + residualP2 + residualP3
        
        Pw = rhow*Temperature*(1.0e0 + reducedRho*residualforP)
                
!        write(*,*)Pressure, Rho
!        read(*,*)
        
        return
        end Subroutine WagnerPcalculation



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  SUBROUTINE zactivityCALCULATION
!  CALCULAT THE ACTIVITY FOR GCMC WITH JOHNSON EQUATION OF STATE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine zactivityCalculation
      USE PHYSICAL_M
      USE MCSETTING_M
      implicit none
      
      integer i,iP
      
      real*8 Ar
      real*8 F
      real*8 sumofA, sumofB
      
      rr = 3.0e0
      do iP = 1, NumOfPre
         F = exp(-rr*(rho(iP)**2))
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  FOR G
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        G(1) = (1-F)/(2.0E0*rr)
        G(2) = -(F*(rho(iP)**2)-2.0E0*G(1))/(2.0E0*rr)
        G(3) = -(F*(rho(iP)**4)-4.0E0*G(2))/(2.0E0*rr)
        G(4) = -(F*(rho(iP)**6)-6.0E0*G(3))/(2.0E0*rr)
        G(5) = -(F*(rho(iP)**8)-8.0E0*G(4))/(2.0E0*rr)
        G(6) = -(F*(rho(iP)**10)-10.0E0*G(5))/(2.0E0*rr)
     
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  SUM OF A & B
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        sumofA = 0.0E0
        sumofB = 0.0E0
        do i = 1, 8
           sumofA = sumofA + a(i)*(rho(iP)**i)/dble(i)
        enddo

        do i = 1, 6
           sumofB = sumofB + b(i)*G(i)
        enddo

        Ar  = sumofA + sumofB

        ResidualChemicalPotential(iP) = Ar + Pressure(iP)/rho(iP) - Temperature
        zactivity(iP) = rho(iP)*exp(ResidualChemicalPotential(iP)/Temperature)
        
     enddo
     
     return 
     end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ZLRC
!      LONG RANGE CORRECTION FOR ENERGY & VIRIAL PER PARTICLE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine ZLRC(rho11, EnergyLRC, VirialLRC)
       USE Constant_M
       USE PHYSICAL_M
       USE MCSETTING_M
       implicit none
       
       real*8 EnergyLRC, VirialLRC,rho11
                     
!      ----------------------------------
!      LENNARD-JONES 12-6 POTENTIAL MODEL
!      ----------------------------------
       EnergyLRC  = (8.E0*pi/9.E0)*rho11*((1.E0/rCutOff)**9 - 3.D0*(1.E0/rCutOff)**3)
       VirialLRC  = (32.E0*pi/9.E0)*rho11*((1.E0/rCutOff)**9 - 1.5D0*(1.E0/rCutOff)**3)
                   
       return
       end   
                  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE POSITIONCHECK
!      THIS IS TO MAKE SURE THAT THE PARTICLE IS IN THE SIMULATION BOX
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine positionCheck(x, y, z)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       implicit none
              
       real*8 x, y, z
              
!         --------------------
!         CHECK THE X VARIABLE
!         --------------------
          if(pbcx)then
             if (x.gt.CarbonLengthx1) then
                x = x - int(x/CarbonLengthx1)*CarbonLengthx1
             elseif (x.lt.0.0E0) then
                x = x - int(x/CarbonLengthx1)*CarbonLengthx1 + CarbonLengthx1
             endif
          endif
!         --------------------
!         CHECK THE Y VARIABLE
!         --------------------
          if(pbcy)then
            if (y.gt.CarbonLengthy1) then
                y = y - int(y/CarbonLengthy1)*CarbonLengthy1
            elseif (y.lt.0.0E0) then
                y = y - int(y/CarbonLengthy1)*CarbonLengthy1 + CarbonLengthy1
            endif
          endif
!         --------------------
!         CHECK THE Z VARIABLE
!         --------------------
          if(pbcz)then
             if (z.gt.BoxLengthZ) then                     
                 z = z - int(z/BoxLengthZ)*BoxLengthZ
             elseif (z.lt.0.0E0) then
                 z = z - int(z/BoxLengthZ)*BoxLengthZ + BoxLengthZ
             endif
          endif
                
       return 
       endsubroutine positionCheck 
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ACCESSIBLE VOLUME
!      THIS SUBROUTINE IS USED TO GET THE TOTAL ACCESSIBLE VOLUME 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine accessiblevolume
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE POREFIGURE_M
       implicit none
       
       integer i,j,m,iInsert, successinsert,iCycle, iiCycle
       integer bflag, new
       integer Nrotate, Maxrotate
       
       real*8 rdn, AccBulkVolume
       real*8 imaxerror(100000), iAccVolume(100000),error, ierror
       real*8 x1,y1,z1, EnergyC
       real*8 LJxnew(20),LJynew(20),LJznew(20) 
       real*8 CLxnew(20),CLynew(20),CLznew(20)

       iInsert       = 0
       successinsert = 0
       iCycle        = 0
       error         = 1.0e0
       iiCycle       = 0
       new           = Npart + 1 
       IF(DecayESF)THEN
          AccBulkVolume = CollectBulkV
       ELSE
          AccBulkVolume = BulkVolume
       ENDIF
       open(unit=21,file='accparticle',status='unknown' )
!       open(unit=20,file='reasonableinsert',status='unknown' )
!       open(unit=19,file='allinsert',status='unknown' )

       do while (error .gt. 1.0e-3) 
            iCycle = iCycle + 1
            imaxerror(iCycle) = 0.0E0
          do i=1,1000                
             iInsert = iInsert + 1            
             call random_number(rdn)
             x1 = rdn * CarbonLengthx1
             call random_number(rdn)
             IF(DecayESF)THEN
                y1 = CollectLY + rdn*(CollectHY-CollectLY)
             ELSE
                y1 = rdn * CarbonLengthy1
             ENDIF
             call random_number(rdn)
             z1 = rdn * BoxLengthZ 
             if(  inclinePore        .and.  &
             &    z1 .gt. IncCornerZ .and.  &
             &    y1 .lt. ((z1-IncCornerZ)/tan(IncAngle)) )then
                CYCLE  
             endif    
             
             if(  surfpore           .and.  &
             &    z1 .gt. CriLengthZ .and.  &
             &    y1 .lt. CarbonLengthy0 .and. &  
             &    y1 .lt. ((z1-CriLengthZ)/tan(surftheta)) )then
                CYCLE  
             endif    
             
            if(BojanSlit)then      
                 bflag = 0
                 do j = 1, NS
                     if( y1.gt.StripLy(j) .AND.y1.le.StripHy(j)  )then
                         if(   StripAngleY(j) .GT. 0.0 .and. &
             &                z1 .lt. InLzLyStripZ(j) .and. &
             &                (InLzLyStripZ(j)-z1)   .GE. (y1-StripLy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT
                         elseif(   StripAngleY(j) .GT. 0.0 .and. &
             &                 z1 .gt. InHzLyStripZ(j) .and. &
             &                 (z1-InHzLyStripZ(j))   .GE. (y1-StripLy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT  
                         elseif(   StripAngleY(j) .LT. 0.0 .and. &
             &                z1 .lt. InLzLyStripZ(j) .and. &
             &                (InLzLyStripZ(j)-z1)   .GE. (y1-StripHy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT
                         elseif(   StripAngleY(j) .LT. 0.0 .and. &
             &                 z1 .gt. InHzLyStripZ(j) .and. &
             &                 (z1-InHzLyStripZ(j))   .GE. (y1-StripHy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT  
                         elseif(   StripAngleY(j) .EQ. 0.0 .and. &
             &                  ( z1 .lt. InLzLyStripZ(j) .or. z1 .gt. InHzLyStripZ(j) ) ) then     
                               bflag = 1
                               EXIT
                         endif
                   endif
                enddo
                if(bflag .eq. 1)CYCLE
             endif    

!--------------------
!     Rotate
!--------------------
           if(NLJ.GT.1)then
              Nrotate = 1
              Maxrotate = 100
           else
              Nrotate = 1
              Maxrotate = 1
           endif
              DO while (Nrotate .le. Maxrotate)
                CALL ROTATE
                do j = 1, NLJ
                   LJxnew(j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j)
                   LJynew(j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j)
                   LJznew(j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j)
                enddo

                do  j = 1, NCL
                  CLxnew(j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) 
                  CLynew(j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) 
                  CLznew(j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) 
                enddo  
                 
                MCx(new) = x1
                MCy(new) = y1
                MCz(new) = z1
                do m = 1, NLJ
                   LJx(new,m) = LJxnew(m) + MCx(new) - MCx0
                   LJy(new,m) = LJynew(m) + MCy(new) - MCy0
                   LJz(new,m) = LJznew(m) + MCz(new) - MCz0
                enddo
         
                do m =1, NCL
                   CLx(new,m) = CLxnew(m) + MCx(new) - MCx0
                   CLy(new,m) = CLynew(m) + MCy(new) - MCy0
                   CLz(new,m) = CLznew(m) + MCz(new) - MCz0
                enddo               
                call energysingleparticleSF(new, EnergyC )
                IF(EnergyC .le. 0)THEN
                     if(GCMC)then
                          if ((y1.ge. AccVYL) .AND. (y1.le.AccVYH))then 
                             successinsert = successinsert + 1
                             write(21,*)x1*scalelength,y1*scalelength,z1*scalelength
                             EXIT
                          endif  
                     elseif(NVT)then
                          if ((y1.ge. AccVYL) .AND. (y1.le.AccVYH) .AND. (z1.le. AccVZH) )then 
                             successinsert = successinsert + 1
                             EXIT
                          endif  
                     else
                          successinsert = successinsert + 1
                          write(21,*)x1*scalelength,y1*scalelength,z1*scalelength
                          EXIT
                     endif 
                ENDIF
                Nrotate = Nrotate+1  
             ENDDO 
          
           if(mod(i,100).eq.0)then
              iiCycle = iiCycle + 1
              iAccVolume(iiCycle) = dble(successinsert)/dble(iInsert)*AccBulkVolume
              if(iiCycle.gt.1)then
                 ierror = abs(iAccVolume(iiCycle)-iAccVolume(iiCycle-1))/iAccVolume(iiCycle)
              endif
              if(ierror.gt.imaxerror(iCycle))then 
                 imaxerror(iCycle)=ierror 
              endif
          endif
       enddo
             
       if(icycle.gt.1)then
          error = abs(imaxerror(iCycle)- imaxerror(iCycle-1))/imaxerror(iCycle)
       endif
       
       enddo
       
       write(*,*)'iInsert, successinsert',iInsert, successinsert
!       write(21,*)'iInsert, successinsert',iInsert, successinsert
!       write(21,*)'friction',dble(successinsert)/dble(iInsert)
       
       AccVolume = dble(successinsert)/dble(iInsert)*AccBulkVolume
!       write(21,*)'BulkVolume, AccVolume',BulkVolume, AccVolume
     
       return
       end
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE rotate
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          subroutine rotate
          USE Constant_M
          USE MCSETTING_M
          implicit none
        
          real*8 theta,phi,psi
          real*8 q0,q1,q2,q3
          real*8 rdn
          
                            
          call random_number(rdn)
          theta = rdn*2.0e0*pi
          call random_number(rdn)
          phi   = rdn*2.0e0*pi
          call random_number(rdn)
          psi   = rdn*2.0e0*pi

          q0 = cos(theta/2.0e0)*cos((phi+psi)/2.0e0)
          q1 = sin(theta/2.0e0)*cos((phi-psi)/2.0e0)
          q2 = sin(theta/2.0e0)*sin((phi-psi)/2.0e0)
          q3 = cos(theta/2.0e0)*sin((phi+psi)/2.0e0)

          Rot(1) = q0*q0 + q1*q1 -q2*q2 - q3*q3
          Rot(2) = 2.0E0*(q1*q2 + q0*q3)
          Rot(3) = 2.0E0*(q1*q3 - q0*q2)  
          Rot(4) = 2.0E0*(q1*q2 - q0*q3)
          Rot(5) = q0*q0 - q1*q1 + q2*q2 - q3*q3 
          Rot(6) = 2.0E0*(q2*q3 + q0*q1)
          Rot(7) = 2.0E0*(q1*q3 + q0*q2) 
          Rot(8) = 2.0E0*(q2*q3 - q0*q1)
          Rot(9) = q0*q0 - q1*q1 - q2*q2 + q3*q3 

          return
          end
            
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE STEELE - POTENTIAL
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine steelepotential(i, Energy)
       USE Constant_M
       USE FLUID_SOLID_M
       USE MCSETTING_M
       implicit none
    
       integer i,j, k 
       
       real*8 Energy, iEnergy, rhosperM3, rho_s
       real*8 mcdistanceZ,distanceZ
       real*8 SigmaSF, WellDepthSF
       real*8 decayf
       
              
       rhosperM3     = 114.0E27
       
       Energy = 0.0e0
       
       rho_s         = rhosperM3*((SCALELENGTH*1.0E-10)**3)
!       if(mcdistanceZ.gt.rcutoff)then
!         Energy = 0.0E0
!         return
!       endif
       
!$OMP  parallel do default(none)                                         &
!$OMP   private (distancez,ienergy,j,k,sigmasf,welldepthsf)             &
!$OMP   shared  (i,graphenelayer,ljz,nlj,nsteele,pi,psteele,rho_s,      &
!$OMP            sigmaff,sigmass,welldepthff,welldepthss)               &
!$OMP   reduction (+:energy)
       DO k = 1, Nsteele
         do j=1,NLJ
            sigmaSF     = (sigmaFF(j) + sigmaSS)/2.0e0
            welldepthSF = sqrt(welldepthFF(j)*welldepthSS) 
            distanceZ = abs(LJz(i,j)-Psteele(k))
!            if(distanceZ.gt.rcutoff)then
!              iEnergy = 0.0e0
!            else
            iEnergy = 4.0E0*pi*WellDepthSF*rho_s*(SigmaSF**2)*GrapheneLayer*((SigmaSF/distanceZ)**10/5.0E0 &
&               - (SigmaSF/distanceZ)**4/2.0E0 - SigmaSF**4/(6.0E0*GrapheneLayer*(distanceZ+0.61e0*GrapheneLayer)**3))
!            endif
          Energy = Energy + iEnergy
         enddo
       ENDDO
       IF(DecayESF)THEN
          call decayfcalculation(i,decayf)
          Energy = Energy*decayf
       ENDIF
       return
       end 
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE RADIUS DISTRIBUTION
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
           subroutine radiusdistribution(switch)
           USE Constant_M
           USE FLUID_SOLID_M
           USE PHYSICAL_M
           USE MCSETTING_M
           USE FLUCTUATION_M
           implicit none
           
           integer bin, i, j, k, ngr, RmaxBin
           integer switch,l,m,n
           integer FLAG1, FLAG2, SFLAG1, SFLAG2
           
           real*8 radiusB,  deltaBin
           real*8 r(1000),rA(1000)
           real*8 vb, nid
           real*8 radN, NpartML, NpartSL
           real*8 mcBin(1000),NNBin(1000),mcNBin(1000)
           real*8 mcSum_Nbin(1000),NNSum_NBin(1000),mcNSum_NBin(1000)
           real*8 Av_mcBin(1000),Av_NNBin(1000),Av_mcNBin(1000)
           real*8 mcxij,mcyij,mczij
           real*8 rmcx,rmcy,rmcz,mcradiusD
           real*8 rNNx,rNNy,rNNz,NNradiusD
           real*8 rmcNx,rmcNy,rmcNz,mcNradiusD
           real*8 mcBin_ML(1000),TotmcBin_ML(1000),Av_mcBin_ML(1000) !! ML = 'Mono Layer', for coordination number in the first layer
           real*8 mcBin_SL(1000),TotmcBin_SL(1000),Av_mcBin_SL(1000) !! SL = 'Second Layer' for coordination number in the 2nd layer
           
           
           common/gr/ngr
           common/FORRD1/mcBin,NNBin,mcNBin
           common/FORRD2/mcSum_Nbin,NNSum_NBin,mcNSum_NBin
           common/FORRD3/Av_mcBin,Av_NNBin,Av_mcNBin
           common/FORRD3/mcBin_ML,TotmcBin_ML,Av_mcBin_ML
           common/FORRD4/mcBin_SL,TotmcBin_SL,Av_mcBin_SL
           
           deltaBin   = 0.1e0/scalelength
           radiusB    = BoxLength/2.0e0  
           RmaxBin     = int(radiusB/deltaBin) 
           MLBin      = int(6.0/scalelength/deltaBin)+1
           
           IF(switch.eq.0)THEN
               open(unit=29, file='radius.txt',status='unknown')
           ENDIF

           IF(switch.eq.1)THEN
               mcSum_NBin   = 0.0e0
               NNSum_NBin   = 0.0e0
               mcNsum_NBin  = 0.0e0
               TotmcBin_ML  = 0.0E0
               TotmcBin_SL  = 0.0E0
               Av_mcBin     = 0.0e0 
               Av_NNBin     = 0.0e0
               Av_mcNBin    = 0.0e0
               Av_mcBin_ML  = 0.0E0
               Av_mcBin_SL  = 0.0E0
               ngr        = 0
               write(29,'(9A16)')'Pressure(pa)','Distance(A)', 'g(r)-mc', 'g(r)-NN', 'g(r)-mcN', 'CoordN_1st','AccCoorN_1st', 'CoordN_2nd','AccCoorN_2nd'
               write(29,*)'================================================================================================'
           ENDIF
         
         IF(switch.eq.2)THEN
              ngr      = ngr + 1
              mcBin    = 0.0e0
              NNBin    = 0.0e0
              mcNBin   = 0.0e0
              radN     = 0.0E0
              mcBin_ML = 0.0E0
              mcBin_SL = 0.0E0
              NpartML  = 0.0E0
              NpartSL  = 0.0E0

              if(Npart .gt. 1)then
                 do i = 1, Npart-1 
                     FLAG1 = 0
                     SFLAG1 = 0
                     IF(mcz(i).LE.SMLimitZ)THEN
                       FLAG1=1
                       NpartML = NpartML + 1
                     ENDIF
                     IF( mcz(i).GT.SMLimitZ .AND. mcz(i).LE.SMLimitZ2)THEN
                       SFLAG1=1
                       NpartSL = NpartSL + 1
                     ENDIF
 !                   if(LineSurf)then
 !                      if(mcy(i).lt. (CarbonLengthy1/2.0e0).OR. &
 !                      &  mcy(i).gt. (CarbonLengthy1/2.0e0 + 1.0e0))cycle
 !                   endif
                    radN = radN + 1.0E0  
                   do j = i+1, Npart
                      FLAG2 = 0
                      SFLAG2 = 0
                      IF(mcz(j).LE.SMLimitZ)FLAG2=1
                      IF( mcz(j).GT.SMLimitZ .AND. mcz(j).LE.SMLimitZ2)SFLAG2=1
 !                    if(LineSurf)then
 !                      if(mcy(i).lt. (CarbonLengthy1/2.0e0).OR. &
 !                      &  mcy(i).gt. (CarbonLengthy1/2.0e0 + 1.0e0))cycle
 !                    endif
                      mcxij = abs(mcx(i) - mcx(j))
                      mcyij = abs(mcy(i) - mcy(j))
                      mczij = abs(mcz(i) - mcz(j))
!------------------------------------------------------
!  Radius distribution for Mass Center
!------------------------------------------------------    
                      rmcx  = mcxij
                      rmcy  = mcyij
                      rmcz  = mczij                            
                      if(pbcx)then
                         if(mcxij .gt. CarbonLengthx1/2.0e0) rmcx = rmcx - CarbonLengthx1
                      endif
                      if(pbcy)then
                        if(mcyij .gt. CarbonLengthy1/2.0e0) rmcy = rmcy - CarbonLengthy1
                      endif
                      if(pbcz)then
                        if(mczij .gt. BoxLengthZ/2.0e0) rmcz = rmcz - BoxLengthZ 
                      endif              
                      mcradiusD = sqrt(rmcx*rmcx + rmcy*rmcy + rmcz*rmcz)
                      if(mcradiusD .le. radiusB)then
                         bin = int(mcradiusD/deltaBin)+1
                         mcBin(bin) = mcBin(bin) + 1.0e0
                         IF(FLAG1.EQ.1 .AND. FLAG2.EQ.1)mcBin_ML(bin) = mcBin_ML(bin) + 1.0          !ML=MonoLayer For the Coordination number in the first layer
                         IF(SFLAG1.EQ.1 .AND. SFLAG2.EQ.1)mcBin_SL(bin) = mcBin_SL(bin) + 1.0 
                      endif
!------------------------------------------------------
!  Radius distribution for N - N
!------------------------------------------------------    
                     do m = 1, NLJ                
                        do n = 1, NLJ                  ! the N site in the different molecular j        
                          rNNx = abs(LJx(i,m)- LJx(j,n))
                          rNNy = abs(LJy(i,m)- LJy(j,n))
                          rNNz = abs(LJz(i,m)- LJz(j,n))
                          if(pbcx)then
                            if(mcxij .gt. CarbonLengthx1/2.0e0) rNNx = rNNx - CarbonLengthx1
                          endif
                          if(pbcy)then
                            if(mcyij .gt. CarbonLengthy1/2.0e0) rNNy = rNNy - CarbonLengthy1
                          endif
                          if(pbcz)then
                            if(mczij .gt. BoxLengthZ/2.0e0) rNNz = rNNz - BoxLengthZ 
                          endif              
                          NNradiusD = sqrt(rNNx*rNNx + rNNy*rNNy + rNNz*rNNz)
                          if(NNradiusD .le. radiusB)then
                             bin = int(NNradiusD/deltaBin)+1
                             NNBin(bin) = NNBin(bin) + 1.0e0
                          endif
                       enddo 
                     enddo
!------------------------------------------------------
!  Radius distribution for mass center - N
!------------------------------------------------------                      
                      do n = 1, NLJ                  ! the N site in the different molecular j        
                         rmcNx = abs(LJx(i,m)- LJx(j,n))
                         rmcNy = abs(LJy(i,m)- LJy(j,n))
                         rmcNz = abs(LJz(i,m)- LJz(j,n))
                         if(pbcx)then
                           if(mcxij .gt. CarbonLengthx1/2.0e0) rmcNx = rmcNx - CarbonLengthx1
                         endif
                         if(pbcy)then
                           if(mcyij .gt. CarbonLengthy1/2.0e0) rmcNy = rmcNy - CarbonLengthy1
                         endif
                         if(pbcz)then
                           if(mczij .gt. BoxLengthZ/2.0e0) rmcNz = rmcNz - BoxLengthZ 
                         endif              
                         mcNradiusD = sqrt(rmcNx*rmcNx + rmcNy*rmcNy + rmcNz*rmcNz)
                         if(mcNradiusD .le. radiusB)then
                            bin = int(mcNradiusD/deltaBin)+1
                            mcNBin(bin) =mcNBin(bin) + 1.0e0
                         endif
                      enddo  
                   
                    enddo
                  enddo
                  
                     IF(mcz(Npart).LE.SMLimitZ)NpartML = NpartML + 1
                     IF( mcz(Npart).GT.SMLimitZ .AND. mcz(Npart).LE.SMLimitZ2)NpartSL = NpartSL + 1
                  do k = 1,Rmaxbin
!                     if(LineSurf)then 
!                       if(radN.gt. 0.0)then
!                       mcsum_NBin(k) = mcsum_NBin(k) + 2.0e0*mcBin(k)/radN
!                       endif
!                     else
                       mcsum_NBin(k) = mcsum_NBin(k) + 2.0e0*mcBin(k)/dble(Npart)
                       NNsum_NBin(k) = NNsum_NBin(k) + 2.0e0*NNBin(k)/dble(NLJ*Npart)
                       mcNsum_NBin(k) = mcNsum_NBin(k) + 2.0e0*mcNBin(k)/dble(Npart)
                       IF(NpartML.GT.0)TotmcBin_ML(k) = TotmcBin_ML(k) + 2.0*mcBin_ML(k)/dble(NpartML)
                       IF(NpartSL.GT.0)TotmcBin_SL(k) = TotmcBin_SL(k) + 2.0*mcBin_SL(k)/dble(NpartSL)
!                     endif
                  enddo              
                endif
         ENDIF
         
         IF(switch.eq.3)THEN 
           Accu_mcBin_ML=0.0
           Accu_mcBin_SL=0.0
           do k = 1, Rmaxbin
              r(k)   = deltaBin*(k-0.5)
              rA(k)  = r(k)*scalelength
              vb     = (k**3 - (k-1.0e0)**3)*deltaBin**3
              IF(GCMC)nid    = (4.0e0/3.0e0)*pi*vb*rho(iPC)
              IF(NVT)nid    = (4.0e0/3.0e0)*pi*vb*rho(iden)
              Av_mcbin(k) =  mcSum_Nbin(k)/(ngr*nid)
              Av_NNbin(k) =  NNSum_Nbin(k)/(ngr*nid)
              Av_mcNbin(k) =  mcNSum_Nbin(k)/(ngr*nid)
              Av_mcBin_ML(k) = TotmcBin_ML(k)/ngr
              Av_mcBin_SL(k) = TotmcBin_SL(k)/ngr
              IF(k.eq.1)THEN
                 Accu_mcBin_ML(k) = Accu_mcBin_ML(1)
                 Accu_mcBin_SL(k) = Accu_mcBin_SL(1)
              ELSE
                 Accu_mcBin_ML(k) = Accu_mcBin_ML(k-1) + Av_mcBin_ML(k)
                 Accu_mcBin_SL(k) = Accu_mcBin_SL(k-1) + Av_mcBin_SL(k)
              ENDIF
              IF(GCMC)write(29,'(9e16.7)')P(iPC),rA(k),Av_mcbin(k),Av_NNbin(k),Av_mcNbin(k),Av_mcBin_ML(k),Accu_mcBin_ML(k),Av_mcBin_SL(k),Accu_mcBin_SL(k)
              IF(NVT)write(29,'(I8,8e16.7)')Npart,rA(k),Av_mcbin(k),Av_NNbin(k),Av_mcNbin(k),Av_mcBin_ML(k),Accu_mcBin_ML(k),Av_mcBin_SL(k),Accu_mcBin_SL(k)
           enddo
         ENDIF
         
         RETURN
         ENDSUBROUTINE radiusdistribution
           
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE DISTANCE DISTRIBUTION
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
             
        
     
        
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   SUBROUTINE COULOMB FORCE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
        subroutine coulombforce(cld, charge1,charge2,CLU,CLV)
        USE Constant_M
        USE MCSETTING_M
        implicit none
        
        real*8 cld,charge1,charge2, CLU,CLV
        
                
!        if((cld.lt.0.4e0) .AND. (charge1*charge2.lt.0.0e0))then
!           CLU = 1.0e10
!           return
!        else       
        CLU = charge1*charge2/(4.0e0*pi*permittivity*Scalelength*1.0e-10*scaleEnergy*kB*cld)
        CLV = CLU/3.0
!        endif
        return
        end
     
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Bojan Potential
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     subroutine BojanPotential(bi,i,BU) 
     USE Constant_M
     USE FLUID_SOLID_M
     USE MCSETTING_M
     USE POREFIGURE_M
     implicit none
     
     integer bi, i, j, k
     
     real*8 y,z,BU,Cz(20),CY
     real*8 totforce
     real*8 BWellDepthSF
     real*8 layerTotal
     real*8 Inctotforce, DisP_L, DisP_CZ, DisCY
     real*8 Clo1totforce,Clo2totforce
     real*8 Line_A1(20), Line_B1(20), Line_C1(20), VLine_A2(20), VLine_B2(20), VLine_C2(20)
     real*8 Project_y

     
!       rhosperM2  = 38.2E18       ! m-2

!       rho_s      = rhosperM2*((SCALELENGTH*1.0E-10)**2)  
       BU = 0.0E0    
       Cz = 0.0e0
       
       IF(bi.le.NS)THEN    ! the z coordinate of different layers  bi: strip number
!          rho_s = rhosperM2(bi)*((SCALELENGTH*1.0E-10)**2) 
          do j = 1, StriplayerN(bi)
             layerTotal = StriplayerN(bi)
             
             if(TopStrip(bi))then
                 Cz(j) = StripZ0(bi)+dble(j-1)*Stripgap(bi)/cos_StripAngleY(bi)
             else
                 Cz(j) = StripZ0(bi)-dble(j-1)*Stripgap(bi)/cos_StripAngleY(bi)
             endif
             IF(StripAngleY(bi).GT.0.0)THEN
                Line_A1(j) = tan_StripAngleY(bi)
                Line_B1(j) = 1.0E0
                Line_C1(j) = -Cz(j)-StripLy(bi)*tan_StripAngleY(bi)
             ELSEIF(StripAngleY(bi).LT.0.0)THEN
                 Line_A1(j) = -tan_StripAngleY(bi)
                 Line_B1(j) = -1.0E0
                 Line_C1(j) = Cz(j)+StripLy(bi)*tan_StripAngleY(bi)
             ENDIF    
             
             
             if(BojanSlit)then
                layerTotal = StriplayerN(bi)*2
                if(StripZ0(bi).gt. 0.0e0)then
                   Cz(j+int(StriplayerN(bi))) = Cz(j)+BoxlengthZ-2.0e0*StripZ0(bi)+dble(StriplayerN(bi)-1)*Stripgap(bi)/cos_StripAngleY(bi)
                else
                   Cz(j+int(StriplayerN(bi))) = Cz(j)+BoxlengthZ+dble(StriplayerN(bi)-1)*Stripgap(bi)/cos_StripAngleY(bi)
                endif
                IF(StripAngleY(bi).GT.0.0)THEN
                    Line_A1(j+int(StriplayerN(bi))) = tan_StripAngleY(bi)
                    Line_B1(j+int(StriplayerN(bi))) = -1.0
                    Line_C1(j+int(StriplayerN(bi))) = Cz(j+int(StriplayerN(bi)))-StripLy(bi)*tan_StripAngleY(bi)
                ELSEIF(StripAngleY(bi).LT.0.0)THEN
                    Line_A1(j+int(StriplayerN(bi))) = -tan_StripAngleY(bi)
                    Line_B1(j+int(StriplayerN(bi))) = 1.0
                    Line_C1(j+int(StriplayerN(bi))) = -Cz(j+int(StriplayerN(bi)))+StripLy(bi)*tan_StripAngleY(bi)
                ENDIF    
             endif
           enddo
        ENDIF

! Lines after !!! contains assignment to a module variable.
!$OMP critical
       do k = 1, NLJ
           IF(bi.le.NS)THEN
              totforce   = 0.0e0             
!!!
              BsigmaSF = (sigmaFF(k)+BsigmaSS(bi))/2.0e0
              BwelldepthSF = sqrt(welldepthFF(k)*BwelldepthSS(bi))
              IF(StripAngleY(bi).eq.0.0)THEN
                  Cy =  LJy(i,k) - StripC(bi)          
!!!
                  Py =  StripW(bi)/2.0e0 - Cy    ! the y coordinates of the right-hand edge and left-hand edge relative to the fluid particle
!!!
                  Ny = -StripW(bi)/2.0e0 - Cy    ! Py: positive y, Ny: negative y, give the location of the edges of the strip relative to the position of the adsorbate
          
                  do j = 1, layerTotal
!!!
                     BdeltaZ = abs(Cz(j)-LJz(i,k)) ! the distance between a LJ site and one of the strips in the z direction
!!!
                     call BojanCalculation          
!!!
                     totforce = totforce + subforce
                  enddo
              ELSE
                  do j = 1, layerTotal
                     VLine_A2(j) = Line_B1(j)/Line_A1(j)
                     VLine_B2(j) = -1.0E0
                     VLine_C2(j) = LJz(i,k)-VLine_A2(j)*LJy(i,k) 
!!!
                     BdeltaZ = abs(Line_A1(j)*LJy(i,k) + Line_B1(j)*LJz(i,k)+Line_C1(j))/sqrt(Line_A1(j)**2.0 + Line_B1(j)**2.0)
                     Project_y = ( Line_B1(j)*VLine_C2(j) - VLine_B2(j)*Line_C1(j) )/ (Line_A1(j)*VLine_B2(j) - VLine_A2(j)*Line_B1(j))
!!!
                     Cy = (  Project_y - StripC(bi) )/cos_StripAngleY(bi)
!!!
                     Py = (  StripW(bi)/2.0e0 )/cos_StripAngleY(bi) - Cy
!!!
                     Ny = ( -StripW(bi)/2.0e0 )/cos_StripAngleY(bi) - Cy
!!!
                     call BojanCalculation          
!!!
                     totforce = totforce + subforce
                  enddo
              ENDIF
       
              BU = BU + 2.0e0*pi*rhosperM2(bi)*(BsigmaSF**2.0e0)*BwelldepthSF*totforce
            ENDIF
 
          IF(bi.eq.NS+1)then   ! For InclinPore
             Inctotforce = 0.0e0
!!!
             BsigmaSF = (sigmaFF(k)+IncsigmaSS)/2.0e0
             BwelldepthSF = sqrt(welldepthFF(k)*IncwelldepthSS)
             DisP_L = abs(Inc_A*LJy(i,k)+Inc_B*LJZ(i,k)+Inc_C)/sqrt(Inc_A**2.0+Inc_B**2)
             DisP_CZ = sqrt(LJy(i,k)**2.0 + (LJZ(i,k)-IncCornerZ)**2.0)
             DisCy = sqrt(DisP_CZ**2.0-DisP_L**2.0)
!!!
             Cy = DisCy - IncMid_y
!!!
             Py =  IncLengthY/2.0e0 - Cy
!!!
             Ny = -IncLengthY/2.0e0 - Cy
             do j = 1, IncLayerN
!!!
                BdeltaZ = DisP_L + (j-1)*Incgap
!!!
                call BojanCalculation
                Inctotforce = Inctotforce + subforce
             enddo
             BU = BU + 2.0e0*pi*IncrhosperM2*(BsigmaSF**2.0e0)*BwelldepthSF*Inctotforce
         ENDIF
         
         IF(bi.eq.NS+2)then   ! For Close1End .OR. Close2Ends
             IF(Close1End)THEN
                 Clo1totforce = 0.0e0
                 IF(Clo1Hard)Then
                    Clo1totforce = 0.0e0
                 ELSE
!!!
                    BsigmaSF = (sigmaFF(k)+Clo1sigmaSS)/2.0e0
                    BwelldepthSF = sqrt(welldepthFF(k)*Clo1welldepthSS)   
!!!
                    Cy =  LJZ(i,k)- Clo1Mid_z
!!!
                    Py =  Clo1LengthZ/2.0e0 - Cy
!!!
                    Ny = -Clo1LengthZ/2.0e0 - Cy
                    do j = 1, Clo1LayerN
!!!
                       BdeltaZ = LJy(i,k)+(j-1)*Clo1gap + Clo1gap ! Add the gap between the pore walls and the closed end
!!!
                       call BojanCalculation
                       Clo1totforce = Clo1totforce + subforce
                    enddo
                 ENDIF
                 BU = BU+2.0e0*pi*Clo1rhosperM2*(BsigmaSF**2.0e0)*BwelldepthSF*Clo1totforce  
             ENDIF
             IF(Close2Ends)THEN
                 Clo2totforce = 0.0e0
!!!
                 BsigmaSF = (sigmaFF(k)+Clo2sigmaSS)/2.0e0
                 BwelldepthSF = sqrt(welldepthFF(k)*Clo2welldepthSS)   
                 Cy =  LJZ(i,k)- Clo2Mid_z
!!!
                 Py =  Clo2LengthZ/2.0e0 - Cy
!!!
                 Ny = -Clo2LengthZ/2.0e0 - Cy
                 do j = 1, Clo2LayerN
!!!
                    BdeltaZ = CarbonLengthy1-LJy(i,k)+(j-1)*Clo2gap + Clo2gap ! Add the gap between the pore walls and the closed end
!!!
                    call BojanCalculation
!!!
                    Clo2totforce = Clo2totforce + subforce
                 enddo
                 BU = BU + 2.0e0*pi*Clo2rhosperM2*(BsigmaSF**2.0e0)*BwelldepthSF*Clo2totforce
             ENDIF
         ENDIF
          
       enddo
!$OMP end critical
    return
    end         
!+++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Bojan Calculation
!
!+++++++++++++++++++++++++++++++++++++              
                subroutine BojanCalculation
                USE FLUID_SOLID_M
                implicit none
                
                real*8 Prepulsion,Pattration,Nrepulsion,Nattration
                real*8 subforce1
                real*8 sub1,sub2,sub3,sub4,sub5,sigmadNy,sigmadPy, zdsigma2
     
                if(  (BdeltaZ .lt. 0.1e0) &                ! (z .lt. StripZ0(bi)) .and.(deltaZ .lt. 0.1e0) &        
               !                                               &   .or. (z .gt. StripZ0(bi) .and. (z-StripZ0(bi)).lt.0.1e0) &
               ! &   .and. (Py .ne. 0.0e0) .and. (Ny .ne. 0.0e0) 
                &   .and. (Py*Ny .gt. 0.0e0)  )then
            
!                  sigmadNy = BsigmaSF(bi)/Ny
!                  sigmadPy = BsigmaSF(bi)/Py
                  sigmadNy = BsigmaSF/Ny
                  sigmadPy = BsigmaSF/Py
                  zdsigma2 = (BdeltaZ/BsigmaSF)**2.0e0
                  sub1 = 6.3e1/1.280e3*(sigmadNy**10.0e0 - sigmadPy**10.0e0) - 3.0e0/1.6e1*(sigmadNy**4.0e0 - sigmadPy**4.0e0) 
                  sub2 = - zdsigma2 * ( 2.31e2/1.024e3*(sigmadNy**12.0e0 - sigmadPy**12.0e0) - 5.0e0/1.6e1*(sigmadNy**6.0e0 - sigmadPy**6.0e0) )
                  sub3 = zdsigma2**2.0e0*( 1.287e3/2.048e3*(sigmadNy**14.0e0 - sigmadPy**14.0e0) - 1.05e2/2.56e2*(sigmadNy**8.0e0 - sigmadPy**8.0e0) )
                  sub4 = -zdsigma2**3.0e0*( 4.5045e4/3.2768e4*(sigmadNy**16.0e0 - sigmadPy**16.0e0) - 6.3e1/1.28e2*(sigmadNy**10.0e0 - sigmadPy**10.0e0) )
                  sub5 = zdsigma2**4.0e0*( 8.5085e4/3.2768e4*(sigmadNy**18.0e0 - sigmadPy**18.0e0) - 1.155e3/2.048e3*(sigmadNy**12.0e0 - sigmadPy**12.0e0) )
               
                  subforce1 = sub1 + sub2 + sub3 + sub4 + sub5 
                  if((Py .gt. 0.0e0) .and. (Ny .gt.0.0e0)) subforce = subforce1
                  if((Py .lt. 0.0e0) .and. (Ny .lt.0.0e0)) subforce = - subforce1
               else
            
                  call Brepatt(Py,BdeltaZ,BsigmaSF,Prepulsion,Pattration)
                  call Brepatt(Ny,BdeltaZ,BsigmaSF,Nrepulsion,Nattration)
                  subforce  = (Prepulsion - Nrepulsion)- (Pattration - Nattration)
               endif
               
               return
               end
    
!+++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Bojan Repulsion and attraction
!
!+++++++++++++++++++++++++++++++++++++++++++++++++             
       pure subroutine Brepatt (y, z, BsigmaSF, Repulsion, Attraction)
       implicit none
       
       intent(in) :: z,y,BsigmaSF
       intent(out) :: Repulsion, Attraction
       real*8 y,z,z2,sumy2z2
       real*8 Repulsion, Attraction
       real*8 BsigmaSF
       
       
            z2 = z*z
            sumy2z2 = y*y + z2
            
            Repulsion = y*(BsigmaSF**10.0e0)* ( 0.2e0/(z2**4.0e0) + 0.1e0/(z2**3.0e0)/sumy2z2    &
 &                      + 3.0e0/4.0e1/((z2*sumy2z2)**2.0e0) + 1.0e0/1.6e1/z2/(sumy2z2**3.0e0)  &
 &                      + 7.0e0/1.28e2/(sumy2z2**4.0e0) )/Z2/sqrt(sumy2z2)
            
            Attraction = y*(BsigmaSF**4.0e0)*(0.5e0/z2 + 0.25e0/sumy2z2)/Z2/sqrt(sumy2z2)
            
            
       return     
       end      
            
            
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE BULK_GCMC
!    ARRANGE THE PARTICLES IN THE CUBIC LATTICE CONFIGURATION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
    subroutine BULK_GCMC
    USE Constant_M
    USE FLUID_SOLID_M
    USE PHYSICAL_M
    USE MCSETTING_M
    implicit none
    
    integer iP
    integer numofinsert, successinsert, numofdelete, successdelete
    integer numberOfMove, numberOfAcceptanceMove
    integer iEqui, iMove, iCycle
    integer i, j, cutoffopt, rcutoffset
        
    real*8 rdn
    real*8 secnds, time
    real*8 numberofequilibrium, EnergyTotal, totAcceptanceRatio
    real*8 Energy, CLEnergy,EnergyC, TEnergy
    real*8 Sum_TEnergy, Sum_Energy, Sum_EnergyC
    real*8 Sum_Npart, Sum_TEN, Sum_EN, Sum_ECN, Sum_NNpart                 
    real*8 numberOfSampling
    real*8 Av_Npart, Ad_Ninpore, Ns_rho, Nb_rho
    real*8 Iso_Energy, Iso_EnergyC, Iso_TEnergy, Iso_Npart
    real*8 Iso_EN, Iso_ECN, Iso_TEN,Iso_NN
    real*8 Iso_HeatFF, Iso_HeatSF,Iso_Heat 
    real*8 Ns_rhoMolPerM2, Nb_rhoMolPerM3, rhoMolPerM3
    real*8 EnergyAverageJoulePerMol
    real*8 Iso_HeatJoulePerMol, Iso_HeatFFJoulePerMol, Iso_HeatSFJoulePerMol                                                                                                  
    real*8 AcceptanceRatio, EnergyAverage 
       
          
      pbcx = .TRUE.
      pbcy = .TRUE.
      pbcz = .TRUE.
      Ncycle = 10000
      Ncycles = 10000
      Nmove = 1000

      Steele = .FALSE.
      Bojan  = .FALSE.
      NC3 = 0      

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    CALCULATION STARTS HERE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!    ------------------------------------------------------  
!    Set rCutOff which can't larger than BoxLength/2
!    ------------------------------------------------------
     if(CarbonLengthx1 .gt. CarbonLengthy1)then
        BoxLength = CarbonLengthy1
        else
        BoxLength = CarbonLengthx1
     endif
     if(BoxLength .gt. BoxLengthZ) BoxLength = BoxLengthZ
     rcutoff    = BoxLength/2.0E0
     rcutoffset = 5.0E0
          
     if(rcutoff.gt.rcutoffset)then
       write(*,*)'------------------------------------- '
       write(*,*)'Which one would choose as the rcutoff:'
       write(*,*)'------------------------------------- '
       write(*,*)'1-', rcutoff
       write(*,*)'2-', rcutoffset
!       read(*,*) cutoffopt
       if(NLJ.GT.1)cutoffopt=1
       cutoffopt = 2
       if(cutoffopt.eq.2)rcutoff = rcutoffset
     endif
     
        
     write(1,*)'rCutOff',rCutOff
     write(*,*)'rCutOff',rCutOff
     r2CutOff = rCutOff*rCutOff
     
      write(*,*)'---------------------------------------------- '
      write(*,*)'ENTER THE MAXIMUM LJ-DISPLACEMENT STEP (0.5) : '
      write(*,*)'---------------------------------------------- '     
      deltaRxInitial = Carbonlengthx1/2.0
      deltaRyInitial = Carbonlengthy1/2.0
      deltaRzInitial = BoxLengthZ/2.0

!    --------------------------------
!    SURFACE AREA AND VOLUME VALUES
!    --------------------------------
     TotSurfaceArea  = CarbonLengthx1*CarbonLengthy1              ! for Surface
     BulkVolume      = CarbonLengthx1*CarbonLengthy1*BoxLengthZ   ! the volume of the simulation box 
     BBulkVolume     = BulkVolume
     
     write(*,*) 'BULKPHASE GCMC'
     write(1,*) 'BULKPHASE GCMC'
     write(*,1130) TotSurfaceArea, BulkVolume
     write(1,1130) TotSurfaceArea, BulkVolume
1130 format( 1x, 'SurfaceArea(-)--------------------------:', f20.10,/, &
  &          1x, 'BulkVolume(-)---------------------------:', f20.10,/)  
     
     Npart = 0 
    do iPC = 1 , NumOfPre
       write(*,*)'---------------------------------------------------------------------------'
       write(*,*)'BULKPHASE The pressure being calculated is _ of _ - _:', iPC,NumOfPre,P(iPC)
       write(*,*)'---------------------------------------------------------------------------'    
!      -------------------------
!      SET THE DISPLACEMENT STEP 
!      -------------------------
       deltaRx = deltaRxInitial
       deltaRy = deltaRyInitial
       deltaRz = deltaRzInitial
       write(*,*)'deltaRx,deltaRy,deltaRz',deltaRx,deltaRy,deltaRz
!    --------------------------------
!    INITIAL POSITION ALL N PARTICLES
!    --------------------------------

!     call Initialization   
!    -------------------------------- 
!    EQUILIBRIUM STEPS FOR BULKPHASE
!    --------------------------------
     numofinsert            = 0
     successinsert          = 0
     numofdelete            = 0
     successdelete          = 0
     numberOfMove           = 0
     numberOfAcceptanceMove = 0
     numberofequilibrium    = 0.0
     sum_npart              = 0.0
     EnergyTotal            = 0.0
     totAcceptanceRatio     = 0.0
     call ConfigEnergy(Energy, CLEnergy,EnergyC)
          TEnergy = Energy + CLEnergy 
          write(*,*)'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy
          write(1,*)'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy
    
     do iEqui = 1, Ncycle
         if (mod(iEqui,1000).eq.0)then
            write(*,*)'================BULK EQULIBRATION==================='
            write(1,*)'================BULK EQULIBRATION==================='
            write(*,645)iEqui, sum_Npart/numberofequilibrium, EnergyTotal/numberofequilibrium
            write(1,645)iEqui, sum_Npart/numberofequilibrium, EnergyTotal/numberofequilibrium
      645   format(1x,'CYCLES:',I8,'  Npart:',f15.6,'  AveEnergy:',E15.6)
            write(*,649)P(iPC),totAcceptanceRatio/dble(iEqui),deltaRx
            write(1,649)P(iPC),totAcceptanceRatio/dble(iEqui),deltaRx
      649   format(1x,'P-pa:',E11.4, ' Ratio:',E15.5, '  deltaRx:',f10.6)

         endif
         
         
         do iMove = 1, Nmove
      
            call random_number(rdn)
            if (rdn .gt. 2.0E0/3.0E0)then
            call mcMove(TEnergy,Energy,CLEnergy,EnergyC, numberOfMove, numberOfAcceptanceMove)
            else
            call mcexc(TEnergy,Energy, CLEnergy,EnergyC, numofinsert,&
            &          successinsert,numofdelete,successdelete)
            endif
                         
            numberofequilibrium = numberofequilibrium + 1.0e0
            sum_npart = sum_npart + npart
            if(Npart .ne. 0)then
            EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
            endif

        enddo
        
        call adjustment(iEqui,numberOfMove, numberOfAcceptanceMove,AcceptanceRatio)
        totAcceptanceRatio = totAcceptanceRatio + AcceptanceRatio
       
     enddo
!    -----------------------------
!    SAMPLING STEPS FOR BULKPHASE
!    -----------------------------
     call ConfigEnergy(Energy, CLEnergy,EnergyC)
     TEnergy = Energy + CLEnergy 
     write(*,*) 'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy
     write(1,*) 'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy
     
     EnergyTotal            = 0.0E0
     numberOfMove           = 0
     numberOfAcceptanceMove = 0
     numberOfSampling       = 0.0E0    
     numofinsert            = 0
     successinsert          = 0
     numofdelete            = 0
     successdelete          = 0 
     totAcceptanceRatio     = 0.0
     
     Sum_TEnergy            = 0.0E0
     Sum_Npart              = 0.0E0
     Sum_TEN                = 0.0E0
     Sum_NNpart             = 0.0E0
     
         
     do iCycle=1, NCycles
     
         if( mod(iCycle,1000).eq.0) then
         
            write(*,*)'==================BULK SAMPLING====================='
            write(1,*)'==================BULK SAMPLING====================='
            write(*,716)iCycle, sum_Npart/numberOfSampling, EnergyTotal/numberOfSampling
            write(1,716)iCycle, sum_Npart/numberOfSampling, EnergyTotal/numberOfSampling
      716   format(1x,'CYCLES:',I8,'  Npart:',f15.6,'  AveEnergy:',e15.6)
            write(*,719)P(iPC),totAcceptanceRatio/dble(icycle),deltaRx
            write(1,719)P(iPC),totAcceptanceRatio/dble(icycle),deltaRx
      719   format(1x,'P-pa:',E11.4, 'Ratio:',e15.5, '  deltaRx:',f10.6)
         endif
           
         do iMove=1, Nmove           
             call random_number(rdn)
             
             if (rdn .gt. 2.0E0/3.0E0)then
                call  mcMove(TEnergy,Energy,CLEnergy,EnergyC, numberOfMove, numberOfAcceptanceMove) 
             else
                call  mcexc(TEnergy,Energy,CLEnergy, EnergyC, numofinsert,&
            &          successinsert,numofdelete,successdelete)
             endif
                       
             numberOfSampling = numberOfSampling + 1.0E0
                
             if(Npart .ne. 0)then
             EnergyTotal = EnergyTotal + TEnergy/dble(Npart)
             endif
                
             Sum_TEnergy  = Sum_TEnergy + TEnergy !!
             Sum_Npart    = Sum_Npart + dble(Npart) !!
             Sum_TEN      = Sum_TEN + TEnergy*dble(Npart) !!
             Sum_NNpart   = Sum_NNpart + dble(Npart)*dble(Npart) !!
             
         enddo
         
!         call adjustment(numberOfMove, numberOfAcceptanceMove,AcceptanceRatio) 
          totAcceptanceRatio = totAcceptanceRatio + AcceptanceRatio 
         successinsert =  0
         successdelete =  0         
     enddo
!    -----------------------------------------
!    FINISH THE SAMPLING; NOW DO THE AVERAGING 
!    -----------------------------------------
     EnergyAverage = EnergyTotal/numberOfSampling
     Av_Npart      = sum_npart / numberOfSampling
     Nb_rho        = Av_Npart / BulkVolume   
     
     Iso_TEnergy   = Sum_TEnergy/numberOfSampling !!
     Iso_Npart     = Sum_Npart/numberOfSampling !!
     Iso_TEN       = Sum_TEN/numberOfSampling !!
     Iso_NN        = Sum_NNpart/numberOfSampling  !!   
     
     Fun_NN_Bulk(iPC)  = Iso_NN - Iso_Npart*Iso_Npart
     N_Bulk(iPC)       = Iso_Npart
     Fun_EN_Bulk(iPC)  = Iso_TEN - Iso_TEnergy*Iso_Npart
     
    enddo
    
    RETURN
    END SUBROUTINE BULK_GCMC  
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE BUCKINGHAM_Exp6
!      Reference: Professor Do's NVT program
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine BuckinghamExp6
       USE FLUID_SOLID_M
       USE  BUCKINGHAM_M
       implicit none
       real*8 a_B, b_B, c_B, d_B
       real*8 aa_B, bb_B, cc_B, dd_B
              
       a_B  = 0.576864e0
       b_B  = 0.047980e0
       c_B  = -2.467966e-3
       d_B  = 4.541375e-5
       
       aa_B = 2.713425e0
       bb_B = -0.387513e0
       cc_B = 0.019765e0
       dd_B = -3.501142e-4
       
       rM = sigmaFF(1)/( a_B + b_B*alpha_B + c_B*alpha_B**2 + d_B*alpha_B**3 )
       rMin = ( aa_B + bb_B*alpha_B + cc_B*alpha_B**2 + dd_B*alpha_B**3 )*rM
       return
       end SUBROUTINE BuckinghamExp6
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE CrowellChangPotential
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine CrowellChangPotential(i,Energy)
       USE Constant_M
       USE FLUID_SOLID_M
       USE MCSETTING_M
       USE BUCKINGHAM_M
       implicit none
    
       integer i,j, k 
       
       real*8 Energy, iEnergy
       real*8 distanceZ
       real*8 reducedX, Fun_pentagamma, term1, term2
              
       Energy = 0.0e0    
       
       do j=1,NLJ
            distanceZ      = abs(LJz(i,j))
            if(distanceZ .lt. GrapheneLayer/5.0e0)then
               iEnergy = 1.0e20
            elseif(distanceZ.gt.rcutoff)then
               iEnergy = 0.0e0
            else
              reducedX       = distanceZ/GrapheneLayer
              Fun_pentagamma = 6.0e0*( 1/reducedX**4.0E0 + 1/3.0E0/(reducedX+0.6E0)**3.0E0 )
              term1          = -ASF_B*Pi*eta/12.0e0/GrapheneLayer**4.0e0
              term2          = 2.0e0*pi*eta*BSF_B*(alphaSF_B*distanceZ+1.0E0)/exp(alphaSF_B*distanceZ)/alphaSF_B**2.0e0
              iEnergy        = term1*Fun_pentagamma + term2
!              if(abs(iEnergy).gt.1e3)then
!              write(*,*)'distanceZ,reducedX,Fun_pentagamma,term1,term2,iEnergy',distanceZ,reducedX,Fun_pentagamma,term1,term2,iEnergy
!              write(*,*)'eta,bsf_b,alphasf_b,alphaSF_B*distanceZ,exp(alphaSF_B*distanceZ)',eta,bsf_b,alphasf_b,alphaSF_B*distanceZ,exp(alphaSF_B*distanceZ)
!              read(*,*)
!              endif
                
            endif
          Energy = Energy + iEnergy
       enddo
 
       return
       end SUBROUTINE CrowellChangPotential


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Fluctuation
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     subroutine fluctuation(switch)
     USE CONSTANT_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE FLUCTUATION_M
     implicit none
     
     integer i, j, k, bin, switch
     integer m, n
     integer bini,binj
     
     real*8 term0, term1, term2, term3, term4
     real*8 distanceZ
     real*8 NinBin_gas
     real*8 VolBin, Ave_N, Ave_NN, Fun_NN
     real*8 AccuU1tobin, AccuG1tobin,  Accuq1tobin
     real*8 AccuqFF1tobin, AccuqSF1tobin,  AccuqFFKL1tobin
     real*8 AccuU1tobinL1, AccuG1tobinL1,  Accuq1tobinL1
     real*8 AccuqFF1tobinL1, AccuqSF1tobinL1,  AccuqFFKL1tobinL1
     real*8 AccuU1tobinL2, AccuG1tobinL2,  Accuq1tobinL2
     real*8 AccuqFF1tobinL2, AccuqSF1tobinL2,  AccuqFFKL1tobinL2
     real*8, Dimension(:),POINTER:: dz, dzA, DDensity,DDensityKmolperM3
     real*8, Dimension(:),POINTER:: sum_Ninbin, sum_UFFinBin, sum_USFinBin
     real*8, Dimension(:),POINTER:: sum_NNinBin, sum_UFFNinBin, sum_USFNinBin
     real*8, Dimension(:),POINTER:: SumUFFKLofBin, LocalEnergy 
     real*8, Dimension(:),POINTER:: UFF1toBin, USF1toBin, UFFout1toBin, N1toBin
     real*8, Dimension(:),POINTER:: sum_SumUFFKLofBin, sum_SumUFFKLNofBin
     real*8, Dimension(:),POINTER:: sum_UFF1toBin, sum_USF1toBin, sum_UFFout1toBin, sum_N1toBin 
     real*8, Dimension(:),POINTER:: sum_UFFN1toBin, sum_USFN1toBin, sum_UFFNout1toBin, sum_NN1toBin
     real*8, Dimension(:),POINTER:: ave_Ninbin, ave_NNinBin
     real*8, Dimension(:),POINTER:: Ave_UFFinBin, Ave_USFinBin, Ave_UFFNinBin, Ave_USFNinBin 
     real*8, Dimension(:),POINTER:: ave_SumUFFKLofBin,  ave_SumUFFKLNofBin
     real*8, Dimension(:),POINTER:: ave_UFF1toBin, ave_USF1toBin, ave_UFFout1toBin, ave_N1toBin
     real*8, Dimension(:),POINTER:: ave_UFFN1toBin, ave_USFN1toBin, ave_UFFNout1toBin, ave_NN1toBin
     real*8, Dimension(:),POINTER:: Fun_NN_bin
     real*8, Dimension(:),POINTER:: Fun_UFFNinBin, Fun_USFNinBin, Fun_UFFKLNofBin
     real*8, Dimension(:),POINTER:: Fun_UFFN1toBin, Fun_USFN1toBin, Fun_UFFNout1toBin, Fun_NN1toBin 
     real*8, Dimension(:),POINTER:: LocalNFluctbin, LocalCompressbin, AccNFluctbin  
     real*8, Dimension(:),POINTER:: LocalEFluctBin, AccEFluct1toBin
     real*8, Dimension(:),POINTER:: Sum_LocalE, Sum_LocalEN, Ave_LocalE, Ave_LocalEN, Fun_LocalEN
     real*8, Dimension(:),POINTER:: Flu_LocalENex, Flu_LocalEUFF, Flu_LocalEUSF, Flu_LocalEUFFKL
     real*8, Dimension(:),POINTER:: Local_ISOH, Local_ISOHUSF, Local_ISOHUFF, Local_ISOHUFFKL
     real*8, Dimension(:),POINTER:: LocalE 
     real*8, Dimension(:,:),POINTER:: UFFKL
     real*8 AccUSF1toBin
     real*8 KJ, KJpermol
     real*8 LayerB1, LayerB2,  sumlocalE,sumlocalEffi,sumlocalEffj,sumlocalEsf
  
     logical steele, Bojan 
         
    
     
       KJ       = ScaleEnergy*kB/1.0e3  !! Used to convert the reduced unit to KJ
       KJpermol = ScaleEnergy*Rg/1.0e3  !! Used to convert the reduced unit to KJ/mol
       
       deltaZ   = 0.1e0/Scalelength
       maxbin   = int(BoxLengthZ/deltaZ) 
       if(mod(BoxLengthZ,deltaZ).gt. 1e-10)maxbin=maxbin+1
       VolBin   = CarbonLengthx1*CarbonLengthy1*deltaZ !!* Volume of each bin(-)
       LayerB1  = BoxlengthZ/2.0e0
       LayerB2  = BoxlengthZ
     if(switch.eq.0)then
       Allocate(dz(maxbin),dzA(maxbin),DDensity(maxbin),DDensityKmolperM3(maxbin))
       Allocate(UFFinBin(maxbin),USFinBin(maxbin))
       Allocate(sum_Ninbin(maxbin), sum_UFFinBin(maxbin), sum_USFinBin(maxbin))
       Allocate(sum_NNinBin(maxbin), sum_UFFNinBin(maxbin), sum_USFNinBin(maxbin))
       Allocate(SumUFFKLofBin(maxbin), LocalEnergy(maxbin))
       Allocate(UFF1toBin(maxbin), USF1toBin(maxbin), UFFout1toBin(maxbin), N1toBin(maxbin))
       Allocate(sum_SumUFFKLofBin(maxbin), sum_SumUFFKLNofBin(maxbin))
       Allocate(sum_UFF1toBin(maxbin), sum_USF1toBin(maxbin), sum_UFFout1toBin(maxbin), sum_N1toBin(maxbin))
       Allocate(sum_UFFN1toBin(maxbin), sum_USFN1toBin(maxbin), sum_UFFNout1toBin(maxbin), sum_NN1toBin(maxbin))
       Allocate(ave_Ninbin(maxbin), ave_NNinBin(maxbin))
       Allocate(Ave_UFFinBin(maxbin), Ave_USFinBin(maxbin), Ave_UFFNinBin(maxbin), Ave_USFNinBin(maxbin))
       Allocate(ave_SumUFFKLofBin(maxbin), ave_SumUFFKLNofBin(maxbin))
       Allocate(ave_UFF1toBin(maxbin), ave_USF1toBin(maxbin), ave_UFFout1toBin(maxbin), ave_N1toBin(maxbin))
       Allocate(ave_UFFN1toBin(maxbin), ave_USFN1toBin(maxbin), ave_UFFNout1toBin(maxbin), ave_NN1toBin(maxbin))
       Allocate(Fun_NN_bin(maxbin))
       Allocate(Fun_UFFNinBin(maxbin), Fun_USFNinBin(maxbin), Fun_UFFKLNofBin(maxbin))
       Allocate(Fun_UFFN1toBin(maxbin), Fun_USFN1toBin(maxbin), Fun_UFFNout1toBin(maxbin), Fun_NN1toBin(maxbin))
       Allocate(LocalNFluctbin(maxbin), LocalCompressbin(maxbin),AccNFluctbin(maxbin))
       Allocate(LocalEFluctBin(maxbin), AccEFluct1toBin(maxbin))
       Allocate(Sum_LocalE(maxbin), Sum_LocalEN(maxbin), Ave_LocalE(maxbin), Ave_LocalEN(maxbin), Fun_LocalEN(maxbin))
       Allocate(Flu_LocalENex(maxbin), Flu_LocalEUFF(maxbin), Flu_LocalEUSF(maxbin), Flu_LocalEUFFKL(maxbin))
       Allocate(Local_ISOH(maxbin), Local_ISOHUSF(maxbin), Local_ISOHUFF(maxbin), Local_ISOHUFFKL(maxbin))
       Allocate(LocalE(maxbin))
       Allocate(UFFKL(maxbin,maxbin))
       NinBin         = 0.0e0
       UFFinBin       = 0.0e0
       USFinBin       = 0.0e0 
       EnergyMatrix   = 0.0e0
        IF(LOCALF)THEN
          open(unit=35, file='LayerISOH.txt',status='unknown')
          rewind(35)
          write(35,'(9A14)')'P(pa)','layerB(1)','layerB(2)','L1ISOH(KJ/mol)','L1ISOHUFF','L1ISOHUSF','L1HUFFKL','L2ISOH(KJ/mol)','L2ISOHUFF','L2ISOHUSF','L2HUFFKL'
        ENDIF
       return
     endif

    if(switch.eq.1)then

       sum_N              = 0.0e0    !!* Total number of particle in simulation box 
       sum_NN             = 0.0E0    !!* N*N of the whole simulation box
       sum_Ninbin         = 0.0e0    !!* 
       sum_UFFinBin       = 0.0e0    !!* F-F interaction in Bin K
       sum_USFinBin       = 0.0e0    !!* S-F interaction of partticles in Bin K 
       sum_NNinBin        = 0.0e0    !!*
       sum_UFFNinBin      = 0.0e0    !!* N is Npart
       sum_USFNinBin      = 0.0e0    !!* N is Npart 
       sum_SumUFFKLofBin  = 0.0e0    !!* The contribution of "externel" F-F interaction of bin K with the other bins
       sum_SumUFFKLNofBin = 0.0e0    !!* N is Npart
 !      sum_UFF1toBin      = 0.0e0
 !      sum_USF1toBin      = 0.0e0
 !      sum_UFFout1toBin   = 0.0e0
       sum_N1toBin        = 0.0e0    !!*
 !      sum_UFFN1toBin     = 0.0e0
 !      sum_USFN1toBin     = 0.0e0
 !      sum_UFFNout1toBin  = 0.0e0
       sum_NN1toBin       = 0.0e0    !!*
       sum_LocalE         = 0.0e0    !!*
       sum_LocalEN        = 0.0e0    !!*
       CountN             = 0        !!*
       return
     endif
     
     if(switch.eq.2)then
        CountN = CountN + 1
        Sum_N  = Sum_N  + dble(Npart)
        Sum_NN = Sum_NN + dble(Npart)*dble(Npart)
        SumUFFKLofBin      = 0.0e0  !!*
        UFFKL              = 0.0e0
        UFF1toBin          = 0.0e0
        USF1toBin          = 0.0e0
        UFFout1toBin       = 0.0e0
        N1toBin            = 0.0e0 !!*
        sumlocalE = 0.0e0
        sumlocalEffi = 0.0e0
        sumlocalEffj = 0.0e0
        sumlocalEsf = 0.0e0
!         write(*,*)'fluctuation 2 CountN',CountN        
        
        do i = 1, Npart
           bini = I_POS(i)
           do j = 1, Npart
              if(j.NE.i)then
                 binj = I_POS(j)
                 if(bini.NE.binj)then
                     SumUFFKLofBin(bini) = SumUFFKLofBin(bini) + 0.5e0*EnergyMatrix(i,j)
                 endif
              endif
           enddo
         enddo
         
 !        do i = 1, Npart-1
 !           bini = I_POS(i)
 !           do j = i+1, Npart
 !              binj = I_POS(j)
 !              if(bini.NE.binj)then
 !                if(bini.LT.binj)then
 !                   UFFKL(bini,binj) = UFFKL(bini,binj) + EnergyMatrix(i,j)
 !                else
 !                   UFFKL(binj,bini) = UFFKL(binj,bini) + EnergyMatrix(i,j)
 !                endif
 !              endif
 !             enddo
 !          enddo  
               
 !        UFF1toBin(1) =  UFFinBin(1) 
 !        USF1toBin(1) =  USFinBin(1) 
         N1toBin(1)   =  NinBin(1)  !!*   
 !        UFFout1toBin(1) = 0.0e0   
 !        do i = 2,maxbin
 !           UFFout1toBin(1) =  UFFout1toBin(1) + UFFKL(1,i)
 !        enddo     
       do i = 2,maxbin
               
 !              UFF1toBin(i) = UFF1toBin(i-1) + UFFinBin(i)
 !              USF1toBin(i) = USF1toBin(i-1) + USFinBin(i)
               N1toBin(i)   = N1toBin(i-1) + NinBin(i)  !!*-*
 !              do j = 1, i-1
 !                 UFF1toBin(i) = UFF1toBin(i) +  UFFKL(j,i)
 !              enddo
               
 !              do j = i+1, maxbin
 !                 UFFout1toBin(i) = UFFout1toBin(i-1) + UFFKL(i,j)
 !              enddo
 !              do m = 1, i-1
 !                 UFFout1toBin(i) = UFFout1toBin(i) - UFFKL(m,i)
 !              enddo
        enddo
        
        do i = 1, maxbin
          LocalE(i)             = UFFinBin(i) + USFinBin(i) + SumUFFKLofBin(i)
          sum_Ninbin(i)         = sum_Ninbin(i)    + dble(Ninbin(i)) !!*
          sum_UFFinBin(i)       = sum_UFFinBin(i)  + UFFinBin(i) !!*
          sum_USFinBin(i)       = sum_USFinBin(i)  + USFinBin(i) !!*
          sum_NNinBin(i)        = sum_NNinBin(i)   + NinBin(i)*NinBin(i)  !!* 
          sum_UFFNinBin(i)      = sum_UFFNinBin(i) + UFFinBin(i)*Npart    !!*
          sum_USFNinBin(i)      = sum_USFNinBin(i) + USFinBin(i)*Npart    !!*    
          sum_SumUFFKLofBin(i)  = sum_SumUFFKLofBin(i) + SumUFFKLofBin(i) !!* 
          sum_SumUFFKLNofBin(i) = sum_SumUFFKLNofBin(i) + SumUFFKLofBin(i)*Npart  !!*
          sum_LocalE(i)         = sum_LocalE(i) + LocalE(i)   !UFFinBin(i) + USFinBin(i) + SumUFFKLofBin(i) !!*
          sum_LocalEN(i)        = sum_LocalEN(i) + LocalE(i)*Npart  !!*
 !         sum_UFF1toBin(i)      = sum_UFF1toBin(i) + UFF1toBin(i) 
 !         sum_USF1toBin(i)      = sum_USF1toBin(i) + USF1toBin(i)
 !         sum_UFFout1toBin(i)   = sum_UFFout1toBin(i) + UFFout1toBin(i)
          sum_N1toBin(i)        = sum_N1toBin(i) + N1toBin(i) !!*
 !         sum_UFFN1toBin(i)     = sum_UFFN1toBin(i) + UFF1toBin(i)*N1toBin(i)
 !         sum_USFN1toBin(i)     = sum_USFN1toBin(i) + USF1toBin(i)*N1toBin(i)
 !         sum_UFFNout1toBin(i)  = sum_UFFNout1toBin(i) + UFFout1toBin(i)*N1toBin(i)
          sum_NN1toBin(i)       = sum_NN1toBin(i) + N1toBin(i)*N1toBin(i) !!*
          sumlocalE = sumlocalE + LocalE(i)
          sumlocalEffi = sumlocalEffi + UFFinBin(i)
          sumlocalEffj = sumlocalEffj + SumUFFKLofBin(i)
          sumlocalEsf = sumlocalEsf + USFinBin(i)
          
        enddo
!        write(*,*)'i,UFF1toBin(i),UFFout1toBin(i),USF1toBin(i)',maxbin, UFF1toBin(maxbin),UFFout1toBin(maxbin),USF1toBin(maxbin) 
!        read(*,*)    
!         write(*,*)'sumlocalE,sumlocalEffi,sumlocalEffj,sumlocalEsf', sumlocalE,sumlocalEffi,sumlocalEffj,sumlocalEsf
!         read(*,*)
     endif
     
     if(switch.eq.3)then
        open(unit=32, file='Fluctuation.txt',status='unknown')
        write(32,'(8A16)')'Pressure(pa)', 'distance(A)', 'LocalD(Kmol/m3)', 'LocNFluct(-)', 'AccuNFluct(-)', 'LocCompres(1/pa)', 'LocalE(KJ)', 'Accu_LocE(KJ)'
        open(unit=33, file='Flu_LocalENex.txt',status='unknown')
        write(33,'(7A16)')'Pressure(pa)','distance(A)', 'LocE_Flu(KJ/mol)','LocUFF_Flu', 'LocUSF_Flu', 'LocUFFKL_Flu', 'Accu_locE_Flu'
        open(unit=34, file='LocalISOH.txt',status='unknown')
   !     write(34,'(6A16)')'Pressure(pa)','distance(A)','LocISOH(KJ/mol)','LocISOHUFF','LocISOHUSF','LocISOHUFFKL'
        write(34,'(7A16,3A18)')'Pressure(pa)','distance(A)','LocISOH(KJ/mol)','LocISOHUFF','LocISOHUSF','LocISOHUFFKL','Accu_locISOH','Accu_LocISOHUFF','Accu_LocISOHUSF','Accu_LocISOHUFFKL'

        NinBin_gas = VolBin*rho(iPC)  !!* (molecular)
        ave_N  = sum_N/dble(CountN)   !!*
        ave_NN = sum_NN/dble(CountN)  !!*
        Fun_NN = ave_NN - ave_N*ave_N !!*
        AccuU1tobin     = 0.0e0
        AccuG1tobin     = 0.0e0
        Accuq1tobin     = 0.0e0
        AccuqFF1tobin   = 0.0e0
        AccuqSF1tobin   = 0.0e0
        AccuqFFKL1tobin = 0.0e0
        AccuU1tobinL1     = 0.0e0
        AccuG1tobinL1     = 0.0e0
        Accuq1tobinL1     = 0.0e0
        AccuqFF1tobinL1   = 0.0e0
        AccuqSF1tobinL1   = 0.0e0
        AccuqFFKL1tobinL1 = 0.0e0
        AccuU1tobinL2     = 0.0e0
        AccuG1tobinL2     = 0.0e0
        Accuq1tobinL2     = 0.0e0
        AccuqFF1tobinL2   = 0.0e0
        AccuqSF1tobinL2   = 0.0e0
        AccuqFFKL1tobinL2 = 0.0e0
        
        do i = 1, maxbin
            dz(i)  = deltaZ*(i-0.5)
            dzA(i) = dz(i)*scalelength
            ave_Ninbin(i)         = sum_Ninbin(i)/dble(CountN) !!*
            ave_UFFinBin(i)       = sum_UFFinBin(i)/dble(CountN) !!*
            ave_USFinBin(i)       = sum_USFinBin(i)/dble(CountN)  !!*
            ave_NNinBin(i)        = sum_NNinBin(i)/dble(CountN) !!*
            ave_UFFNinBin(i)      = sum_UFFNinBin(i)/dble(CountN) !!*
            ave_USFNinBin(i)      = sum_USFNinBin(i)/dble(CountN) !!*
            ave_SumUFFKLofBin(i)  = sum_SumUFFKLofBin(i)/dble(CountN) !!*
            ave_SumUFFKLNofBin(i) = sum_SumUFFKLNofBin(i)/dble(CountN) !!*
!            ave_UFF1toBin(i)      = sum_UFF1toBin(i)/dble(CountN)
!            ave_USF1toBin(i)      = sum_USF1toBin(i)/dble(CountN)
!            ave_UFFout1toBin(i)   = sum_UFFout1toBin(i)/dble(CountN)
            ave_N1toBin(i)        = sum_N1toBin(i)/dble(CountN) !!
!            ave_UFFN1toBin(i)     = sum_UFFN1toBin(i)/dble(CountN)
!            ave_USFN1toBin(i)     = sum_USFN1toBin(i)/dble(CountN)
!            ave_UFFNout1toBin(i)  = sum_UFFNout1toBin(i)/dble(CountN)
            ave_NN1toBin(i)       = sum_NN1toBin(i)/dble(CountN) !!*
            ave_LocalE(i)         = sum_LocalE(i)/dble(CountN)  !!*
            ave_LocalEN(i)        = sum_LocalEN(i)/dble(CountN)  !!*

            DDensity(i)           = ave_Ninbin(i)/(CarbonLengthx1*CarbonLengthy1*deltaZ) !!*
            DDensityKmolperM3(i)  = DDensity(i)/6.0230E23/(Scalelength*1.0e-10)**3/1.0e3 !!* 
            Fun_NN_bin(i)         = ave_NNinBin(i) - ave_Ninbin(i)*ave_Ninbin(i) !!*
            Fun_UFFNinBin(i)      = Ave_UFFNinBin(i) - Ave_UFFinBin(i)*ave_N   !!*
            Fun_USFNinBin(i)      = Ave_USFNinBin(i) - Ave_USFinBin(i)*ave_N   !!*
            Fun_UFFKLNofBin(i)    = ave_SumUFFKLNofBin(i) - ave_SumUFFKLofBin(i)*ave_N   !!*
  !          Fun_UFFN1toBin(i)     = ave_UFFN1toBin(i) - ave_UFF1toBin(i)*ave_N1toBin(i)
  !          Fun_USFN1toBin(i)     = ave_USFN1toBin(i) - ave_USF1toBin(i)*ave_N1toBin(i)  
  !          Fun_UFFNout1toBin(i)  = ave_UFFout1toBin(i) - ave_UFFout1toBin(i)*ave_N1toBin(i) 
            Fun_NN1toBin(i)       = ave_NN1toBin(i)-ave_N1toBin(i)*ave_N1toBin(i) !!*
            Fun_LocalEN(i)        = ave_LocalEN(i) - ave_LocalE(i)*ave_N  !!* f(Uk,N)
            LocalEnergy(i)        = ave_UFFinBin(i) + ave_USFinBin(i) + ave_SumUFFKLofBin(i)  !!* Local Energy of Interaction -- Uk
            Flu_LocalENex(i)      = -Fun_LocalEN(i)/(Fun_NN- CFun_NN_BULK(iPC)) !!* Fluctuation of local energy with the change in the total excess number --- Gk
            Flu_LocalEUFF(i)      = -Fun_UFFNinBin(i)/(Fun_NN- CFun_NN_BULK(iPC)) !!* Fluctuation of local energy with the change in the total excess number-UFF
            Flu_LocaLEUSF(i)      = -Fun_USFNinBin(i)/(Fun_NN- CFun_NN_BULK(iPC)) !!* Fluctuation of local energy with the change in the total excess number-USF
            Flu_LocalEUFFKL(i)    = -Fun_UFFKLNofBin(i)/(Fun_NN- CFun_NN_BULK(iPC)) !!* Fluctuation of local energy with the change in the total excess number-UFFKL 
            Local_ISOH(i)         = Flu_LocalENex(i) + NinBin_gas*Temperature/CFun_NN_BULK(iPC) !!* Local Isosteric heat -- qk
            
            Local_ISOHUFF(i)      = Flu_LocalEUFF(i) + NinBin_gas*Temperature/CFun_NN_BULK(iPC) !!* Local Isosteric heat -- UFF
            Local_ISOHUSF(i)      = Flu_LocalEUSF(i)   !!* Local Isosteric heat -- USF
            Local_ISOHUFFKL(i)    = Flu_LocalEUFFKL(i) !!* Local Isosteric heat -- UFFKL
            
           
 !         write(32,'(5E16.6)')ave_Ninbin(i), ave_UFFinBin(i), ave_USFinBin(i),ave_N1toBin(i),ave_NN1toBin(i)  
            if(ave_Ninbin(i).GT. 0.e0)then
               LocalNFluctbin(i)    = Fun_NN_bin(i)/ave_Ninbin(i)  !!*
               LocalCompressbin(i)  = Fun_NN_bin(i)/(Rg*T*DDensityKmolperM3(i))*1.0E-3/ave_Ninbin(i)    !(m3/j=1/pa) !!*
 !              LocalEFluctBin(i)    = (Fun_UFFNinBin(i)+Fun_UFFKLNofBin(i)+Fun_USFNinBin(i))/Fun_NN_bin(i)
            else
               LocalNFluctbin(i)    = 1.0e0 !!*
               LocalCompressbin(i)  = 1.0e0/P(iPC)    !(m3/j=1/pa) !!*
 !              LocalEFluctBin(i)    = ave_USFinBin(i)
            endif
          
            AccNFluctBin(i)    = Fun_NN1toBin(i)/ave_N1toBin(i)
!            AccEFluct1toBin(i) = (Fun_UFFN1toBin(i)+Fun_UFFNout1toBin(i)+ Fun_USFN1toBin(i) )/Fun_NN1toBin(i)
            if(ave_N1toBin(i) .EQ. 0.e0)then
!              AccUSF1toBin = 0.0e0
!              do j = 1, i
!                 AccUSF1toBin = AccUSF1toBin + ave_USFinBin(j)
!              enddo
!              AccEFluct1toBin(i) = AccUSF1toBin/dble(i)
              AccNFluctBin(i)    = 1.0e0
            endif
           
          !---------------------------------------------------
          !
          ! Accumulative of Uk, Gk, qk from 1st bin to ith bin 
          !
          !---------------------------------------------------
                AccuU1tobin     =  AccuU1tobin     +  LocalEnergy(i)
                AccuG1tobin     =  AccuG1tobin     +  Flu_LocalENex(i)
                Accuq1tobin     =  Accuq1tobin     +  Local_ISOH(i) 
                AccuqFF1tobin   =  AccuqFF1tobin   +  Local_ISOHUFF(i)
                AccuqSF1tobin   =  AccuqSF1tobin   +  Local_ISOHUSF(i)
                AccuqFFKL1tobin =  AccuqFFKL1tobin +  Local_ISOHUFFKL(i) 
              
                
          !---------------------------------------
          !
          ! Accumulative of Uk, Gk, qk for layers
          !
          !--------------------------------------
             if(dz(i).le.layerB(1))then
                AccuU1tobinL1     =  AccuU1tobinL1     +  LocalEnergy(i)
                AccuG1tobinL1     =  AccuG1tobinL1     +  Flu_LocalENex(i)
                Accuq1tobinL1     =  Accuq1tobinL1     +  Local_ISOH(i) 
                AccuqFF1tobinL1   =  AccuqFF1tobinL1   +  Local_ISOHUFF(i)
                AccuqSF1tobinL1   =  AccuqSF1tobinL1   +  Local_ISOHUSF(i)
                AccuqFFKL1tobinL1 =  AccuqFFKL1tobinL1 +  Local_ISOHUFFKL(i) 
             elseif(dz(i).le.layerB(2))then
                AccuU1tobinL2     =  AccuU1tobinL2     +  LocalEnergy(i)
                AccuG1tobinL2     =  AccuG1tobinL2     +  Flu_LocalENex(i)
                Accuq1tobinL2     =  Accuq1tobinL2     +  Local_ISOH(i) 
                AccuqFF1tobinL2   =  AccuqFF1tobinL2   +  Local_ISOHUFF(i)
                AccuqSF1tobinL2   =  AccuqSF1tobinL2   +  Local_ISOHUSF(i)
                AccuqFFKL1tobinL2 =  AccuqFFKL1tobinL2 +  Local_ISOHUFFKL(i) 
             endif
             
          
            write(32,'(8E16.6)')P(iPC), dzA(i),DDensityKmolperM3(i),LocalNFluctbin(i),AccNFluctBin(i), LocalCompressbin(i),LocalEnergy(i)*KJpermol, AccuU1tobin*KJpermol
            write(33,'(7E16.6)')P(iPC), dzA(i),Flu_LocalENex(i)*KJpermol, Flu_LocalEUFF(i)*KJpermol, Flu_LocaLEUSF(i)*KJpermol, Flu_LocalEUFFKL(i)*KJpermol, AccuG1tobin*KJpermol   
            write(34,'(10E16.6)')P(iPC), dzA(i),Local_ISOH(i)*KJpermol, Local_ISOHUFF(i)*KJpermol, Local_ISOHUSF(i)*KJpermol, Local_ISOHUFFKL(i)*KJpermol , Accuq1tobin*KJpermol, &
           &                   AccuqFF1tobin*KJpermol, AccuqSF1tobin*KJpermol, AccuqFFKL1tobin*KJpermol   

          enddo
            write(35,'(11E14.4)')P(iPC), layerB(1)*scalelength,layerB(2)*scalelength,Accuq1tobinL1*KJpermol, AccuqFF1tobinL1*KJpermol, AccuqSF1tobinL1*KJpermol, AccuqFFKL1tobinL1*KJpermol, &
           &                  Accuq1tobinL2*KJpermol, AccuqFF1tobinL2*KJpermol, AccuqSF1tobinL2*KJpermol, AccuqFFKL1tobinL2*KJpermol  
!          do i = 1,maxbin
!          write(32,'(8E16.7)')Fun_NN_bin(i), Fun_UFFNinBin(i),Fun_USFNinBin(i),Fun_UFFKLNofBin(i), Fun_UFFN1toBin(i),Fun_USFN1toBin(i),Fun_UFFNout1toBin(i), Fun_NN1toBin(i)
!          enddo
          
!          do i = 1,maxbin
!          write(32,'(7E16.7)')P(iPC), dzA(i),DDensityKmolperM3(i),LocalNFluctbin(i), LocalCompressbin(i),LocalEFluctBin(i),AccEFluct1toBin(i) 
!          enddo
        endif   ! endif (swith.eq.3)
        
        if(switch.eq.4)then
           DEAllocate(dz,dzA,DDensity,DDensityKmolperM3)
           DEAllocate(UFFinBin,USFinBin)
           DEAllocate(sum_Ninbin, sum_UFFinBin, sum_USFinBin)
!           DEAllocate(sum_NNinBin, sum_UFFNinBin, sum_USFNinBin)
           DEAllocate(SumUFFKLofBin,LocalEnergy)
           DEAllocate(UFF1toBin, USF1toBin, UFFout1toBin, N1toBin)
           DEAllocate(sum_SumUFFKLofBin, sum_SumUFFKLNofBin)
           DEAllocate(sum_UFF1toBin, sum_USF1toBin, sum_UFFout1toBin, sum_N1toBin)
           DEAllocate(sum_UFFN1toBin, sum_USFN1toBin, sum_UFFNout1toBin, sum_NN1toBin)
           DEAllocate(ave_Ninbin, ave_NNinBin)
           DEAllocate(Ave_UFFinBin, Ave_USFinBin, Ave_UFFNinBin, Ave_USFNinBin )
           DEAllocate(ave_SumUFFKLofBin,  ave_SumUFFKLNofBin)
           DEAllocate(ave_UFF1toBin, ave_USF1toBin, ave_UFFout1toBin, ave_N1toBin)
           DEAllocate(ave_UFFN1toBin, ave_USFN1toBin, ave_UFFNout1toBin, ave_NN1toBin)
           DEAllocate(Fun_NN_bin)
           DEAllocate(Fun_UFFNinBin, Fun_USFNinBin, Fun_UFFKLNofBin)
           DEAllocate(Fun_UFFN1toBin, Fun_USFN1toBin, Fun_UFFNout1toBin, Fun_NN1toBin )
           DEAllocate(LocalNFluctbin, LocalCompressbin, AccNFluctbin)
           DEAllocate(LocalEFluctBin, AccEFluct1toBin)
           DEAllocate(Sum_LocalE, Sum_LocalEN, Ave_LocalE, Ave_LocalEN, Fun_LocalEN)
           DEAllocate(Flu_LocalENex, Flu_LocalEUFF, Flu_LocalEUSF, Flu_LocalEUFFKL)
           DEAllocate(Local_ISOH, Local_ISOHUSF, Local_ISOHUFF, Local_ISOHUFFKL)
           DEAllocate(LocalE)
           DEAllocate(UFFKL)
        endif
     
   return
   end
