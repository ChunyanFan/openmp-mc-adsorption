!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE kMC_VLE
!    FOR THE kMC VAPOR-LIQUID EQUILIBRIUM
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       subroutine kMC_VLE
       USE Constant_M
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       USE ROTATE_M
       implicit none
    
       integer NumberOfMove, NumberOfAcceptanceMove
       integer iEqui, iMove, iCycle
       integer i, j, fSub
       integer numberofequilibrium
       integer iSub, iS, iRot,iSRot
       integer FLAG1, FLAG2, FLAG3
    
       real*8 secnds, time,rdn, sumpair
       real*8 EnergyTotal,totAcceptanceRatio
       real*8 Energy, EnergyC, CLEnergy, TEnergy, EnergyAverage
       real*8 numberOfSampling
       real*8 rhoMolPerM3, PressurePa, EnergyAverageJoule
       real*8 AcceptanceRatio, EnergyLRC, VirialLRC
       real*8 VirialTotal 
       real*8 PGas1, PGas2, PLiq, PLRC(3)  
       real*8 PGas1Pa, PGas2Pa, PLiqPa
       real*8 AverhoMolPerM3(3), VirialPropTotal(3), VirialPropAverage(3)
       real*8 totNinPropSub(3),AveNinPropSub(3),PropVol(3)
       real*8 EnergyPropTotal(3),EnergyPropAverage(3)
       real*8 EvapHeat1,EvapHeat2,EvapHeatKJPerMole1,EvapHeatKJPerMole2
       real*8 FirstDerivETotal, FirstDerivEAve, SurfaceTension, SurfaceTensionJPerM2
       real*8 FirstDerivETotal_Ale, FirstDerivEAve_Ale, SurfaceTension_Ale, SurfaceTensionJPerM2_Ale
       real*8 TEnergyNew, EnergyNew, CLEnergyNew
       real*8 IdChemPkVLE, ExChemPkVLE, ChemPkVLE, ChemPkVLEJPerMole
       real*8 TEnergyOld, EnergyOld 
       real*8, Dimension(:),POINTER::PairEnergyOld, PairEnergyNew
       
       common/SAMPSUM1/TEnergy,numberOfSampling, EnergyTotal,FirstDerivETotal, FirstDerivETotal_Ale
       common/SAMPSUM2/totNinPropSub, VirialPropTotal, EnergyPropTotal   
    
      
!    ------------------------------------------------------  
!    Set rCutOff which can't larger than BoxLength/2
!    ------------------------------------------------------
         rcutoff =  VLEBoxlengthX/2.0E0
         if(rcutoff.GT.VLEBoxlengthY/2.0E0)then
            rcutoff = VLEBoxlengthY/2.0E0
         endif
         if(rcutoff.GT.VLEBoxlengthZ/2.0E0)then
            rcutoff = VLEBoxlengthZ/2.0E0
         endif
         if(rcutoff.gt.5.0)rcutoff=5.0e0
         r2CutOff = rCutOff*rCutOff
         
         do i =1,3
           rCutoffprop(i)=VLEBoxlengthX/2.0D0
           if(rCutoffprop(i).gt.(PropHy(i)-PropLy(i))/2.0D0)rCutoffprop(i)=(PropHy(i)-PropLy(i))/2.0D0
           if(rCutoffprop(i).gt.VLEBoxlengthZ/2.0D0)rCutoffprop(i)=VLEBoxlengthZ/2.0D0
           r2Cutoffprop(i) = rCutoffprop(i)*rCutoffprop(i)
         enddo

         BulkVolume  = VLEBoxlengthX*VLEBoxlengthY*VLEBoxlengthZ
         TotSurfaceArea = VLEBoxlengthX*VLEBoxlengthZ
    
         write(62,1130) TotSurfaceArea, BulkVolume
         write(51,1130) TotSurfaceArea, BulkVolume
1130     format( 1x, 'SurfaceArea(-)--------------------------:', f20.10,/, &
  &              1x, 'BulkVolume(-)---------------------------:', f20.10,/)

       
       open(unit=25, file='kVLEFile.txt',status='unknown')
       rewind(25)
       write(25,*)'------------------------------------------------------------'
       write(25,'(4x,A18,1A14)')'NumofEqulibrim','rho_SBox(i)'
       write(25,'(4x,2A14)')'(-)', '(mol/m3)'   
       write(25,*)'------------------------------------------------------------'
       write(25,111)SubBox,(i,i=1,SubBox)
111    FORMAT(1x,I15,200I8) 
!       ------------------------
!       STARTING CONFIGURATION
!       ------------------------
        NinSBox = 0
        SeqNinSBox = 0
        IndNinSBox = 0
        NinPropSub = 0 
        SeqNinPropSub = 0
        indNinPropSub = 0
        IF(SolidSlab)THEN
           CALL IniSolidSlab
        ELSEIF(HomoStart)THEN
           if(Npart.eq.0)Npart = Int(LiqDens*VLEBoxlengthX*VLEBoxlengthY*VLEBoxlengthZ)+1
           CALL InsertVLE
        ELSE
           if(Npart.eq.0)Npart = Int(LiqDens*VLEBoxlengthX*(LiqHy-LiqLy)*VLEBoxlengthZ)+1
           CALL InsertVLE
        ENDIF
        ALLOCATE(rate(Npart), AccuR(Npart), Energyi(Npart))
        ALLOCATE(PairEnergyOld(Npart),PairEnergyNew(Npart),PairEnergyi(Npart))  
        ALLOCATE(EnergyMatrix(Npart,Npart),VirialMatrix(Npart,Npart),DerivEMatrix(Npart,Npart))
        ALLOCATE(DerivEMatrix_Ale(Npart,Npart))
        
        EnergyMatrix = 0.0E0
        VirialMatrix = 0.0E0
        DerivEMatrix = 0.0E0
        DerivEMatrix_Ale = 0.0E0        
!       ------------
!       START TIMING
!       ------------
        time = secnds(0.0)

!       -----------------
!       EQUILIBRIUM STEPS
!       -----------------
        numberOfMove           = 0
        numberOfAcceptanceMove = 0
        numberofequilibrium    = 0
        VirialTotal            = 0.0
        totAcceptanceRatio     = 0.0 
        EnergyTotal            = 0.0
        totNinSBox             = 0.0
        TotTime                = 0.0
        
        call ConfigEnergy(Energy, CLEnergy)
        TEnergy = Energy 
        write(62,*)'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy
        write(51,*)'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy
     

        DO iEqui = 1, Ncycle
            if (mod(iEqui,1000).eq.0)then
               write(62,*)'====================kVLE EQULIBRATION==================='
               write(51,*)'====================kVLE EQULIBRATION==================='
               write(62,646)iEqui, EnergyTotal/TotTime/dble(Npart),Npart
               write(51,646)iEqui, EnergyTotal/TotTime/dble(Npart),Npart
      646      format(1x,'CYCLES:',I8, '  AveEnergy:',E15.6, '   Npart:',I5)
            endif
            

!           -----------------------------------------------------
!           CHECK THE CENTER OF THE FLUID AND RESET THE PARTICLE 
!           ----------------------------------------------------            
            if(mod(iEqui,1000).eq.0)then
               call LiqPosCheck
               call ConfigEnergy(Energy, CLEnergy)
               TEnergy = Energy 
               IF(ChangeXZ)THEN
                  call ShapeChange(TEnergy, Energy,CLEnergy,TEnergyNew, EnergyNew,CLEnergyNew)
                  TEnergy = TEnergyNew
                  Energy = EnergyNew
                  call random_number(rdn)
                  InterTime = log(1.0/rdn)/TotRate   
                  TotTime = TotTime + InterTime
                  EnergyTotal = EnergyTotal + TEnergy*InterTime          
                  do i = 1, SubBox
                      totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
                  enddo 
               ENDIF
            endif
!           ---------------------------
!           Randomly pick one particle
!           ---------------------------
            call select(iS)
            CALL ThermalPropertyChange(1,iS)
                        
!              -------------------------
!              Pick a position randomly 
!              -------------------------
               call random_number(rdn)
               MCx(iS) = rdn*VLEBoxlengthX
               call random_number(rdn)
               MCy(iS) = rdn*VLEBoxlengthY
               call random_number(rdn)
               MCz(iS) = rdn*VLEBoxlengthZ 
            
               CALL ROTATE
               do j = 1, NLJ
                  LJx(iS,j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j) + MCx(iS) - MCx0
                  LJy(iS,j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j) + MCy(iS) - MCy0
                  LJz(iS,j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j) + MCz(iS) - MCz0
               enddo

               do j = 1, NCL
                  CLx(iS,j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) + MCx(iS) - MCx0
                  CLy(iS,j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) + MCy(iS) - MCy0
                  CLz(iS,j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) + MCz(iS) - MCz0
               enddo     

                   !------------------------------------------
                   ! Update the location information in Bins
                   !------------------------------------------ 
                    fSub = indNinSBox(iS)
                    do i =1, NinSBox(fSub)
                       if(SeqNinSBox(fSub,i).eq.iS)then
                         SeqNinSBox(fSub,i)= SeqNinSBox(fSub,NinSBox(fSub))
                         EXIT
                       endif
                    enddo
                    NinSBox(fSub) = NinSBox(fSub)-1
                 
                    do iSub=1, SubBox
                        if((MCy(iS).GE.LayerYL(iSub)).AND.(MCy(iS).LT.LayerYH(iSub)))then
                           NinSBox(iSub)= NinSBox(iSub)+1
                           SeqNinSBox(iSub,NinSBox(iSub))=iS
                           indNinSBox(iS)=iSub
                           EXIT
                        endif
                    enddo
                 
                   !--------------------
                   ! Update the Energy
                   !--------------------
                    TEnergy = TEnergy - Energyi(iS) 
                    Energy  = Energy - Energyi(iS)             
                    CALL PairEnergy(iS)
                    TEnergy = TEnergy + Energyi(iS) 
                    Energy  = Energy + Energyi(iS)             
                    CALL ThermalPropertyChange(2,iS)
            
            call random_number(rdn)
            InterTime = log(1.0/rdn)/TotRate   
            TotTime = TotTime + InterTime
            EnergyTotal = EnergyTotal + TEnergy*InterTime
                        
            do i = 1, SubBox
                 totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
             enddo 
        ENDDO !iEqui
             do i = 1, SubBox
                 ave_NinSBox(i) = totNinSBox(i)/TotTime
             enddo 
        write(25,'(5X,100E16.6)')(ave_NinSBox(i),i=1,SubBox)
        
        write(62,*) 'L710 TEnergy, Energy, CLEnergy',TEnergy, Energy, CLEnergy
        write(51,*) 'L710 TEnergy, Energy, CLEnergy',TEnergy, Energy, CLEnergy

!       --------------
!       SAMPLING STEPS
!       --------------
        call ConfigEnergy(Energy, CLEnergy)
        TEnergy = Energy  
        write(62,*) 'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy
        write(51,*) 'TEnergy, Energy, CLEnergy', TEnergy, Energy, CLEnergy

     
        EnergyTotal            = 0.0E0
        numberOfSampling       = 0.0E0    
        VirialPropTotal        = 0
        EnergyPropTotal        = 0
        ExChemPVLETotal        = 0.0
        IdChemPVLETotal        = 0.0 
        ChemPVLETotal          = 0.0
        FirstDerivETotal       = 0.0
        FirstDerivETotal_Ale   = 0.0
        TotTime                = 0.0
  
        IF(VLEradius)call VLEradiusdistribution(1)
        IF(VLEdistance)call VLEdistancedistribution(1)
 !       call FLUCTVLE(0)
 !       call FLUCTVLE(1)
        
         totNinSBox     = 0.0
         totNinPropSub  = 0.0
         Nmove=SubBox*(SubBox-1)/2
         
         do iCycle = 1, Ncycles
              if (mod(iCycle,1000).eq.0)then
                 write(62,*)'================kVLE SAMPLING==================='
                 write(51,*)'================kVLE SAMPLING==================='
                 write(62,711)iCycle, EnergyTotal/TotTime/dble(Npart), Npart
                 write(51,711)iCycle, EnergyTotal/TotTime/dble(Npart), Npart
      711        format(1x,'CYCLES:',I8, '  AveEnergy:',E15.6, '   Npart:',I5)
              endif
!             -----------------------------------------------------
!             CHECK THE CENTER OF THE FLUID AND RESET THE PARTICLE 
!             ----------------------------------------------------            
              if (mod(iCycle,1000).eq.0)then
                 call LiqPosCheck
                 call ConfigEnergy(Energy, CLEnergy)
                 TEnergy = Energy 
                 IF(ChangeXZ)THEN
                    call ShapeChange(TEnergy, Energy,CLEnergy,TEnergyNew, EnergyNew,CLEnergyNew)
                    TEnergy = TEnergyNew
                    Energy = EnergyNew
                    CALL ENSEMBLESUM   
                 ENDIF
              endif           
!             ----------------------------------
!             Pick one particle with Rosenbluth
!             ----------------------------------
              CALL select(iS)
              CALL ThermalPropertyChange(1,iS)
              
!               -------------------------
!               Pick a position randomly 
!               -------------------------
                 call random_number(rdn)
                 MCx(iS) = rdn*VLEBoxlengthX
                 call random_number(rdn)
                 MCy(iS) = rdn*VLEBoxlengthY
                 call random_number(rdn)
                 MCz(iS) = rdn*VLEBoxlengthZ 
            
                 CALL ROTATE
                 do j = 1, NLJ
                    LJx(iS,j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j) + MCx(iS) - MCx0
                    LJy(iS,j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j) + MCy(iS) - MCy0
                    LJz(iS,j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j) + MCz(iS) - MCz0
                 enddo

                 do j = 1, NCL
                    CLx(iS,j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) + MCx(iS) - MCx0
                    CLy(iS,j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) + MCy(iS) - MCy0
                    CLz(iS,j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) + MCz(iS) - MCz0
                 enddo     

                    !------------------------------------------
                    ! Update the location information in Bins
                    !------------------------------------------ 
                     fSub = indNinSBox(iS)
                     do i =1, NinSBox(fSub)
                        if(SeqNinSBox(fSub,i).eq.iS)then
                          SeqNinSBox(fSub,i)= SeqNinSBox(fSub,NinSBox(fSub))
                          EXIT
                        endif
                     enddo
                     NinSBox(fSub) = NinSBox(fSub)-1
                 
                     do iSub=1, SubBox
                        if((MCy(iS).GE.LayerYL(iSub)).AND.(MCy(iS).LT.LayerYH(iSub)))then
                            NinSBox(iSub)= NinSBox(iSub)+1
                          SeqNinSBox(iSub,NinSBox(iSub))=iS
                          indNinSBox(iS)=iSub
                          EXIT
                       endif
                     enddo              
                    !--------------------
                    ! Update the Energy
                    !--------------------
                    TEnergy = TEnergy - Energyi(iS) 
                    Energy  = Energy - Energyi(iS)             
                    CALL PairEnergy(iS)
                    TEnergy = TEnergy + Energyi(iS) 
                    Energy  = Energy + Energyi(iS)             
                    CALL ThermalPropertyChange(2,iS)


              CALL ENSEMBLESUM

              
              IF(VLEradius)call VLEradiusdistribution(2)
              IF(VLEdistance)call VLEdistancedistribution(2)   
 !             call FLUCTVLE(2)
         enddo  !iCycle
!      -----------------------------------------
!      FINISH THE SAMPLING; NOW DO THE AVERAGING 
!      -----------------------------------------
       EnergyAverage = EnergyTotal/TotTime
       FirstDerivEAve  = FirstDerivETotal/TotTime
       SurfaceTension  = FirstDerivEAve/2.0/(VLEBoxlengthX*VLEBoxlengthZ)
       SurfaceTensionJPerM2 = SurfaceTension*SCALEENERGY*kB/((SCALELENGTH*1.0E-10)**2.0)
       FirstDerivEAve_Ale  = FirstDerivETotal_Ale/TotTime
       SurfaceTension_Ale  = FirstDerivEAve_Ale/2.0/(VLEBoxlengthX*VLEBoxlengthZ)
       SurfaceTensionJPerM2_Ale = SurfaceTension_Ale*SCALEENERGY*kB/((SCALELENGTH*1.0E-10)**2.0)
       
       do i=1,3
          VirialPropAverage(i) = VirialPropTotal(i)/TotTime
          AveNinPropSub(i) = totNinPropSub(i)/TotTime
          EnergyPropAverage(i) = EnergyPropTotal(i)/AveNinPropSub(i)/TotTime
          write(25,*)'9131 i VirialPropTotal(i),VirialPropAverage(i)',i, VirialPropTotal(i),VirialPropAverage(i)
       enddo
       
       do i=1, SubBox
          ave_NinSBox(i) = totNinSBox(i)/TotTime
          rho_SBox(i)  = ave_NinSBox(i)/SubVol(i) 
          rho_SBoxMolPerM3(i) = rho_SBox(i)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
       enddo
       
       IdChemPkVLE = Temperature*LOG(deBro3/(VLEBoxlengthX*VLEBoxlengthY*VLEBoxlengthZ))
       ExChemPkVLE = Temperature*LOG(numberOfSampling/TotTime)
       ChemPkVLE = IdChemPkVLE + ExChemPkVLE
       ChemPkVLEJPerMole = ChemPkVLE * SCALEENERGY*Rg
       
       
       if(WIDOM_VLE)then
         do i=1, SubBox
             ExChemPVLEAve(i) = ExChemPVLETotal(i)/dble(Ncycles)
             IdChemPVLEAve(i) = IdChemPVLETotal(i)/dble(Ncycles)
             ChemPVLEAve(i)   = ChemPVLETotal(i)/dble(Ncycles)
             ExChemPVLEJPerMole(i)=ExChemPVLEAve(i)*SCALEENERGY*Rg
             IdChemPVLEJPerMole(i) = IdChemPVLEAve(i)*SCALEENERGY*Rg
             ChemPVLEJPerMole(i) = ChemPVLEAve(i)*SCALEENERGY*Rg
         enddo
       endif
  
       PropVol(1) = VLEBoxlengthX*(PropHy(1)-PropLy(1))*VLEBoxlengthZ
       AveProprho(1) = AveNinPropSub(1)/PropVol(1)
       PGas1 = AveProprho(1)*Temperature + VirialPropAverage(1)/PropVol(1)
       call ZLRC(AveProprho(1), EnergyLRC, VirialLRC)
       PLRC(1) = AveProprho(1)*VirialLRC
       PGas1 = PGas1 + PLRC(1)
       AverhoMolPerM3(1) = AveProprho(1)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
       PGAS1Pa = PGas1*SCALEENERGY*kB/(SCALELENGTH*1.0E-10)**3
       PLRC(1) = PLRC(1)*SCALEENERGY*kB/(SCALELENGTH*1.0E-10)**3
       
       
       PropVol(2) = VLEBoxlengthX*(PropHy(2)-PropLy(2))*VLEBoxlengthZ
       AveProprho(2) = AveNinPropSub(2)/PropVol(2)
       PLiq = AveProprho(2)*Temperature + VirialPropAverage(2)/PropVol(2)
       call ZLRC(AveProprho(2), EnergyLRC, VirialLRC)
       PLRC(2) = AveProprho(2)*VirialLRC
       PLiq = PLiq + PLRC(2)
       AverhoMolPerM3(2) = AveProprho(2)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
       PLiqPa = PLiq*SCALEENERGY*kB/(SCALELENGTH*1.0E-10)**3
       PLRC(2) = PLRC(2)*SCALEENERGY*kB/(SCALELENGTH*1.0E-10)**3

       
       PropVol(3) = VLEBoxlengthX*(PropHy(3)-PropLy(3))*VLEBoxlengthZ
       AveProprho(3) = AveNinPropSub(1)/PropVol(3)
       PGas2 = AveProprho(3)*Temperature + VirialPropAverage(3)/PropVol(3)
       call ZLRC(AveProprho(3), EnergyLRC, VirialLRC)
       PLRC(3) = AveProprho(3)*VirialLRC
       PGas2 = PGas2 + PLRC(3)
       AverhoMolPerM3(3) = AveProprho(3)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
       PGAS2Pa = PGas2*SCALEENERGY*kB/(SCALELENGTH*1.0E-10)**3
       PLRC(3) = PLRC(3)*SCALEENERGY*kB/(SCALELENGTH*1.0E-10)**3
       
       IF(VLEradius)call VLEradiusdistribution(3)
       IF(VLEdistance)call VLEdistancedistribution(3)
       
       EvapHeat1 = EnergyPropAverage(1)- EnergyPropAverage(2) + PGas1*(1.0/AveProprho(1) - 1.0/AveProprho(2))
       EvapHeat2 = EnergyPropAverage(3)- EnergyPropAverage(2) + PGas2*(1.0/AveProprho(3) - 1.0/AveProprho(2))
       
       EvapHeatKJPerMole1 = EvapHeat1*SCALEENERGY*Rg/1.0E3
       EvapHeatKJPerMole2 = EvapHeat2*SCALEENERGY*Rg/1.0E3
            
       write(25,*)'------BOX DIMENSIONS(A)-----------------------------------------'
       write(25,'(5X,3E16.6)')VLEBoxlengthX0,VLEBoxlengthY0,VLEBoxlengthZ0
       write(25,'(5X,3E16.6)')VLEBoxlengthX,VLEBoxlengthY,VLEBoxlengthZ
       write(25,'(5X,3E16.6)')VLEBoxlengthX*SCALELENGTH,VLEBoxlengthY*SCALELENGTH,VLEBoxlengthZ*SCALELENGTH
       write(25,*)'-------------------sampling results-----------------------------'  
       write(25,'(5X,100E16.6)')((LayerYL(i)+(LayerYH(i)-LayerYL(i))/2.0)*SCALELENGTH, i=1,SubBox)
       write(25,'(5X,100E16.6)')(rho_SBoxMolPerM3(i),i=1,SubBox)
       write(25,'(5X,100E16.6)')(ave_NinSBox(i),i=1,SubBox)
       write(25,*)'-------------------CHEMICAL POTENTIAL kMC-----------------------' 
       write(25,'(5X,4E16.6)')ChemPkVLEJPerMole,IdChemPkVLE*SCALEENERGY*Rg,ExChemPkVLE* SCALEENERGY*Rg, (ChemPkVLE-Temperature*LOG(deBro3))/Temperature
       write(25,*)'-------------------CHEMICAL POTENTIAL widom-----------------------' 
       write(25,554)Ncycles,(ChemPVLEJPerMole(i),i=1,SubBox)
       write(25,*)'----------------IDEAL CHEMICAL POTENTIAL widom--------------------' 
       write(25,554)Ncycles,(IdChemPVLEJPerMole(i),i=1,SubBox)
       write(25,*)'----------------EXCESS CHEMICAL POTENTIAL widom-------------------' 
       write(25,554)Ncycles,(ExChemPVLEJPerMole(i),i=1,SubBox)
       write(25,*)'-------------------1st Gas Phase--------------------------------'  
       write(25,'(4E18.6)')AverhoMolPerM3(1), PGAS1Pa, PGAS1Pa-PLRC(1), PLRC(1)
       write(25,*)'-------------------Liquid Phase---------------------------------'  
       write(25,'(4E18.6)')AverhoMolPerM3(2), PLiqPa, PLiqPa-PLRC(2), PLRC(2)
       write(25,*)'-------------------2nd Gas Phase--------------------------------'  
       write(25,'(4E18.6)')AverhoMolPerM3(3), PGAS2Pa, PGAS2Pa-PLRC(3), PLRC(3)
       write(25,*)'-------------------Evaporation Heat(KJ/mol)--------------------'  
       write(25,'(2E18.6)')EvapHeatKJPerMole1,EvapHeatKJPerMole2
       write(25,*)'-------------------Surface tension(- & J/m2)--------------------'  
       write(25,'(2E18.6)')SurfaceTension,SurfaceTensionJPerM2
       write(25,*)'---------------Surface tension_Ale(- & J/m2)--------------------'  
       write(25,'(2E18.6)')SurfaceTension_Ale,SurfaceTensionJPerM2_Ale

554    format(I20,200E18.6)       
       
!       call FLUCTVLE(3)
!       call FLUCTVLE(4)
!      ------------------
!      DISPLAY THE OUTPUT
!      ------------------
       time = secnds(0.0) - time
        
       write(25,50) time            
50     format(1x,'Time (sec) --------------------------------- : ', f20.10/)
!        -------------------------------
!        STORE THE POSITION OF PARTICLES
!        -------------------------------
         write(52,*) 'Npart=', Npart
         do i = 1, Npart
            write(2,436) mcx(i), mcy(i), mcz(i)
436         format(3f16.7)
         enddo  
         write(52,*) 'LJSITES' 
         do i = 1, Npart
            do j = 1, NLJ
               write(52,436) LJX(i,j),LJY(i,j),LJZ(i,j)
            enddo
          enddo
         write(52,*) 'COULOMBSITES' 
         do i = 1, Npart
            do j = 1, NCL
               write(52,436) CLX(i,j),CLY(i,j),CLZ(i,j)
            enddo
         enddo
  
        DEALLOCATE(LayerZL,LayerZH, SubVol)
        DEALLOCATE(IniStepX,IniStepY,IniStepZ)
        DEALLOCATE(SStepX,SStepY,SStepZ)
        DEALLOCATE(SeqNinSBox,totNinSBox,ave_NinSBox)
        DEALLOCATE(AtpinMove, SucinMove,SucMovein,SucMoveout)  
        DEALLOCATE(LayerXL,LayerXH,LayerYL,LayerYH)
        DEALLOCATE(rho_SBox,rho_SBoxMolPerM3)
        DEALLOCATE(rhoFW,ChemicalPVLE, ExcessChemicalPVLE, IdealChemicalPVLE)         
        DEALLOCATE(ChemPVLETotal, ExChemPVLETotal, IdChemPVLETotal)
        DEALLOCATE(ChemPVLEAve, ExChemPVLEAve, IdChemPVLEAve)
        DEALLOCATE(ChemPVLEJPerMole, ExChemPVLEJPerMole, IdChemPVLEJPerMole)  
        DEALLOCATE(rate, AccuR, Energyi)
        DEALLOCATE(PairEnergyOld,PairEnergyNew,PairEnergyi)
        DEALLOCATE(EnergyMatrix,VirialMatrix,DerivEMatrix, DerivEMatrix_Ale)
        

     RETURN
     END SUBROUTINE kMC_VLE
     
 
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE MC_SURFACE
!    FOR THE Kinetic MC for adsorpiton on surface
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       subroutine MC_SURFACE
       USE Constant_M
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       USE ROTATE_M
       implicit none
    
       integer NumberOfMove, NumberOfAcceptanceMove
       integer iMove, iCycle
       integer i, j, fSub, cutoffopt
       integer numberofequilibrium
       integer iSub, iS, iRot,iSRot
       integer FLAG1, FLAG2, FLAG3
       integer iSOld, OverZ
    
       real*8 secnds, time,rdn, iEqui
       real*8 EnergyTotal
       real*8 Energy, EnergyC, CLEnergy, TEnergy, EnergyAverage
       real*8 numberOfSampling, BulkVirialTotal,BulkEnergyTotal
       real*8 BulkVirialAve, SEinputuMolPerM2
       real*8 rhoMolPerM3, PressurePa, EnergyAverageJoule
       real*8 EnergyLRC, VirialLRC
       real*8 VirialTotal, PLRC
       real*8 TEnergyNew, EnergyNew, CLEnergyNew
       real*8 IdChemPkVLE, ExChemPkVLE, ChemPkVLE, ChemPkVLEJPerMole
       real*8 TEnergyOld, EnergyOld 
       real*8 rcutoffset
       real*8 EinBulkAve, CNBulk, CEBulk, CEBulkJoule  
       real*8 TotMTtime, InterTimeOld, Ave_MTtime, PresentTime
       real*8 LowUsys0, HighUsys0, kineticE
                    
!    ------------------------------------------------------  
!    Set rCutOff which can't larger than BoxLength/2
!    ------------------------------------------------------
          minBoxLength = BoxLengthX
         if(minBoxLength .gt. BoxLengthY)minBoxLength = BoxLengthY
         if(minBoxLength .gt. BoxLengthZ)minBoxLength = BoxLengthZ
         if(pbcx) rcutoff =  BoxLengthX/2.0E0
         if(pbcy .AND. rcutoff.GT.BoxLengthY/2.0E0)rcutoff = BoxLengthY/2.0E0
         if(pbcz .AND. rcutoff.GT.BoxLengthZ/2.0E0)rcutoff = BoxLengthZ/2.0E0
     
          rcutoffset = 5.0E0
          
          if(rcutoff.gt.rcutoffset)then
            write(62,*)'---------------------------------------- '
            write(62,*)'Which one would be chosen as the rcutoff:'
            write(62,*)'---------------------------------------- '
            write(62,*)'1-', rcutoff
            write(62,*)'2-', rcutoffset
            if(NLJ.GT.1)then
               cutoffopt=1
            else
               cutoffopt = 2
            endif
            if(cutoffopt.eq.2)rcutoff = rcutoffset
          endif
     
          IF((.NOT. pbcx) .and. (.NOT. pbcy) .and. (.NOT. pbcz))rcutoff = rcutoffset
          IF(Bojan)rcutoff = BoxLengthX/2.0E0
     
          write(62,*)'rCutOff=',rCutOff   
          write(51,*)'rCutOff',rCutOff
          r2CutOff = rCutOff*rCutOff

          BulkVolume  = BoxlengthX*BoxlengthY*BoxlengthZ
          TotSurfaceArea = BoxlengthX*BoxlengthY
          if(Nsteele .gt. 1)then
             TotSurfaceArea  = Nsteele*BoxlengthX*BoxlengthY
          endif
          write(62,*)'------------------------------------- '
          write(62,*)'Calculating the Access Volume:        '
          write(62,*)'------------------------------------- '  
          call accessiblevolume
          AccVolumeM3  = AccVolume*(SCALELENGTH*1.0E-10)**3
    
          write(62,1130) TotSurfaceArea, BulkVolume, AccVolume
          write(51,1130) TotSurfaceArea, BulkVolume, AccVolume
1130      format( 1x, 'SurfaceArea(-)--------------------------:', f20.10,/, &
  &              1x, 'BulkVolume(-)---------------------------:', f20.10,/,&
  &              1x, 'AccVolume(-)----------------------------:', f20.10,/)  
       

        open(unit=57, file='SurfIsotherm.txt',status='unknown')
        rewind(57)
        write(57,*)'------------------------------------------------------------------------------------------------'
        write(57,'(A4,14A14)') 'no.','TotPickN','OutZPickN','rCutoff','rho','P','Npart','SurfE','N-Nbulk','U-Ubulk','ChemicalP','S','S/kB','Ave_NinBin', 'ChemP_Bin'
        write(57,'(5x,14A14)')'(-)','(-)','(-)','(mol/m3)','(pa)','(-)','(umol/m2)','(-)','(Joule)','(J/mol)','(-)','J/(Kmol)','(-)','J/(Kmol)'
        write(57,*)'-------------------------------------------------------------------------------------------------'
!       ------------------------
!       STARTING CONFIGURATION
!       ------------------------
        ALLOCATE(rate(5000), AccuR(5000), Energyi(5000))
        ALLOCATE(EnergyMatrix(5000,5000),VirialMatrix(5000,5000))
        Npart = Npart !+ int( SEinput(1)*TotSurfaceArea)+1
        NinSBox = 0
        SeqNinSBox = 0
        IndNinSBox = 0
        CALL InsertSurf(0,Npart)
        IF(DistanceD)call distancedistribution(0)
        
        DO iden = 1, NumOfrho
            TotPickN = 0
            OutZPickN = 0
            WRITE(62,*)'INSERT PARTICLES'
            CALL InsertSurf(Npart,AddN)

            Npart = Npart + AddN   
            write(51,*)'Npart=',Npart
            write(62,*)'Npart=',Npart
            EnergyMatrix = 0.0E0
            VirialMatrix = 0.0E0
            IF(MassTSFER)THEN
              UpFlux = 0.0
              DownFlux = 0.0
              TotMTtime = 0.0
            ENDIF
!           ------------
!           START TIMING
!           ------------
            time = secnds(0.0)

!           -----------------
!           EQUILIBRIUM STEPS
!           -----------------
            numberofequilibrium    = 0
            VirialTotal            = 0.0
            EnergyTotal            = 0.0
            totNinSBox             = 0.0
            TotTime                = 0.0
            PresentTime            = 0.0
            WRITE(62,*)'CALCULATE THE CONFIGENERGY'
            call ConfigEnergyADS(Energy,CLEnergy, EnergyC)
            TEnergy = Energy + EnergyC
            write(62,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
            write(51,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
            call random_number(rdn)
            InterTimeOld = log(1.0/rdn)/TotRate  
            IF(ENTROPY)THEN
              LowUsys0 = LowUi*Npart
              HighUsys0 = HighUi*Npart
              LowUsys = TEnergy
              HighUsys = TEnergy
            ENDIF

            iSOld = 0
            OverZ = -1
            DO iEqui = 1, Ncycle
                if (mod(iEqui,100000.0).eq.0)then
                   write(62,*)'====================kMCSurf EQULIBRATION==================='
                   write(51,*)'====================kMCSurf EQULIBRATION==================='
                   write(62,646)int(iEqui), EnergyTotal/TotTime/dble(Npart),Npart
                   write(51,646)int(iEqui), EnergyTotal/TotTime/dble(Npart),Npart
      646          format(1x,'CYCLES:',I12, '  AveEnergy:',E15.6, '   Npart:',I5)
                endif
                
!               -----------------------------
!               ADJUST THE SIMULATION BOX
!               -----------------------------         
                if(mod(iEqui,100.0).eq.0)then
                   IF(ChangeXY)THEN
                      call ShapeChange(TEnergy, Energy,CLEnergy,TEnergyNew, EnergyNew,CLEnergyNew)
                      TEnergy = TEnergyNew
                      call random_number(rdn)
                      InterTime = log(1.0/rdn)/TotRate   
                      TotTime = TotTime + InterTime
                      EnergyTotal = EnergyTotal + TEnergy*InterTime          
                      do i = 1, SubBox
                          totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
                      enddo 
                   ENDIF
                endif
 !               DO iMove = 1,SubBox    
!                  ---------------------------
!                  Randomly pick one particle
!                  ---------------------------
                   IF(NBOURBIN)THEN
                      call NBOURSelect(iS)
                   ELSE
                      if((iSOld.NE.0).and. (OverZ.EQ.1))then
                         iS = iSOld
                      else
                         call select(iS)
                      endif
                   ENDIF
                   CALL TherPropChangeADS(1,iS)
                        
!                  ---------------------
!                  Pick a NEW position  
!                  ---------------------
                 IF(MeanFPath)THEN
                      call MFPMCpick(iS,OverZ)
                      IF(OverZ.EQ.1)iSOld = iS
                      IF(OverZ.EQ.-1)iSOld = 0
                 ELSE
                      call SurfMCpick(iS)
                 ENDIF
                 
                      CALL ROTATE
                      do j = 1, NLJ
                         LJx(iS,j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j) + MCx(iS) - MCx0
                         LJy(iS,j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j) + MCy(iS) - MCy0
                         LJz(iS,j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j) + MCz(iS) - MCz0
                      enddo

                      do j = 1, NCL
                         CLx(iS,j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) + MCx(iS) - MCx0
                         CLy(iS,j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) + MCy(iS) - MCy0
                         CLz(iS,j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) + MCz(iS) - MCz0
                      enddo     

                          !------------------------------------------
                          ! Update the location information in Bins
                          !------------------------------------------ 
                           fSub = indNinSBox(iS)
                           do i =1, NinSBox(fSub)
                              if(SeqNinSBox(fSub,i).eq.iS)then
                                SeqNinSBox(fSub,i)= SeqNinSBox(fSub,NinSBox(fSub))
                                EXIT
                              endif
                           enddo
                           NinSBox(fSub) = NinSBox(fSub)-1
                           do iSub=1, SubBox
                               if((MCz(iS).GE.LayerZL(iSub)).AND.(MCz(iS).LT.LayerZH(iSub)))then
                                  NinSBox(iSub)= NinSBox(iSub)+1
                                  SeqNinSBox(iSub,NinSBox(iSub))=iS
                                  indNinSBox(iS)=iSub
                                  EXIT
                               endif
                           enddo
                          !--------------------
                          ! Update the Energy
                          !--------------------
                           TEnergy = TEnergy - Energyi(iS) 
                           CALL PairEnergyADS(iS)
                           TEnergy = TEnergy + Energyi(iS) 
                           CALL TherPropChangeADS(2,iS)
            
                   call random_number(rdn)
                   InterTime = log(1.0/rdn)/TotRate   
                   TotTime = TotTime + InterTime
                   EnergyTotal = EnergyTotal + TEnergy*InterTime
                   IF(MassTSFER)THEN  
                      IF(fSub.NE.iSub)CALL FLUXCALCULATE(fSub, iSub)
                      Ave_MTtime = (InterTimeOld + InterTime)/2.0
                      TotMTtime = TotMTtime + Ave_MTtime
                      InterTimeOld = InterTime
                   ENDIF
                   IF(ENTROPY)THEN
                      if(TEnergy.gt.LowUsys)LowUsys = TEnergy
                      if(TEnergy.lt.HighUsys)HighUsys = TEnergy
                   ENDIF    
                   do i = 1, SubBox
                     totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
                   enddo 
    !           ENDDO !iMove 
                !    IF(MassTSFER .AND. mod(iEqui,100).EQ.0 .AND. (iEqui.le.1000000))THEN
                   IF(MassTSFER .AND. modulo(LOG10(iEqui),1.0).EQ.0)THEN               
                      do i =1,NFluxSurf
                         NetFlux(i) = (DownFlux(i)-UpFlux(i))/TotMTtime
                      enddo
                      do i = 1, SubBox
                         ave_NinSBox(i) = totNinSBox(i)/TotTime
                      enddo 
                      PresentTime = PresentTime + TotTime/2.0e0
                      write(11,'(I10,500E16.6)')int(iEqui),PresentTime,(ave_NinSBox(i)/SubVol(i),i=1,SubBox) !, (NetFlux(i),i=1,NFluxSurf)
                      PresentTime = PresentTime + TotTime/2.0e0
                      DownFlux = 0.0
                      UpFlux   = 0.0
                      TotMTtime  = 0.0
                      totNinSBox = 0.0
                      TotTime = 0.0                     
                   ENDIF                   
            ENDDO !iEqui
                 do i = 1, SubBox
                     ave_NinSBox(i) = totNinSBox(i)/TotTime
                 enddo 

            write(62,*) 'EndEqulibrium TEnergy',TEnergy
            write(51,*) 'EndEqulibrium TEnergy',TEnergy
!           --------------
!           SAMPLING STEPS
!           --------------
            call ConfigEnergyADS(Energy, CLEnergy, EnergyC)
            TEnergy = Energy + EnergyC  
            write(62,*) 'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
            write(51,*) 'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC

     
            TotPickN = 0
            OutZPickN = 0
            EnergyTotal            = 0.0E0
            numberOfSampling       = 0.0E0    
            BulkVirialTotal        = 0.0E0
            BulkEnergyTotal        = 0.0E0
            ExChemPVLETotal        = 0.0
            IdChemPVLETotal        = 0.0 
            ChemPVLETotal          = 0.0
            TotTime                = 0.0
            rateTSBox              = 0.0
            timeSBox               = 0.0
            
            IF(ENTROPY)THEN
              write(62,*)'1476',LowUsys0,HighUsys0
              write(62,*)'1477',LowUsys,HighUsys
              write(51,*)'1476',LowUsys0,HighUsys0
              write(51,*)'1477',LowUsys,HighUsys
              IF(LowUsys.GT.LowUsys0) LowUsys=LowUsys0
              IF(HighUsys.LT.HighUsys0) HighUsys=HighUsys0
              write(62,*)'1482',LowUsys,HighUsys
    !          CALL EntropyCalculation(1)
            ENDIF
  
!            IF(VLEradius)call VLEradiusdistribution(1)
            IF(DistanceD)call distancedistribution(1)
        
             totNinSBox     = 0.0
             Nmove=SubBox*(SubBox-1)/2
             
             iSOld = 0
             OverZ = -1

             do iCycle = 1, Ncycles
                  if (mod(iCycle,100000).eq.0)then
                     write(62,*)'================kMCSurf SAMPLING==================='
                     write(51,*)'================kMCSurf SAMPLING==================='
                     write(62,711)iCycle, EnergyTotal/TotTime/dble(Npart), Npart
                     write(51,711)iCycle, EnergyTotal/TotTime/dble(Npart), Npart
      711            format(1x,'CYCLES:',I12, '  AveEnergy:',E15.6, '   Npart:',I5)
                  endif
!                 -----------------------------------------------------
!                 CHECK THE CENTER OF THE FLUID AND RESET THE PARTICLE 
!                 ----------------------------------------------------            
                  if (mod(iCycle,100).eq.0)then
                     IF(ChangeXY)THEN
                        call ShapeChange(TEnergy, Energy,CLEnergy,TEnergyNew, EnergyNew,CLEnergyNew)
                        TEnergy = TEnergyNew
                        Energy = EnergyNew
                        numberOfSampling  = numberOfSampling + 1.0
                        call random_number(rdn)
                        InterTime = log(1.0/rdn)/TotRate   
                        TotTime = TotTime + InterTime
                        EnergyTotal = EnergyTotal + TEnergy*InterTime 
                        BulkVirialTotal = BulkVirialTotal + BulkVirial*InterTime 
                        BulkEnergyTotal = BulkEnergyTotal + EinBulk*InterTime     
                        do i = 1, SubBox
                           totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
                        enddo 
              !          IF(ENTROPY)CALL EntropyCalculation(2,TEnergy)
                     ENDIF
                  endif           
    !              DO iMove =1, SubBox
!                    ----------------------------------
!                    Pick one particle with Rosenbluth
!                    ----------------------------------
                     IF(NBOURBIN)THEN
                        call NBOURSelect(iS)
                     ELSE
                        if((iSOld.NE.0).and. (OverZ.EQ.1))then
                           iS = iSOld
                        else
                           call select(iS)
                        endif
                     ENDIF

                     CALL TherPropChangeADS(1,iS)
              
!                      -------------------------
!                      Pick a position randomly 
!                      -------------------------
                        IF(MeanFPath)THEN
                            call MFPMCpick(iS,OverZ)
                            IF(OverZ.EQ.1)iSOld = iS
                            IF(OverZ.EQ.-1)iSOld = 0
                        ELSE
                            call SurfMCpick(iS)
                        ENDIF
            
                        CALL ROTATE
                        do j = 1, NLJ
                           LJx(iS,j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j) + MCx(iS) - MCx0
                           LJy(iS,j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j) + MCy(iS) - MCy0
                           LJz(iS,j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j) + MCz(iS) - MCz0
                        enddo

                        do j = 1, NCL
                           CLx(iS,j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) + MCx(iS) - MCx0
                           CLy(iS,j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) + MCy(iS) - MCy0
                           CLz(iS,j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) + MCz(iS) - MCz0
                        enddo     

                           !------------------------------------------
                           ! Update the location information in Bins
                           !------------------------------------------ 
                            fSub = indNinSBox(iS)
                            do i =1, NinSBox(fSub)
                               if(SeqNinSBox(fSub,i).eq.iS)then
                                 SeqNinSBox(fSub,i)= SeqNinSBox(fSub,NinSBox(fSub))
                                 EXIT
                               endif
                            enddo
                            NinSBox(fSub) = NinSBox(fSub)-1
                 
                            do iSub=1, SubBox
                               if((MCz(iS).GE.LayerZL(iSub)).AND.(MCz(iS).LT.LayerZH(iSub)))then
                                 NinSBox(iSub)= NinSBox(iSub)+1
                                 SeqNinSBox(iSub,NinSBox(iSub))=iS
                                 indNinSBox(iS)=iSub
                                 EXIT
                              endif
                            enddo              
                           !--------------------
                           ! Update the Energy
                           !--------------------
                           TEnergy = TEnergy - Energyi(iS) 
                           CALL PairEnergyADS(iS)
                           TEnergy = TEnergy + Energyi(iS) 
                           CALL TherPropChangeADS(2,iS)

                           numberOfSampling  = numberOfSampling + 1.0
                           call random_number(rdn)
                           InterTime = log(1.0/rdn)/TotRate   
                           TotTime = TotTime + InterTime
                           EnergyTotal = EnergyTotal + TEnergy*InterTime
                           BulkVirialTotal = BulkVirialTotal + BulkVirial*InterTime  
                           BulkEnergyTotal = BulkEnergyTotal + EinBulk*InterTime             
                           do i = 1, SubBox
                              totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
                              rateTSBox(i) = rateTSBox(i) + rateSBox(i)*InterTime
                              timeSBox(i) = timeSBox(i) + InterTime
                           enddo 
     !                      IF(ENTROPY)CALL EntropyCalculation(2,TEnergy) 
 !                    IF(VLEradius)call VLEradiusdistribution(2)
                     IF(DistanceD)call distancedistribution(2)   
 !               ENDDO ! iMove 
             enddo  !iCycle
!          -----------------------------------------
!          FINISH THE SAMPLING; NOW DO THE AVERAGING 
!          -----------------------------------------
           EnergyAverage = EnergyTotal/TotTime   !/dble(Npart)
           BulkVirialAve = BulkVirialTotal/TotTime
           EinBulkAve = BulkEnergyTotal/TotTime
           
              
           do i=1, SubBox
              ave_NinSBox(i) = totNinSBox(i)/TotTime
              rho_SBox(i)  = ave_NinSBox(i)/SubVol(i) 
              rho_SBoxMolPerM3(i) = rho_SBox(i)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
           enddo

           write(11,'(A10,500E16.6)')'Sampling',TotTime,(rho_SBox(i),i=1,SubBox) 
       
           IF(NLJ.GT.1)THEN
              IdChemPkVLE = Temperature*( LOG(deBro3) + LOG(2*thetar/T) + thetav/2.0/T + LOG(1-exp(-thetav/T)) - ChemPDe/kB/T ) + Temperature*LOG(1/(BoxlengthX*BoxlengthY*BoxlengthZ))
              WRITE(1,*)LOG(deBro3) , LOG(2*thetar/T) , thetav/2.0/T , LOG(1-exp(-thetav/T)) , ChemPDe/kB/T
           ELSE
              IdChemPkVLE = Temperature*LOG(deBro3/(BoxlengthX*BoxlengthY*BoxlengthZ))
           ENDIF
           ExChemPkVLE = Temperature*LOG(numberOfSampling/TotTime)
           ChemPkVLE = IdChemPkVLE + ExChemPkVLE
           ChemPkVLEJPerMole = ChemPkVLE * SCALEENERGY*Rg
           do i=1, SubBox
              ChemPSBox(i) = Temperature*LOG(deBro3*rateTSBox(i)/timeSBox(i)/SubVol(i)) * SCALEENERGY*Rg
           enddo

           
           rho(iden) = rho_SBox(SubFGas)
           IF((.NOT. Steele) .or. SigmaSS.EQ.0.0)rho(iden)=dble(Npart)/BulkVolume
           P(iden) = rho(iden)*Temperature + BulkVirialAve/SubVol(SubFGas)
           CNBulk = rho(iden)*AccVolume
           CEBulk = EinBulkAve*AccVolume/SubVol(SubFGas)
           CEBulkJoule = CEBulk*SCALEENERGY*kB
           call ZLRC(rho(iden), EnergyLRC, VirialLRC)
           PLRC = rho(iden)*VirialLRC
           P(iden) = P(iden) + PLRC
           SEinput(iden) = (Npart-rho(iden)*AccVolume)/totSurfaceArea
!           IF(ENTROPY)CALL EntropyCalculation(3) 
           IF(NLJ.GT.1)THEN
              kineticE = dble(Npart)*Temperature*(2.5E0 + thetav/T/2.0 + thetav/T/(exp(thetav/T) - 1)  - ChemPDe/kB/T)
              WRITE(1,*)thetav/T/2.0, thetav/T/(exp(thetav/T) - 1), ChemPDe/kB/T
              EntropySys = (kineticE+EnergyAverage+P(iden)*BulkVolume-dble(Npart)*ChemPkVLE)/Temperature/dble(Npart)*Rg  !(J/(K*mol))
           ELSE
              kineticE = 1.5*dble(Npart)*Temperature
              EntropySys = (1.5*dble(Npart)*Temperature+EnergyAverage+P(iden)*BulkVolume-dble(Npart)*ChemPkVLE)/Temperature/dble(Npart)*Rg  !(J/(K*mol))
           ENDIF
 !           EntropySys = (1.5*dble(Npart)*Temperature*kB*ScaleEnergy + EnergyAverage*SCALEENERGY*kB + P(iden)*SCALEENERGY*kB*BulkVolume &
 !        &  -dble(Npart)*ChemPkVLE*SCALEENERGY*kB)/(Temperature*ScaleEnergy)/dble(Npart)*AvogadroNumber
           
!          ---------------------- 
!          DIMENSIONAL QUANTITIES
!          ----------------------
           rhoMolPerM3              = rho(iden)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
           SEinputuMolPerM2         = SEinput(iden)/(SCALELENGTH*1.0E-10)**2/AvogadroNumber*1.0e6
           PressurePa               = P(iden)*SCALEENERGY*kB/(SCALELENGTH*1E-10)**3
           EnergyAverageJoule       = EnergyAverage*SCALEENERGY*kB 
     
 !          IF(VLEradius)call VLEradiusdistribution(3)
            IF(DistanceD)call distancedistribution(3)
              
!          ------------------
!          DISPLAY THE OUTPUT
!          ------------------
           time = secnds(0.0) - time
           write(62,50)  AccVolume,Npart,P(iden),&
           &          PLRC, EnergyAverage,&
           &          SEinputuMolPerM2,PressurePa, AccVolumeM3,EnergyAverageJoule,&
           &          time            
           write(51,50)  AccVolume,Npart,P(iden),&
           &          PLRC, EnergyAverage,&
           &          SEinputuMolPerM2,PressurePa, AccVolumeM3,EnergyAverageJoule,&
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
           &      1x, 'Ns_Density (uMol/m2) ------------------------ : ', e20.10,/, &
           &      1x, 'Pressure (Pa) ------------------------------ : ', e20.10,/, &
           &      1x, 'AccessVolume (m3) -------------------------- : ', e20.10,/, &
           &      1x, 'Average Energy per Particle (J) ------------ : ', e20.10,/, &
           &      1x, 'Time (sec) --------------------------------- : ', f20.10/)
!            -------------------------------
!            STORE THE POSITION OF PARTICLES
!            -------------------------------
             write(62,*) 'Npart=', Npart
             do i = 1, Npart
                write(52,436) mcx(i), mcy(i), mcz(i)
436             format(3f16.7)
             enddo  
             write(52,*) 'LJSITES' 
             do i = 1, Npart
                do j = 1, NLJ
                   write(52,436) LJX(i,j),LJY(i,j),LJZ(i,j)
                enddo
              enddo
             write(52,*) 'COULOMBSITES' 
             do i = 1, Npart
                do j = 1, NCL
                   write(52,436) CLX(i,j),CLY(i,j),CLZ(i,j)
                enddo
             enddo
!         ----------------------------------
!         STORE THE INFORMATION FOR ISOTHERM
!         ----------------------------------
           write(57,716)iden,TotPickN,OutZPickN,rCutOff, rhoMolPerM3,PressurePa, Npart, SEinputuMolPerM2,Npart-CNBulk, EnergyAverageJoule-CEBulkJoule,ChemPkVLEJPerMole,EntropySys,EntropySys/Rg, (ave_NinSBox(i), i=1,SubBox), (ChemPSBox(i), i=1,SubBox)
716        FORMAT(1x,I4,2I12,3E18.8,I8,200E18.8)  
           
flush(51)
flush(52)
flush(57)      

        ENDDO
        
        IF(DistanceD)call distancedistribution(4)
        DEALLOCATE(rate, AccuR, Energyi)
        DEALLOCATE(EnergyMatrix,VirialMatrix)
        IF(MassTSFER)DEALLOCATE(FluxSurf,UpFlux, DownFlux,NetFlux)
        

     RETURN
     END SUBROUTINE MC_SURFACE
 
     

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ENERGYSINGLEPARTICLE
!      THIS SUBROUTINE CALCULATES THE ENERGY OF ONE PARTICLE WITH THE OTHER 
!      PARTICLES AND THE SOLID
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine energysingleparticle(i, Energy, CLEnergy)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M 
       implicit none
       
       integer i, j,k,m

       real*8  Energy, CLEnergy
       real*8  d2, dc2,cld
       real*8  U, V, CLU, CLV
       real*8  xlj, ylj, zlj, xclj, yclj, zclj,xmclj, ymclj, zmclj
       real*8  xcl, ycl, zcl
       real*8  sigma,welldepth
       real*8  xmcn,ymcn,zmcn,mcd2,smcy
       real*8  xcmcn,ycmcn,zcmcn,cmcd2
       real*8  surfljd2, U1st, surfcld2, CLU1st
       real*8  PairE, PairV, pairDerivE
       
!      --------------------------------
!      INITIALISE THE ENERGY 
!      --------------------------------
       Energy   = 0.0E0
       CLEnergy = 0.0E0
       Virial   = 0.0E0
       VirialF  = 0.0E0
       EnergyF  = 0.0E0
       DerivE   = 0.0E0
       DerivE_Ale = 0.0E0
!      ------------------------------------
!      INTERACTION ENERGY BETWEEN PARTICLES
!      ------------------------------------
       DO j = 1, Npart !! 1-do
           IF(i .ne. j)THEN !! 5-if       
              Energy = Energy + EnergyMatrix(i,j)
              Virial = Virial + VirialMatrix(i,j)
              DerivE = DerivE + DerivEMatrix(i,j)
              DerivE_Ale = DerivE_Ale + DerivEMatrix_Ale(i,j)
              if(indNinPropSub(j).GT.0)then
                 VirialF(indNinPropSub(j)) = VirialF(indNinPropSub(j)) + VirialMatrix(i,j)
                 EnergyF(indNinPropSub(j)) = EnergyF(indNinPropSub(j)) + EnergyMatrix(i,j)
               endif
           ENDIF  !! 5-if
         ENDDO  !! 1-do
         
       RETURN
       ENDSUBROUTINE energysingleparticle              
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ENERGYSINGLEPARTICLE
!      THIS SUBROUTINE CALCULATES THE ENERGY OF ONE PARTICLE WITH THE OTHER 
!      PARTICLES AND THE SOLID
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine energysingleparticleADS(i, Energy, CLEnergy)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M 
       implicit none
       
       integer i, j,k,m

       real*8  Energy, CLEnergy
       real*8  d2, dc2,cld
       real*8  U, V, CLU, CLV
       real*8  xlj, ylj, zlj, xclj, yclj, zclj,xmclj, ymclj, zmclj
       real*8  xcl, ycl, zcl
       real*8  sigma,welldepth
       real*8  xmcn,ymcn,zmcn,mcd2,smcy
       real*8  xcmcn,ycmcn,zcmcn,cmcd2
       real*8  surfljd2, U1st, surfcld2, CLU1st
       real*8  PairE, PairV, pairDerivE
external omp_set_num_threads

     if((npart .gt. 50) .and. (npart .le. 150)) then
         call omp_set_num_threads(2)
     endif
    if((npart .gt. 150) .and. (npart .le. 250)) then
            call omp_set_num_threads(3)
    endif
    if((npart .gt. 250) .and. (npart .le. 350)) then
           call omp_set_num_threads(4)
     endif
    if((npart .gt. 350) .and. (npart .le. 450)) then
          call omp_set_num_threads(5)
     endif
   if(npart .gt. 450) then
           call omp_set_num_threads(6)
     endif
       
!      --------------------------------
!      INITIALISE THE ENERGY 
!      --------------------------------
       Energy   = 0.0E0
       CLEnergy = 0.0E0
       Virial   = 0.0E0
       BulkVirialF  = 0.0E0
       EinBulkF   = 0.0E0
!      ------------------------------------
!      INTERACTION ENERGY BETWEEN PARTICLES
!      ------------------------------------
!$OMP  parallel  if(npart .gt. 50)                                  &
!$OMP& shared(i,energymatrix,indninsbox,npart,subfgas,virialmatrix) &
!$OMP& private(j,energy, virial, bulkvirialf,einbulkf)

!$OMP do schedule(static)

!!!takes more time--> reduction(+:energy, virial, bulkvirialf,einbulkf)
       DO j = 1, Npart !! 1-do
           IF(i .ne. j)THEN !! 5-if       
              Energy = Energy + EnergyMatrix(i,j)
              Virial = Virial + VirialMatrix(i,j)
              if(indNinSBox(j).EQ.SubFGas)then
                 BulkVirialF = BulkVirialF + VirialMatrix(i,j)
                 EinBulkF    = EinBulkF + EnergyMatrix(i,j)
               endif
           ENDIF  !! 5-if
         ENDDO  !! 1-do
!$OMP end do nowait
!$OMP end parallel         
       RETURN
       ENDSUBROUTINE energysingleparticleADS                    


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE PairEnergy 
!      THIS SUBROUTINE CALCULATES THE ENERGY BETWEEN ONE SELECTED PARTICLE WITH
!      EVERY OTHER PARTICLE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine PairEnergy(i)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M 
       implicit none
       
       integer i, j,k,m, FLAG, MiniL

       real*8  Energy, EnergyC, CLEnergy
       real*8  d2, dc2,cld, distanceZ
       real*8  U, V, CLU, CLV
       real*8  xlj, ylj, zlj, xclj, yclj, zclj,xmclj, ymclj, zmclj
       real*8  xcl, ycl, zcl
       real*8  sigma,welldepth
       real*8  xmcn,ymcn,zmcn,mcd2
       real*8  surfljd2, U1st, surfcld2, CLU1st
       real*8  PairEij, PairVij, PairDerivEij, PairDerivEij_Ale
       real*8  MCVectorX, MCVectorY, MCVectorZ, MCVectorX0, MCVectorY0, MCVectorZ0
       real*8  LJVectorX, LJVectorY, LJVectorZ, LJVectorX0, LJVectorY0, LJVectorZ0
       real*8  CLVectorX, CLVectorY, CLVectorZ, CLVectorX0, CLVectorY0, CLVectorZ0
       real*8  VectorProduct
       real*8  U1st_Ale, CLU1st_Ale
         

!      ------------------------
!      INITIALISE THE ENERGY 
!      ------------------------
       PairEnergyi   = 0.0E0
!      ------------------------------------
!      INTERACTION ENERGY BETWEEN PARTICLES
!      ------------------------------------
       DO j = 1, Npart !! 1-do
          IF(i .ne. j)THEN !! 5-if
               MiniL = 0
               Energyi(i) = Energyi(i) - EnergyMatrix(i,j) ! Minus the old pariwise Energy 
               Energyi(j) = Energyi(j) - EnergyMatrix(i,j) ! Minus the old pariwise Energy
               PairEij = 0.0e0
               PairVij = 0.0e0
               PairDerivEij = 0.0e0 
               PairDerivEij_Ale = 0.0e0
               EnergyMatrix(i,j)=0.0e0
               EnergyMatrix(j,i)=0.0e0 
               VirialMatrix(i,j)=0.0e0
               VirialMatrix(j,i)=0.0e0
               DerivEMatrix(i,j)=0.0e0
               DerivEMatrix(j,i)=0.0e0
               DerivEMatrix_Ale(i,j)=0.0e0
               DerivEMatrix_Ale(j,i)=0.0e0
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
                     if(xmclj > VLEBoxlengthX/2.0D0)then
                        xmcn = xmclj - VLEBoxlengthX
                        IF(MCVectorX0 .GT. 0.0)MCVectorX = MCVectorX0 - VLEBoxlengthX
                        IF(MCVectorX0 .LT. 0.0)MCVectorX = MCVectorX0 + VLEBoxlengthX
                     endif
                endif
                if(pbcy)then
                     if(ymclj > VLEBoxlengthY/2.0D0)then
                        ymcn = ymclj - VLEBoxlengthY
                        IF(MCVectorY0 .GT. 0.0)MCVectorY = MCVectorY0 - VLEBoxlengthY
                        IF(MCVectorY0 .LT. 0.0)MCVectorY = MCVectorY0 + VLEBoxlengthY
                     endif
                endif
                if(pbcz)then
                     if(zmclj > VLEBoxlengthZ/2.0D0)then
                        zmcn = zmclj - VLEBoxlengthZ
                        IF(MCVectorZ0 .GT. 0.0)MCVectorZ = MCVectorZ0 - VLEBoxlengthZ
                        IF(MCVectorZ0 .LT. 0.0)MCVectorZ = MCVectorZ0 + VLEBoxlengthZ
                     endif
                endif
                mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
                surfljd2 = xmcn*xmcn + zmcn*zmcn - 2.0*ymcn*ymcn
                if(mcd2.gt.r2Cutoff)then !! 4-if
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
                            if(xmclj > VLEBoxlengthX/2.0D0)then
                               xlj = xlj - VLEBoxlengthX
                               IF(LJVectorX0.GT.0.0)LJVectorX = LJVectorX0 - VLEBoxlengthX
                               IF(LJVectorX0.LT.0.0)LJVectorX = LJVectorX0 + VLEBoxlengthX
                            endif
                      endif
                      if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0)then
                               ylj = ylj - VLEBoxlengthY
                               IF(LJVectorY0.GT.0.0)LJVectorY = LJVectorY0 - VLEBoxlengthY
                               IF(LJVectorY0.LT.0.0)LJVectorY = LJVectorY0 + VLEBoxlengthY
                            endif
                      endif
                      if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0)then
                               zlj = zlj - VLEBoxlengthZ
                               IF(LJVectorZ0.GT.0.0)LJVectorZ = LJVectorZ0 - VLEBoxlengthZ
                               IF(LJVectorZ0.LT.0.0)LJVectorZ = LJVectorZ0 + VLEBoxlengthZ
                            endif
                      endif
                
!                  -----------------------
!                  INTER-PARTICLE DISTANCE
!                  -----------------------
                      d2 = xlj*xlj + ylj*ylj + zlj*zlj
                      IF(sqrt(d2).LT.Dist_Limit*sigma)THEN
                         MiniL = 1
                         PairEij = PairEff_Limit*Temperature
                         PairVij = PairV_Limit
                         PairDerivEij = DerivE_Limit
                         PairDerivEij_Ale = DerivE_Limit
                         EXIT
                      ELSE
!                  ----------------
!                  POTENTIAL ENERGY
!                  ----------------
                          call potentialEnergy(d2, sigma, welldepth,U, V, U1st)

                          VectorProduct = MCVectorX*LJVectorX + MCVectorY*LJVectorY + MCVectorZ*LJVectorZ
                          V = V*VectorProduct/d2
                          U1st_Ale = U1st*(VectorProduct - 3.0*abs(ymcn)*abs(ylj))/2.0/sqrt(d2)
                          U1st = U1st*VectorProduct/sqrt(d2)
                          PairEij  = PairEij + U
                          PairVij  = PairVij + V
                          PairDerivEij = PairDerivEij + surfljd2/2.0/mcd2*U1st
                          PairDerivEij_Ale = PairDerivEij_Ale + U1st_Ale
                      ENDIF
                   enddo ! do m
                IF(MiniL.EQ.1)EXIT
                enddo ! do k

           !!------------------
           !! Coulombs force
           !!------------------
               IF(MiniL.EQ.0)THEN
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
                               if(xmclj > VLEBoxlengthX/2.0D0)then
                                  xcl = xcl - VLEBoxlengthX
                                  IF(CLVectorX0 .GT. 0.0)CLVectorX = CLVectorX0 - VLEBoxlengthX
                                  IF(CLVectorX0 .LT. 0.0)CLVectorX = CLVectorX0 + VLEBoxlengthX
                               endif
                         endif
                         if(pbcy)then
                               if(ymclj > VLEBoxlengthY/2.0D0)then
                                  ycl = ycl - VLEBoxlengthY
                                  IF(CLVectorY0.GT.0.0)CLVectorY = CLVectorY0 - VLEBoxlengthY
                                  IF(CLVectorY0.LT.0.0)CLVectorY = CLVectorY0 + VLEBoxlengthY
                               endif
                         endif
                         if(pbcz)then
                               if(zmclj > VLEBoxlengthZ/2.0D0)then
                                  zcl = zcl - VLEBoxlengthZ
                                  IF(CLVectorZ0.GT.0.0)CLVectorZ = CLVectorZ0 - VLEBoxlengthZ
                                  IF(CLVectorZ0.LT.0.0)CLVectorZ = CLVectorZ0 + VLEBoxlengthZ
                               endif
                         endif
                         cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                         IF(cld.LT.Dist_Limit)THEN
                            MiniL = 1
                            PairEij = PairEff_Limit*Temperature
                            PairVij = PairV_Limit
                            PairDerivEij = DerivE_Limit
                            PairDerivEij_Ale = DerivE_Limit
                            EXIT
                         ELSE
                            call coulombforce(cld, charges(k),charges(m),CLU,CLV, CLU1st)
                            VectorProduct = MCVectorX*CLVectorX + MCVectorY*CLVectorY + MCVectorZ*CLVectorZ
                            CLV = CLV*VectorProduct/cld**2.0
                            CLU1st_Ale = CLU1st*(VectorProduct - 3.0*abs(ymcn)*abs(ycl))/2.0/cld
                            CLU1st = CLU1st*VectorProduct/cld
                            PairEij  = PairEij + CLU
                            PairVij  = PairVij + CLV
                            PairDerivEij = PairDerivEij + surfljd2/2.0/mcd2*CLU1st
                            PairDerivEij_Ale = PairDerivEij_Ale + CLU1st
                         ENDIF
                       enddo ! do m
                    IF(MiniL.EQ.1)EXIT
                    enddo ! do k 
                ENDIF 
                    
                    IF(PairEij/Temperature.GT.PairEff_Limit)THEN
                        PairEij = PairEff_Limit*Temperature
                        PairVij = PairV_Limit
                        PairDerivEij = DerivE_Limit
                        PairDerivEij_Ale = DerivE_Limit
                    ENDIF
                    
                    EnergyMatrix(i,j) = PairEij
                    EnergyMatrix(j,i) = PairEij
                    VirialMatrix(i,j) = PairVij
                    VirialMatrix(j,i) = PairVij
                    DerivEMatrix(i,j) = PairDerivEij
                    DerivEMatrix(j,i) = PairDerivEij
                    DerivEMatrix_Ale(i,j) = PairDerivEij_Ale
                    DerivEMatrix_Ale(j,i) = PairDerivEij_Ale
              endif !! 4-if
              Energyi(i) = Energyi(i) + EnergyMatrix(i,j) ! Add the new pariwise Energy 
              Energyi(j) = Energyi(j) + EnergyMatrix(i,j) ! Add the new pariwise Energy
         ENDIF !! 5-if
      ENDDO  !! 1-do
       
         TotRate = 0.0
         DO j = 1, Npart
            rate(j) = exp(Energyi(j)/Temperature)
            TotRate = TotRate + rate(j)
         ENDDO
 
             
     RETURN
     ENDSUBROUTINE PairEnergy        

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE PairEnergyADS 
!      THIS SUBROUTINE CALCULATES THE ENERGY BETWEEN ONE SELECTED PARTICLE WITH
!      EVERY OTHER PARTICLE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine PairEnergyADS(i)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       USE Constant_M
       implicit none
       
       double precision:: ostart, oend 
       integer i, j,k,m, FLAG, MiniL
       integer FLAG1, FLAG2

       real*8  d6, d2, dc2,cld, distanceZ
       real*8  U, V, CLU, CLV, UC, BU, PEnergyC
       real*8  xlj, ylj, zlj, xclj, yclj, zclj,xmclj, ymclj, zmclj
       real*8  xcl, ycl, zcl
       real*8  sigma,welldepth
       real*8  xmcn,ymcn,zmcn,mcd2
       real*8  U1st, surfcld2, CLU1st
       real*8  PairEij, PairVij
       real*8  MCVectorX, MCVectorY, MCVectorZ, MCVectorX0, MCVectorY0, MCVectorZ0
       real*8  LJVectorX, LJVectorY, LJVectorZ, LJVectorX0, LJVectorY0, LJVectorZ0
       real*8  CLVectorX, CLVectorY, CLVectorZ, CLVectorX0, CLVectorY0, CLVectorZ0
       real*8  VectorProduct    
       
external omp_set_num_threads

!     if((npart .gt. 50) .and. (npart .le. 150)) then
!         call omp_set_num_threads(2)
!     endif
!     if((npart .gt. 150) .and. (npart .le. 250)) then
!           call omp_set_num_threads(3)
!    endif
!     if((npart .gt. 250) .and. (npart .le. 350)) then
!           call omp_set_num_threads(4)
!     endif
!     if((npart .gt. 350) .and. (npart .le. 450)) then
!          call omp_set_num_threads(5)
!     endif
!    if(npart .gt. 450) then
!           call omp_set_num_threads(6)
!     endif

!      ------------------------------------
!      INTERACTION ENERGY BETWEEN PARTICLES
!      ------------------------------------
       FLAG1=indNinSBox(i)
! Directive inserted by Cray Reveal.  May be incomplete.
!!$OMP  parallel if(npart .gt. 50)                                        &
!!$OMP&   private (cld,clu,clu1st,clv,clvectorx,clvectorx0,clvectory,     &
!!$OMP&            clvectory0,clvectorz,clvectorz0,d2,d6, flag2,j,         &
!!$OMP&            ljvectorx,ljvectorx0,ljvectory,ljvectory0,ljvectorz,   &
!!$OMP&            ljvectorz0,mcd2,mcvectorx,mcvectorx0,mcvectory,      &
!!$OMP&            mcvectory0,mcvectorz,mcvectorz0,minil,sigma, &
!!$OMP&            u,u1st,vectorproduct,welldepth,xmclj,xmcn,   &
!!$OMP&            ymclj,ymcn,zmclj,zmcn)& 
!!!!!!!paireij, pairvij, xlj,ylj,zlj,xcl,ycl,zcl,V,minil,k,m)     &
!!$OMP&   shared  (i,boxlengthx,boxlengthy,boxlengthz,charges,clx,cly,    &
!!$OMP&            clz,dist_limit,extrabox,flag1,indninsbox,ljx,ljy,ljz,  &
!!$OMP&            mcx,mcy,mcz,ncl,nlj,npart,paireff_limit,pairv_limit,   &
!!$OMP&            pbcx,pbcy,pbcz,pore,r2cutoff,sigmaff,subbox,subfgas,   &
!!$OMP&            temperature,welldepthff,energyi,energymatrix,virialmatrix)
!
!!$OMP do schedule(static) reduction(+:paireij, pairvij) reduction(-:xlj,ylj,zlj,xcl,ycl,zcl) reduction(*:V)
       DO j = 1, Npart !! DO j
          IF(i .ne. j)THEN !! 5-if
               MiniL = 0
               FLAG2=indNinSBox(j)
!!$OMP critical
               Energyi(i) = Energyi(i) - EnergyMatrix(i,j)   ! Minus the old pariwise Energy 
               Energyi(j) = Energyi(j) - EnergyMatrix(i,j) ! Minus the old pariwise Energy
!!$OMP end critical
               PairEij = 0.0e0
               PairVij = 0.0e0
               EnergyMatrix(i,j)=0.0e0
               EnergyMatrix(j,i)=0.0e0 
               VirialMatrix(i,j)=0.0e0
               VirialMatrix(j,i)=0.0e0

               IF(PORE.AND.ExtraBox)THEN
                  IF((FLAG1.EQ.SubBOX) .AND. (FLAG1.NE.FLAG2))CYCLE
                  IF((FLAG2.EQ.SubBOX) .AND. (FLAG1.NE.FLAG2))CYCLE
               ENDIF

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
                     if(xmclj > BoxlengthX/2.0D0)then
                        xmcn = xmclj - BoxlengthX
                        IF(MCVectorX0 .GT. 0.0)MCVectorX = MCVectorX0 - BoxlengthX
                        IF(MCVectorX0 .LT. 0.0)MCVectorX = MCVectorX0 + BoxlengthX
                     endif
                endif
                if(pbcy)then
                     if(ymclj > BoxlengthY/2.0D0)then
                        ymcn = ymclj - BoxlengthY
                        IF(MCVectorY0 .GT. 0.0)MCVectorY = MCVectorY0 - BoxlengthY
                        IF(MCVectorY0 .LT. 0.0)MCVectorY = MCVectorY0 + BoxlengthY
                     endif
                endif
                if(pbcz .OR. (PORE .AND. (FLAG1.EQ.SubFGas) .AND.(FLAG2.EQ.SubFGas)) )then
                     if(zmclj > BoxlengthZ/2.0D0)then
                        zmcn = zmclj - BoxlengthZ
                        IF(MCVectorZ0 .GT. 0.0)MCVectorZ = MCVectorZ0 - BoxlengthZ
                        IF(MCVectorZ0 .LT. 0.0)MCVectorZ = MCVectorZ0 + BoxlengthZ
                     endif
                endif
                mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
                if(mcd2.gt.r2Cutoff)then !! 4-if
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
                            if(xmclj > BoxlengthX/2.0D0)then
                               xlj = xlj - BoxlengthX
                               IF(LJVectorX0.GT.0.0)LJVectorX = LJVectorX0 - BoxlengthX
                               IF(LJVectorX0.LT.0.0)LJVectorX = LJVectorX0 + BoxlengthX
                            endif
                      endif
                      if(pbcy)then
                            if(ymclj > BoxlengthY/2.0D0)then
                               ylj = ylj - BoxlengthY
                               IF(LJVectorY0.GT.0.0)LJVectorY = LJVectorY0 - BoxlengthY
                               IF(LJVectorY0.LT.0.0)LJVectorY = LJVectorY0 + BoxlengthY
                            endif
                      endif
                      if(pbcz .OR. (PORE .AND. (FLAG1.EQ.SubFGas) .AND.(FLAG2.EQ.SubFGas)) )then
                            if(zmclj > BoxlengthZ/2.0D0)then
                               zlj = zlj - BoxlengthZ
                               IF(LJVectorZ0.GT.0.0)LJVectorZ = LJVectorZ0 - BoxlengthZ
                               IF(LJVectorZ0.LT.0.0)LJVectorZ = LJVectorZ0 + BoxlengthZ
                            endif
                      endif
                
!                  -----------------------
!                  INTER-PARTICLE DISTANCE
!                  -----------------------
                      d2 = xlj*xlj + ylj*ylj + zlj*zlj
                      IF(sqrt(d2).LT.Dist_Limit*sigma)THEN
                         MiniL = 1
                         PairEij = PairEff_Limit*Temperature
                         PairVij = PairV_Limit
                         EXIT
                      ELSE
!                  ----------------
!                  POTENTIAL ENERGY
!                  ----------------
                          !call potentialEnergy(d2, sigma, welldepth,U, V, U1st)
                            d6  = (sigma**2/d2)**3
           U   =  4.0E0*welldepth*d6*(d6 - 1.0E0)
           V   =  16.0E0*welldepth*d6*(d6 - 0.5E0)
           U1st = -24.0E0*welldepth*(sigma**6.0)*( 2.0*(sigma**6)/(d2**6.5) - 1.0/(d2**3.5) )
!---
           VectorProduct = MCVectorX*LJVectorX + MCVectorY*LJVectorY + MCVectorZ*LJVectorZ
                          V = V*VectorProduct/d2
                          PairEij  = PairEij + U
                          PairVij  = PairVij + V
                      ENDIF
                   enddo ! do m
                IF(MiniL.EQ.1)EXIT
                enddo ! do k

           !!------------------
           !! Coulombs force
           !!------------------
               IF(MiniL.EQ.0)THEN
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
                               if(xmclj > BoxlengthX/2.0D0)then
                                  xcl = xcl - BoxlengthX
                                  IF(CLVectorX0 .GT. 0.0)CLVectorX = CLVectorX0 - BoxlengthX
                                  IF(CLVectorX0 .LT. 0.0)CLVectorX = CLVectorX0 + BoxlengthX
                               endif
                         endif
                         if(pbcy)then
                               if(ymclj > BoxlengthY/2.0D0)then
                                  ycl = ycl - BoxlengthY
                                  IF(CLVectorY0.GT.0.0)CLVectorY = CLVectorY0 - BoxlengthY
                                  IF(CLVectorY0.LT.0.0)CLVectorY = CLVectorY0 + BoxlengthY
                               endif
                         endif
                         if(pbcz .OR. (PORE .AND. (FLAG1.EQ.SubFGas) .AND.(FLAG2.EQ.SubFGas)) )then
                               if(zmclj > BoxlengthZ/2.0D0)then
                                  zcl = zcl - BoxlengthZ
                                  IF(CLVectorZ0.GT.0.0)CLVectorZ = CLVectorZ0 - BoxlengthZ
                                  IF(CLVectorZ0.LT.0.0)CLVectorZ = CLVectorZ0 + BoxlengthZ
                               endif
                         endif
                         cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                         IF(cld.LT.Dist_Limit)THEN
                            MiniL = 1
                            PairEij = PairEff_Limit*Temperature
                            PairVij = PairV_Limit
                            EXIT
                         ELSE
                            !call coulombforce(cld, charges(k),charges(m),CLU,CLV, CLU1st)
        CLU = charges(k)*charges(m)/(4.0e0*pi*permittivity*Scalelength*1.0e-10*scaleEnergy*kB*cld)
        CLV = CLU/3.0
        CLU1st = -CLU/cld
!--
                            VectorProduct = MCVectorX*CLVectorX + MCVectorY*CLVectorY + MCVectorZ*CLVectorZ
                            CLV = CLV*VectorProduct/cld**2.0
                            PairEij  = PairEij + CLU
                            PairVij  = PairVij + CLV
                         ENDIF
                       enddo ! do m
                    IF(MiniL.EQ.1)EXIT
                    enddo ! do k 
                ENDIF 
                    
                    IF(PairEij/Temperature.GT.PairEff_Limit)THEN
                        PairEij = PairEff_Limit*Temperature
                        PairVij = PairV_Limit
                    ENDIF
                    
                    EnergyMatrix(i,j) = PairEij
                    EnergyMatrix(j,i) = PairEij
                    VirialMatrix(i,j) = PairVij
                    VirialMatrix(j,i) = PairVij
              endif !! 4-if
!!$OMP critical
             Energyi(i) = Energyi(i) + EnergyMatrix(i,j) ! Add the new pariwise Energy 
              Energyi(j) = Energyi(j) + EnergyMatrix(i,j) ! Add the new pariwise Energy
!!$OMP end critical
         ENDIF !! 5-if
      ENDDO  !! DO j
!!$OMP end do nowait
!!$OMP end parallel       
!         --------------------------------------
!         INTERACTION BETWEEN PARTICLE AND SOLID
!         --------------------------------------         
          Energyi(i) = Energyi(i) - EnergyMatrix(i,i)
          IF( SURFACE .OR. (PORE .AND. (MCy(i).GE.AccVYL) .AND. (MCy(i).LE.AccVYH)) )THEN
             if(steele)then
                  IF(PORE.AND. ExtraBox)THEN
                       IF(FLAG1.eq.SubFGas)THEN
                          UC=0.0E0
                       ELSE
                          call steelepotential(i,UC)
                       ENDIF
                   ELSE
                       call steelepotential(i,UC)
                   ENDIF
                   EnergyMatrix(i,i) = UC
              endif

              if(Bojan)then
                    PEnergyC = 0.0
                    do j = 1, NS
                      call BojanPotential(j, i,BU)
                      PEnergyC = PEnergyC + BU
                    enddo
                    IF(Close1End.OR.Close2Ends)then
                       call BojanPotential(NS+2, i,BU) 
                       PEnergyC = PEnergyC + BU 
                    ENDIF
                    IF(PEnergyC/Temperature.GT.PairEsf_Limit)PEnergyC = PairEsf_Limit*Temperature
                    EnergyMatrix(i,i) = PEnergyC
              endif
           ELSE
              EnergyMatrix(i,i) = 0.0E0
          ENDIF
          Energyi(i) = Energyi(i) + EnergyMatrix(i,i)
        
         TotRate = 0.0
         rateSBox = 0.0
         DO j = 1, Npart
            rate(j) = exp(Energyi(j)/Temperature)
            TotRate = TotRate + rate(j)
            rateSBox(indNinSBox(j)) = rateSBox(indNinSBox(j)) + rate(j)
         ENDDO
                      
     RETURN
     ENDSUBROUTINE PairEnergyADS        

   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE CONFIGENERGY
!      THIS SUBROUTINE CALCULATES THE CONFIGURATION ENERGY OF N PARTICLES 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ConfigEnergy(Energy,CLEnergy)
      USE FLUID_SOLID_M
      USE PHYSICAL_M
      USE MCSETTING_M
      USE SUBBOX_M
      implicit none
 
!     ------------------------
!     DECLARATION OF VARIABLES
!     ------------------------
      integer i, j, N, MiniL
      integer k, m, FLAG1, FLAG2
      
      real*8 sigma, welldepth               
      real*8 Energy,CLEnergy,EnergyC 
      real*8 U,V,CLU, CLV
      real*8 xlj, ylj, zlj, d2
      real*8 xcl,ycl,zcl,cld
      real*8 xclj, yclj, zclj, dc2
      real*8 xmclj,ymclj,zmclj
      real*8 xmcn,ymcn,zmcn,mcd2
      real*8 xcmcn,ycmcn,zcmcn,cmcd2
      real*8 surfljd2, U1st, surfcld2,CLU1st
      real*8 PairE, PairV,  pairFirstDerivE, pairFirstDerivE_Ale
      real*8 MCVectorX, MCVectorY, MCVectorZ, MCVectorX0, MCVectorY0, MCVectorZ0
      real*8 LJVectorX, LJVectorY, LJVectorZ, LJVectorX0, LJVectorY0, LJVectorZ0
      real*8 CLVectorX, CLVectorY, CLVectorZ, CLVectorX0, CLVectorY0, CLVectorZ0
      real*8 VectorProduct
      real*8 U1st_Ale, CLU1st_Ale

!     -----------------------------------
!     INITIALIZATION OF ENERGY 
!     -----------------------------------
       Energy        = 0.0E0
       CLEnergy      = 0.0E0
       Virial        = 0.0E0
       VirialCL      = 0.0E0
       N             = Npart  
       VirialProp    = 0.0
       EnergyProp    = 0.0
       FirstDerivE   = 0.0
       FirstDerivE_Ale = 0.0
       TotRate       = 0.0
       rate          = 0.0
       AccuR         = 0.0
       Energyi       = 0.0
       
       EnergyMatrix = 0.0E0
       VirialMatrix = 0.0E0
       DerivEMatrix = 0.0E0
       DerivEMatrix_Ale = 0.0E0
       
       if(N.eq.0) return
         
!     ---------------------------------------------------------------------
!     SUMMING THE PAIRWISE INTERACTION ENERGIES FOR PROPERTIES CALCULATION
!     ---------------------------------------------------------------------
        DO i=1, N-1
           FLAG1=indNinPropSub(i)
           DO j=i+1, N
              MiniL = 0
              FLAG2=indNinPropSub(j)
              PairE = 0.0
              PairV = 0.0
              pairFirstDerivE = 0.0e0
              pairFirstDerivE_Ale = 0.0e0
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
                       if(xmclj > VLEBoxlengthX/2.0D0)then
                          xmcn = xmclj - VLEBoxlengthX
                          IF(MCVectorX0.GT.0.0)MCVectorX = MCVectorX0 - VLEBoxlengthX
                          IF(MCVectorX0.LT.0.0)MCVectorX = MCVectorX0 + VLEBoxlengthX
                       endif
                 endif
                 if(pbcy)then
                       if(ymclj > VLEBoxlengthY/2.0D0)then
                           ymcn = ymclj - VLEBoxlengthY
                           IF(MCVectorY0.GT.0.0)MCVectorY = MCVectorY0 - VLEBoxlengthY
                           IF(MCVectorY0.LT.0.0)MCVectorY = MCVectorY0 + VLEBoxlengthY
                       endif
                 endif
                 if(pbcz)then
                      if(zmclj > VLEBoxlengthZ/2.0D0)then
                           zmcn = zmclj - VLEBoxlengthZ
                           IF(MCVectorZ0.GT.0.0)MCVectorZ = MCVectorZ0 - VLEBoxlengthZ
                           IF(MCVectorZ0.LT.0.0)MCVectorZ = MCVectorZ0 + VLEBoxlengthZ 
                      endif
                 endif 
                 mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
                 surfljd2 = xmcn*xmcn + zmcn*zmcn - 2.0*ymcn*ymcn
              if(mcd2.gt.r2cutoff)then
                 CYCLE
              else
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
                            if(xmclj > VLEBoxlengthX/2.0D0)then
                               xlj = xlj - VLEBoxlengthX
                               IF(LJVectorX0 .GT. 0.0)LJVectorX = LJVectorX0 - VLEBoxlengthX
                               IF(LJVectorX0 .LT. 0.0)LJVectorX = LJVectorX0 + VLEBoxlengthX
                            endif
                       endif
                       if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0)then
                               ylj = ylj - VLEBoxlengthY
                               IF(LJVectorY0 .GT. 0.0)LJVectorY = LJVectorY0 - VLEBoxlengthY
                               IF(LJVectorY0 .LT. 0.0)LJVectorY = LJVectorY0 + VLEBoxlengthY
                            endif 
                       endif
                       if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0)then
                                zlj = zlj - VLEBoxlengthZ
                                IF(LJVectorZ0.GT.0.0)LJVectorZ = LJVectorZ0 - VLEBoxlengthZ
                                IF(LJVectorZ0.LT.0.0)LJVectorZ = LJVectorZ0 + VLEBoxlengthZ
                            endif
                       endif 
!                -----------------------
!                INTER-PARTICLE DISTANCE
!                -----------------------
                       d2 = xlj*xlj + ylj*ylj + zlj*zlj
                       IF(sqrt(d2) .LT. Dist_Limit*sigma)THEN
                          MiniL = 1
                          PairE = PairEff_Limit*Temperature
                          PairV = PairV_Limit
                          pairFirstDerivE = DerivE_Limit
                          pairFirstDerivE_Ale = DerivE_Limit
                          EXIT
                       ELSE
!                ----------------
!                POTENTIAL ENERGY
!                ----------------
                          call potentialEnergy(d2,sigma,welldepth,U,V,U1st)
                           VectorProduct = MCVectorX*LJVectorX + MCVectorY*LJVectorY + MCVectorZ*LJVectorZ
                           V = V*VectorProduct/d2
                           U1st_Ale = U1st*(VectorProduct-3.0*abs(ymcn)*abs(ylj))/2.0/sqrt(d2)
                           U1st = U1st*VectorProduct/sqrt(d2)
                           PairE = PairE + U
                           PairV = PairV + V
                           pairFirstDerivE = pairFirstDerivE + surfljd2/2.0/mcd2*U1st
                           pairFirstDerivE_Ale = pairFirstDerivE_Ale + U1st_Ale
                       ENDIF
                    enddo ! do m
                  IF(MiniL.EQ.1)EXIT
                  enddo  ! do k
             !!---------------------------
             !! COULOMB FORCE
             !!--------------------------- 
              IF(MiniL.EQ.0)THEN
                  do k = 1, NCL
                     do m = 1, NCL
                      sigma = (sigmaFF(k)+sigmaFF(m))/2.0e0
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
                         if(xmclj > VLEBoxlengthX/2.0D0)then
                             xcl = xcl - VLEBoxlengthX
                             IF(CLVectorX0.GT.0.0)CLVectorX = CLVectorX0 - VLEBoxlengthX
                             IF(CLVectorX0.LT.0.0)CLVectorX = CLVectorX0 + VLEBoxlengthX
                          endif
                      endif
                      if(pbcy)then
                          if(ymclj > VLEBoxlengthY/2.0D0)then
                              ycl = ycl - VLEBoxlengthY
                              IF(CLVectorY0.GT.0.0)CLVectorY = CLVectorY0 - VLEBoxlengthY
                              IF(CLVectorY0.LT.0.0)CLVectorY = CLVectorY0 + VLEBoxlengthY
                           endif
                      endif
                      if(pbcz)then
                          if(zmclj > VLEBoxlengthZ/2.0D0)then
                              zcl = zcl - VLEBoxlengthZ
                              IF(CLVectorZ0.GT.0.0)CLVectorZ = CLVectorZ0 - VLEBoxlengthZ
                              IF(CLVectorZ0.LT.0.0)CLVectorZ = CLVectorZ0 + VLEBoxlengthZ
                          endif 
                      endif
                      cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                        IF(cld .LT. Dist_Limit*sigma)THEN
                           MiniL = 1
                           PairE = PairEff_Limit*Temperature
                           PairV = PairV_Limit
                           pairFirstDerivE = DerivE_Limit
                           pairFirstDerivE_Ale = DerivE_Limit
                           EXIT
                        ELSE
                           call coulombforce(cld, charges(k),charges(m),CLU,CLV,CLU1st)
                           VectorProduct = MCVectorX*CLVectorX + MCVectorY*CLVectorY + MCVectorZ*CLVectorZ
                           CLV = CLV*VectorProduct/cld**2.0
                           CLU1st_Ale = CLU1st*(VectorProduct-3.0*abs(ymcn)*abs(ycl))/2.0/cld
                           CLU1st = CLU1st*VectorProduct/cld
                           CLEnergy = CLEnergy + CLU
                           VirialCL = VirialCL + CLV
                           PairE = PairE + CLU
                           PairV = PairV + CLV
                           pairFirstDerivE = pairFirstDerivE + surfljd2/2.0/mcd2*CLU1st 
                           pairFirstDerivE_Ale = pairFirstDerivE_Ale + CLU1st_Ale  
                        ENDIF    
                     enddo ! do m
                  IF(MiniL.EQ.1)EXIT
                  enddo ! do k
                ENDIF
                       IF(PairE/Temperature.GT.PairEff_Limit)THEN
                          PairE = PairEff_Limit*Temperature
                          PairV = PairV_Limit
                          pairFirstDerivE = DerivE_Limit
                          pairFirstDerivE_Ale = DerivE_Limit
                       ENDIF
                       
                       EnergyMatrix(i,j) = PairE
                       EnergyMatrix(j,i) = PairE  
                       VirialMatrix(i,j) = PairV
                       VirialMatrix(j,i) = PairV
                       DerivEMatrix(i,j) = pairFirstDerivE
                       DerivEMatrix(j,i) = pairFirstDerivE
                       DerivEMatrix_Ale(i,j) = pairFirstDerivE_Ale
                       DerivEMatrix_Ale(j,i) = pairFirstDerivE_Ale
                           
                       Energy = Energy + PairE
                       Virial = Virial + PairV
                       FirstDerivE = FirstDerivE + pairFirstDerivE 
                       FirstDerivE_Ale = FirstDerivE_Ale + pairFirstDerivE_Ale
                        
                        if((FLAG1.GT.0).AND.(FLAG1.EQ.FLAG2))then
                             VirialProp(FLAG1)= VirialProp(FLAG1)+PairV
                             EnergyProp(FLAG1)= EnergyProp(FLAG1)+PairE
                        endif
                        if((FLAG1.GT.0).AND.(FLAG2.EQ.0))then
                             VirialProp(FLAG1)= VirialProp(FLAG1)+PairV/2.0
                             EnergyProp(FLAG1)= EnergyProp(FLAG1)+PairE/2.0
                        endif
                        if((FLAG1.EQ.0).AND.(FLAG2.GT.0))then
                             VirialProp(FLAG2)= VirialProp(FLAG2)+PairV/2.0
                             EnergyProp(FLAG2)= EnergyProp(FLAG2)+PairE/2.0
                        endif
                        if((FLAG1.GT.0).AND.(FLAG2.GT.0).AND.(FLAG1.NE.FLAG2))then
                             VirialProp(FLAG1)= VirialProp(FLAG1)+PairV/2.0
                             VirialProp(FLAG2)= VirialProp(FLAG2)+PairV/2.0
                             EnergyProp(FLAG1)= EnergyProp(FLAG1)+PairE/2.0
                             EnergyProp(FLAG2)= EnergyProp(FLAG2)+PairE/2.0
                        endif         
              endif               
            ENDDO    ! DO j
         ENDDO  ! Do i
        
!     --------------------------------------------------------
!     CALCULATE THE ENERGY OF EACH PARTICLE WITH THE REST
!     --------------------------------------------------------
         DO i = 1, N
            do j = 1, N !! j-do
              IF(i .ne. j)THEN !! 5-if
                 Energyi(i) = Energyi(i) + EnergyMatrix(i,j)
              ENDIF
            enddo
            rate(i) = exp(Energyi(i)/Temperature)
            TotRate = TotRate + rate(i)
         ENDDO
                  
     RETURN
     END SUBROUTINE ConfigEnergy  
    
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE CONFIGENERGYADS
!      THIS SUBROUTINE CALCULATES THE CONFIGURATION ENERGY OF N PARTICLES 
!      AND WITH SOLID
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ConfigEnergyADS(Energy,CLEnergy,EnergyC)
      USE FLUID_SOLID_M
      USE PHYSICAL_M
      USE MCSETTING_M
      USE SUBBOX_M
      implicit none
 
!     ------------------------
!     DECLARATION OF VARIABLES
!     ------------------------
      integer i, j, N, MiniL
      integer k, m, FLAG1, FLAG2
      
      real*8 sigma, welldepth               
      real*8 Energy,CLEnergy, EnergyC, PEnergyC 
      real*8 U,V,CLU, CLV, UC, BU
      real*8 xlj, ylj, zlj, d2
      real*8 xcl,ycl,zcl,cld
      real*8 xclj, yclj, zclj, dc2
      real*8 xmclj,ymclj,zmclj
      real*8 xmcn,ymcn,zmcn,mcd2
      real*8 xcmcn,ycmcn,zcmcn,cmcd2
      real*8 U1st, CLU1st
      real*8 PairE, PairV
      real*8 MCVectorX, MCVectorY, MCVectorZ, MCVectorX0, MCVectorY0, MCVectorZ0
      real*8 LJVectorX, LJVectorY, LJVectorZ, LJVectorX0, LJVectorY0, LJVectorZ0
      real*8 CLVectorX, CLVectorY, CLVectorZ, CLVectorX0, CLVectorY0, CLVectorZ0
      real*8 VectorProduct
      real*8 U1st_Ale, CLU1st_Ale

!     -----------------------------------
!     INITIALIZATION OF ENERGY 
!     -----------------------------------
       N             = Npart  
       Energy        = 0.0E0
       CLEnergy      = 0.0E0
       EnergyC       = 0.0E0
       EinBulk       = 0.0E0
       Virial        = 0.0E0
       BulkVirial    = 0.0E0
       TotRate       = 0.0E0
       rate          = 0.0E0
       AccuR         = 0.0E0
       Energyi       = 0.0E0
       rateSBox      = 0.0E0
       
       EnergyMatrix = 0.0E0
       VirialMatrix = 0.0E0
       
       if(N.eq.0) return
         
!     ---------------------------------------------------------------------
!     SUMMING THE PAIRWISE INTERACTION ENERGIES FOR PROPERTIES CALCULATION
!     ---------------------------------------------------------------------
        DO i=1, N-1
           FLAG1=indNinSBox(i)
           DO j=i+1, N
              MiniL = 0
              FLAG2=indNinSBox(j)
              PairE = 0.0
              PairV = 0.0

              IF(PORE.AND.ExtraBox)THEN
                 if((FLAG1.EQ.SubFGas).AND.(FLAG1.NE.FLAG2))CYCLE
                 if((FLAG2.EQ.SubFGas).AND.(FLAG1.NE.FLAG2))CYCLE
              ENDIF

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
                       if(xmclj > BoxlengthX/2.0D0)then
                          xmcn = xmclj - BoxlengthX
                          IF(MCVectorX0.GT.0.0)MCVectorX = MCVectorX0 - BoxlengthX
                          IF(MCVectorX0.LT.0.0)MCVectorX = MCVectorX0 + BoxlengthX
                       endif
                 endif
                 if(pbcy)then
                       if(ymclj > BoxlengthY/2.0D0)then
                           ymcn = ymclj - BoxlengthY
                           IF(MCVectorY0.GT.0.0)MCVectorY = MCVectorY0 - BoxlengthY
                           IF(MCVectorY0.LT.0.0)MCVectorY = MCVectorY0 + BoxlengthY
                       endif
                 endif
                 if(pbcz .OR. (PORE .AND. (FLAG1.EQ.SubFGas) .AND.(FLAG2.EQ.SubFGas)) )then
                      if(zmclj > BoxlengthZ/2.0D0)then
                           zmcn = zmclj - BoxlengthZ
                           IF(MCVectorZ0.GT.0.0)MCVectorZ = MCVectorZ0 - BoxlengthZ
                           IF(MCVectorZ0.LT.0.0)MCVectorZ = MCVectorZ0 + BoxlengthZ 
                      endif
                 endif 
                 mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
              if(mcd2.gt.r2cutoff)then
                 CYCLE
              else
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
                            if(xmclj > BoxlengthX/2.0D0)then
                               xlj = xlj - BoxlengthX
                               IF(LJVectorX0 .GT. 0.0)LJVectorX = LJVectorX0 - BoxlengthX
                               IF(LJVectorX0 .LT. 0.0)LJVectorX = LJVectorX0 + BoxlengthX
                            endif
                       endif
                       if(pbcy)then
                            if(ymclj > BoxlengthY/2.0D0)then
                               ylj = ylj - BoxlengthY
                               IF(LJVectorY0 .GT. 0.0)LJVectorY = LJVectorY0 - BoxlengthY
                               IF(LJVectorY0 .LT. 0.0)LJVectorY = LJVectorY0 + BoxlengthY
                            endif 
                       endif
                       if(pbcz .OR. (PORE .AND. (FLAG1.EQ.SubFGas) .AND.(FLAG2.EQ.SubFGas)) )then
                            if(zmclj > BoxlengthZ/2.0D0)then
                                zlj = zlj - BoxlengthZ
                                IF(LJVectorZ0.GT.0.0)LJVectorZ = LJVectorZ0 - BoxlengthZ
                                IF(LJVectorZ0.LT.0.0)LJVectorZ = LJVectorZ0 + BoxlengthZ
                            endif
                       endif 
!                -----------------------
!                INTER-PARTICLE DISTANCE
!                -----------------------
                       d2 = xlj*xlj + ylj*ylj + zlj*zlj
                       IF(sqrt(d2) .LT. Dist_Limit*sigma)THEN
                          MiniL = 1
                          PairE = PairEff_Limit*Temperature
                          PairV = PairV_Limit
                          EXIT
                       ELSE
!                ----------------
!                POTENTIAL ENERGY
!                ----------------
                          call potentialEnergy(d2,sigma,welldepth,U,V,U1st)
                           VectorProduct = MCVectorX*LJVectorX + MCVectorY*LJVectorY + MCVectorZ*LJVectorZ
                           V = V*VectorProduct/d2
                           PairE = PairE + U
                           PairV = PairV + V
                       ENDIF
                    enddo ! do m
                  IF(MiniL.EQ.1)EXIT
                  enddo  ! do k
             !!---------------------------
             !! COULOMB FORCE
             !!--------------------------- 
              IF(MiniL.EQ.0)THEN
                  do k = 1, NCL
                     do m = 1, NCL
                      sigma = (sigmaFF(k)+sigmaFF(m))/2.0e0
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
                         if(xmclj > BoxlengthX/2.0D0)then
                             xcl = xcl - BoxlengthX
                             IF(CLVectorX0.GT.0.0)CLVectorX = CLVectorX0 - BoxlengthX
                             IF(CLVectorX0.LT.0.0)CLVectorX = CLVectorX0 + BoxlengthX
                          endif
                      endif
                      if(pbcy)then
                          if(ymclj > BoxlengthY/2.0D0)then
                              ycl = ycl - BoxlengthY
                              IF(CLVectorY0.GT.0.0)CLVectorY = CLVectorY0 - BoxlengthY
                              IF(CLVectorY0.LT.0.0)CLVectorY = CLVectorY0 + BoxlengthY
                           endif
                      endif
                      if(pbcz.OR. (PORE .AND. (FLAG1.EQ.SubFGas) .AND.(FLAG2.EQ.SubFGas)) )then
                          if(zmclj > BoxlengthZ/2.0D0)then
                              zcl = zcl - BoxlengthZ
                              IF(CLVectorZ0.GT.0.0)CLVectorZ = CLVectorZ0 - BoxlengthZ
                              IF(CLVectorZ0.LT.0.0)CLVectorZ = CLVectorZ0 + BoxlengthZ
                          endif 
                      endif
                      cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                        IF(cld .LT. Dist_Limit*sigma)THEN
                           MiniL = 1
                           PairE = PairEff_Limit*Temperature
                           PairV = PairV_Limit
                           EXIT
                        ELSE
                           call coulombforce(cld, charges(k),charges(m),CLU,CLV,CLU1st)
                           VectorProduct = MCVectorX*CLVectorX + MCVectorY*CLVectorY + MCVectorZ*CLVectorZ
                           CLV = CLV*VectorProduct/cld**2.0
                           CLEnergy = CLEnergy + CLU
                           PairE = PairE + CLU
                           PairV = PairV + CLV
                        ENDIF    
                     enddo ! do m
                  IF(MiniL.EQ.1)EXIT
                  enddo ! do k
                ENDIF
                       IF(PairE/Temperature.GT.PairEff_Limit)THEN
                          PairE = PairEff_Limit*Temperature
                          PairV = PairV_Limit
                       ENDIF
                       
                       EnergyMatrix(i,j) = PairE
                       EnergyMatrix(j,i) = PairE  
                       VirialMatrix(i,j) = PairV
                       VirialMatrix(j,i) = PairV
                           
                       Energy = Energy + PairE
                       Virial = Virial + PairV
                        
                       
                           if((FLAG1.EQ.SubFGas).AND.(FLAG1.EQ.FLAG2))then
                                BulkVirial = BulkVirial + PairV
                                EinBulk = EinBulk + PairE
                           endif
                           if((FLAG1.EQ.SubFGas).AND.(FLAG1.NE.FLAG2))then
                                BulkVirial = BulkVirial + PairV/2.0
                                EinBulk = EinBulk + PairE/2.0
                           endif
                           if((FLAG1.NE.SubFGas).AND.(FLAG2.EQ.SubFGas))then
                                BulkVirial = BulkVirial + PairV/2.0
                                EinBulk = EinBulk + PairE/2.0
                           endif
                     
              endif               
            ENDDO    ! DO j
         ENDDO  ! Do i
 

!       ------------------------------------------------------
!       THE PAIRWISE INTERACTION BETWEEN FLUID AND SOLID
!       ------------------------------------------------------
        DO i = 1, N
           IF( SURFACE .OR. (PORE .AND. (MCy(i).GE. AccVYL) .AND. (MCy(i).LE. AccVYH) ) )THEN
              if(steele)then
                   FLAG1=indNinSBox(i)
                   IF(PORE.AND. ExtraBox)THEN
                       IF(FLAG1.eq.SubFGas)THEN
                          UC=0.0E0
                       ELSE
                          call steelepotential(i,UC)
                       ENDIF
                   ELSE
                       call steelepotential(i,UC)
                   ENDIF
                       EnergyC = EnergyC + UC 
                       EnergyMatrix(i,i) = UC
                       if(SURFACE .AND. (indNinSBox(i).eq.SubFGas))EinBulk = EinBulk + UC                 
               endif
       
               if(Bojan)then
                  PEnergyC = 0.0
                  do j = 1, NS
                     call BojanPotential(j, i,BU)
                     PEnergyC = PEnergyC + BU
                  enddo
                  IF(Close1End.OR.Close2Ends)then
                      call BojanPotential(NS+2, i,BU) 
                      PEnergyC = PEnergyC + BU 
                  ENDIF                   
                  IF(PEnergyC/Temperature.GT.PairEsf_Limit)PEnergyC = PairEsf_Limit*Temperature
                  EnergyC = EnergyC + PEnergyC 
                  EnergyMatrix(i,i) = PEnergyC
                  if(SURFACE .AND. (indNinSBox(i).eq.SubFGas))EinBulk = EinBulk + PEnergyC
               endif
            ELSE
               EnergyMatrix(i,i) = 0.0E0
            ENDIF
       ENDDO
!     --------------------------------------------------------
!     CALCULATE THE ENERGY OF EACH PARTICLE WITH THE REST
!     --------------------------------------------------------
         DO i = 1, N
            do j = 1, N !! j-do
              IF(i .ne. j)THEN !! 5-if
                 Energyi(i) = Energyi(i) + EnergyMatrix(i,j)
              ENDIF
            enddo
            Energyi(i) = Energyi(i) + EnergyMatrix(i,i)
            rate(i) = exp(Energyi(i)/Temperature)
            TotRate = TotRate + rate(i)
            rateSBox(indNinSBox(i)) = rateSBox(indNinSBox(i)) + rate(i)
         ENDDO
                  
     RETURN
     END SUBROUTINE ConfigEnergyADS     
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE POTENTIALENERGY
!      THIS SUBROUTINE CALCULATE THE POTENTIALENERGY AMONG THE N PARTICLES
!      USE LENNARD-JONES 12-6 POTENTIAL MODEL 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine potentialEnergy(d2, sigma1, welldepth1, U, V, U1st)
      USE Constant_M
      USE MCSETTING_M
      implicit none
      
      real*8 d1, d2, d6, U, V, U1st
      real*8 sigma1, welldepth1
      real*8 term1, term2, term3
      
!          ----------------------------------
!          LENNARD-JONES 12-6 POTENTIAL MODEL
!          ----------------------------------
           d6  = (sigma1**2/d2)**3
           U   =  4.0E0*welldepth1*d6*(d6 - 1.0E0)
           V   =  16.0E0*welldepth1*d6*(d6 - 0.5E0)
           U1st = -24.0E0*welldepth1*(sigma1**6.0)*( 2.0*(sigma1**6)/(d2**6.5) - 1.0/(d2**3.5) )
     
     return
     end
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE LIMITPOTENTIALENERGY
!      THIS SUBROUTINE CALCULATE THE LIMIT POTENTIALENERGY WHEN THE DISTANCE
!      BETWEEN TWO PARTICLES REACH THE DISTANCE_LIMIT, E.G. 0.01
!      USE LENNARD-JONES 12-6 POTENTIAL MODEL 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine LimitPotentialEnergy(V, U1st)
      USE Constant_M
      USE MCSETTING_M
      USE PHYSICAL_M
      implicit none
      
      real*8 drevs6, equiD, V, U1st
      
!          ----------------------------------
!          LENNARD-JONES 12-6 POTENTIAL MODEL
!          ----------------------------------
           drevs6 = (1.0 + sqrt(1.0+PairEff_Limit*Temperature))/2.0
           equiD = (1.0/drevs6)**(1.0/6.0)
           V   =  16.0E0*(1.0/(equiD**12.0) - 0.5E0/(equiD**6.0))
           U1st = -12.0E0*( 2.0/(equiD**13.0) - 1.0/(equiD**7.0) )
     return
     end



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE rotate
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          SUBROUTINE rotate
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

          RETURN
          ENDSUBROUTINE rotate
            
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE VLE RADIUS DISTRIBUTION
!      Calculate the radius distribution in the gas and liuquid phase in VLE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
           SUBROUTINE VLEradiusdistribution(switch)
           USE Constant_M
           USE FLUID_SOLID_M
           USE PHYSICAL_M
           USE MCSETTING_M
           USE SUBBOX_M
           implicit none
           
           integer bin, i, j, k, ngr, maxBinVLE(3)
           integer switch,l,m,n
           integer iPSub, IDi
           
           real*8 radiusBVLE(3),  deltaBin
           real*8 r(1000),rA(1000)
           real*8 vb, nid
           real*8 mcBin(3,1000), mcSum_Nbin(3,1000), Av_mcBin(3,1000)
           real*8 mcxij,mcyij,mczij
           real*8 rmcx,rmcy,rmcz,mcradiusD
           
           
           common/gr/ngr

           
           deltaBin   = 0.1e0/scalelength
           do iPSub=1,3
              radiusBVLE(iPSub) = VLEBoxLengthX/2.0e0
              if(radiusBVLE(iPSub).LT.5.0)radiusBVLE(iPSub)=5.0
              maxBinVLE(iPSub)     = int(radiusBVLE(iPSub)/deltaBin) 
           enddo
           
           if(switch.eq.1)then
              mcSum_NBin   = 0.0e0
              Av_mcBin     = 0.0e0 
              ngr        = 0
              open(unit=21, file='radiusVLE',status='unknown')
              write(21,*)'Section#   Distance(A)   g(r)-mc'
           endif
         
           if(switch.eq.2)then
              ngr      = ngr + 1
              mcBin    = 0.0e0
              
              DO iPSub = 1,3
                 do IDi = 1, NinPropSub(iPSub)-1 
                    i = SeqNinPropSub(iPSub,IDi)  
                    do j = 1, Npart
                      IF(j.NE.i)THEN
                         mcxij = abs(mcx(i) - mcx(j))
                         mcyij = abs(mcy(i) - mcy(j))
                         mczij = abs(mcz(i) - mcz(j))
         !-----------------------------------------
         !  Radius distribution for Mass Center
         !-----------------------------------------
                         rmcx  = mcxij
                         rmcy  = mcyij
                         rmcz  = mczij                            
                         if(pbcx)then
                            if(mcxij .gt. VLEBoxLengthX/2.0e0) rmcx = rmcx - VLEBoxLengthX
                         endif
                         if(pbcy)then
                            if(mcyij .gt. VLEBoxLengthY/2.0e0) rmcy = rmcy - VLEBoxLengthY
                         endif
                         if(pbcz)then
                            if(mczij .gt. VLEBoxLengthZ/2.0e0) rmcz = rmcz - VLEBoxLengthZ 
                         endif              
                         mcradiusD = sqrt(rmcx*rmcx + rmcy*rmcy + rmcz*rmcz)
                         if(mcradiusD .le. radiusBVLE(iPSub))then
                            bin = int(mcradiusD/deltaBin)+1
                            mcBin(iPSub,bin) = mcBin(iPSub,bin) + 1.0e0
                         endif  
                      ENDIF               
                   enddo
                enddo

                do k = 1,maxbinVLE(iPSub)
                   if(NinPropSub(iPSub).gt.0)mcsum_NBin(iPSub,k) = mcsum_NBin(iPSub,k) + mcBin(iPSub,k)/dble(NinPropSub(iPSub))
                enddo   
             ENDDO           
           endif

         
           if(switch.eq.3)then 
             DO iPSub=1,3
               do k = 1, maxbinVLE(iPSub)
                  r(k)   = deltaBin*(k-0.5)
                  rA(k)  = r(k)*scalelength
                  vb     = (k**3 - (k-1.0e0)**3)*deltaBin**3
                  nid    = (4.0e0/3.0e0)*pi*vb*AveProprho(iPSub)
                  Av_mcbin(iPSub,k) =  mcSum_Nbin(iPSub,k)/(ngr*nid)
                  write(21,'(I10,2e16.7)')iPSub,rA(k),Av_mcbin(iPSub,k)
               enddo
             ENDDO
           endif
         
         RETURN
         END SUBROUTINE VLEradiusdistribution         
           

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   SUBROUTINE COULOMB FORCE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
        subroutine coulombforce(cld, charge1,charge2,CLU,CLV,CLU1st)
        USE Constant_M
        USE MCSETTING_M
        implicit none
        
        real*8 cld,charge1,charge2, CLU,CLV,CLU1st
        
                
!        if((cld.lt.0.4e0) .AND. (charge1*charge2.lt.0.0e0))then
!           CLU = 1.0e10
!           return
!        else       
        CLU = charge1*charge2/(4.0e0*pi*permittivity*Scalelength*1.0e-10*scaleEnergy*kB*cld)
        CLV = CLU/3.0
        CLU1st = -CLU/cld
!        endif
        return
        end
     
           
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE InsertVLE
!    Iinitialize the constant Npart in the liquid section
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE InsertVLE   
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer i,j, k, iSub
      integer N, N1, N3
      integer iCount, iSelected 
      integer iPSub, FLAG

      real*8 rdn
      real*8 BoxLength1, BoxLength2, BoxLength3
      real*8 iMCx0(50000), iMCy0(50000), iMCz0(50000)
     
    
     N = Npart
     if(N .eq. 0)return     
!    --------------------------------------
!    NUMBER OF POSITION IN LINEAR DIMENSION
!    --------------------------------------
     N1         = int(N**(1.E0/3.E0)) + 1
     BoxLength1 = (VLEBoxlengthX)/N1
     IF(HomoStart)THEN
        BoxLength2 = (VLEBoxlengthY)/N1
     ELSE
        BoxLength2 = (LiqHy-LiqLy)/N1
     ENDIF
     BoxLength3 = (VLEBoxlengthZ)/N1
!    -------------------------------------------
!    NUMBER OF AVAILABLE POSITIONS IN CUBE IS N3
!    -------------------------------------------
     N3 = N1*N1*N1
!     if(N3.gt.N)N=N3
!    ------------------------------------------
!    ASSIGN THE COORDINATES OF THE N3 POSITIONS
!    ------------------------------------------
     do k=1, N1
        do j=1, N1
              do i=1, N1
                 iCount = ( (k-1)*N1+ (j-1))*N1 + i
                 iMCx0(iCount) = BoxLength1/2.0 + (i-1)*Boxlength1 
                 IF(HomoStart)THEN
                   iMCy0(iCount) = BoxLength2/2.0 + (j-1)*Boxlength2
                 ELSE
                   iMCy0(iCount) = BoxLength2/2.0 + (j-1)*Boxlength2 + LiqLy
                 ENDIF
                 iMCz0(iCount) = BoxLength3/2.0 + (k-1)*Boxlength3 
              enddo
        enddo
     enddo
!    ------------------------------------------------------------
!    NOW RANDOMLY DISTRIBUTE N PARTICLES IN N3 POSSIBLE LOCATIONS
!    ------------------------------------------------------------
     DO i=1, N
        call random_number(rdn)
         iSelected       = int(rdn*N3) + 1
         if(iSelected.gt.N3) iSelected = N3
          MCx(i)        = iMCx0(iSelected)
          MCy(i)        = iMCy0(iSelected)
          MCz(i)        = iMCz0(iSelected)
                    
          CALL ROTATE
          do j = 1, NLJ
             LJx(i,j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j) + MCx(i) - MCx0
             LJy(i,j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j) + MCy(i) - MCy0
             LJz(i,j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j) + MCz(i) - MCz0
          enddo
          do j = 1, NCL
             CLx(i,j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) + MCx(i) - MCx0
             CLy(i,j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) + MCy(i) - MCy0
             CLz(i,j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) + MCz(i) - MCz0
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
           
                
        do iSub=1, SubBox
          if((MCy(i).GE.LayerYL(iSub)).AND.(MCy(i).LT.LayerYH(iSub)))then
             NinSBox(iSub)= NinSBox(iSub)+1
             SeqNinSBox(iSub,NinSBox(iSub))=i
             indNinSBox(i)=iSub
             EXIT
          endif
        enddo
        
        FLAG=0
        do iPSub = 1, 3
          if((MCy(i).GE.PropLy(iPSub)).AND.(MCy(i).LT.PropHy(iPSub)))then
            NinPropSub(iPSub)= NinPropSub(iPSub) + 1
            SeqNinPropSub(iPSub,NinPropSub(iPSub))= i
            indNinPropSub(i) = iPSub
            FLAG=1
            EXIT
          endif
        enddo
        IF(FLAG.EQ.0)indNinPropSub(i) = 0
     ENDDO
     
     return
     END SUBROUTINE InsertVLE
     

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE IniSolidSlab
!    Construct the solid slab in the simulation box
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE IniSolidSlab  
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer i,j, k, m
      integer N, iSub
      integer iPSub, FLAG

      real*8 rdn
     
    
     Npart = 0   
     N = 0
!    ---------------------------------------------------------
!    Locate the pariticles in the plane parallel to YZ surface
!    --------------------------------------------------------
     
     do i = 1, NunitX
        do j = 1, NunitY
           do k = 1, NunitZ
              N = (i-1)*(NunitY*NunitZ) + (j-1)*NunitZ + k
              MCx(N) = (i-1)*LunitX
              MCy(N) = (j-1)*LunitY + LunitY/2.0 + LiqLy
              MCz(N) = (k-1)*LunitZ + LunitZ/2.0
              do m = 1, NLJ
                 LJx(N,m) = LJx0(m) + MCx(N) - MCx0
                 LJy(N,m) = LJy0(m) + MCy(N) - MCy0
                 LJz(N,m) = LJz0(m) + MCz(N) - MCz0
              enddo
         
              do m =1, NCL
                CLx(N,m) = CLx0(m) + MCx(N) - MCx0
                CLy(N,m) = CLy0(m) + MCy(N) - MCy0
                CLz(N,m) = CLz0(m) + MCz(N) - MCz0
             enddo
           enddo
        enddo
     enddo
     Npart = N
!    ---------------------------------------------------------
!    Locate the pariticles in the plane parallel to XZ surface
!    --------------------------------------------------------
     do j = 1, NunitY
        do i = 1, NunitX 
           do k = 1, NunitZ
              N = Npart + (j-1)*(NunitX*NunitZ) + (i-1)*NunitZ + k
              MCx(N) = (i-1)*LunitX + LunitX/2.0
              MCy(N) = (j-1)*LunitY + LiqLy
              MCz(N) = (k-1)*LunitZ + LunitZ/2.0
              do m = 1, NLJ
                 LJx(N,m) = LJx0(m) + MCx(N) - MCx0
                 LJy(N,m) = LJy0(m) + MCy(N) - MCy0
                 LJz(N,m) = LJz0(m) + MCz(N) - MCz0
              enddo
         
              do m =1, NCL
                CLx(N,m) = CLx0(m) + MCx(N) - MCx0
                CLy(N,m) = CLy0(m) + MCy(N) - MCy0
                CLz(N,m) = CLz0(m) + MCz(N) - MCz0
             enddo
           enddo
        enddo
     enddo
     Npart = N

!    ---------------------------------------------------------
!    Locate the pariticles in the plane parallel to XY surface
!    --------------------------------------------------------
     do k = 1, NunitZ
        do i = 1, NunitX 
           do j = 1, NunitY
              N = Npart + (k-1)*(NunitX*NunitY) + (i-1)*NunitY + j
              MCx(N) = (i-1)*LunitX + LunitX/2.0
              MCy(N) = (j-1)*LunitY + LunitY/2.0 + LiqLy
              MCz(N) = (k-1)*LunitZ 
              do m = 1, NLJ
                 LJx(N,m) = LJx0(m) + MCx(N) - MCx0
                 LJy(N,m) = LJy0(m) + MCy(N) - MCy0
                 LJz(N,m) = LJz0(m) + MCz(N) - MCz0
              enddo
         
              do m =1, NCL
                CLx(N,m) = CLx0(m) + MCx(N) - MCx0
                CLy(N,m) = CLy0(m) + MCy(N) - MCy0
                CLz(N,m) = CLz0(m) + MCz(N) - MCz0
             enddo
           enddo
        enddo
     enddo
     Npart = N   
      
!    ---------------------------------------------------------
!    Locate the pariticles in the corner of the cubic
!    --------------------------------------------------------
     do i = 1, NunitX
        do j = 1, NunitY 
           do k = 1, NunitZ
              N = Npart + (i-1)*(NunitY*NunitZ) + (j-1)*NunitZ + k
              MCx(N) = (i-1)*LunitX 
              MCy(N) = (j-1)*LunitY + LiqLy
              MCz(N) = (k-1)*LunitZ 
           enddo
        enddo
     enddo
     Npart = N
     DO i = 1, Npart
          CALL ROTATE
          do j = 1, NLJ
             LJx(i,j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j) + MCx(i) - MCx0
             LJy(i,j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j) + MCy(i) - MCy0
             LJz(i,j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j) + MCz(i) - MCz0
          enddo

          do j = 1, NCL
             CLx(i,j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) + MCx(i) - MCx0
             CLy(i,j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) + MCy(i) - MCy0
             CLz(i,j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) + MCz(i) - MCz0
          enddo  
     
        do iSub=1, SubBox
          if((MCy(i).GE.LayerYL(iSub)).AND.(MCy(i).LT.LayerYH(iSub)))then
             NinSBox(iSub)= NinSBox(iSub)+1
             SeqNinSBox(iSub,NinSBox(iSub))=i
             indNinSBox(i)=iSub
             EXIT
          endif
        enddo
        
        FLAG=0
        do iPSub = 1, 3
          if((MCy(i).GE.PropLy(iPSub)).AND.(MCy(i).LT.PropHy(iPSub)))then
            NinPropSub(iPSub)= NinPropSub(iPSub) + 1
            SeqNinPropSub(iPSub,NinPropSub(iPSub))= i
            indNinPropSub(i) = iPSub
            FLAG=1
            EXIT
          endif
        enddo
        IF(FLAG.EQ.0)indNinPropSub(i) = 0
     ENDDO
     
     return
     ENDSUBROUTINE IniSolidSlab



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE LiqPosCheck
!    Check the position of the fluid, make sure the center is in the middle of 
!    the simulation box
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       SUBROUTINE LiqPosCheck
       USE Constant_M
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       IMPLICIT NONE
    
       integer i, j, iSub, FLAG, iPSub
    
       real*8 SumMCY, AveMCY, DifMCY, MCyold
    
       SumMCY = 0.0  
       do i=1, Npart
          SumMCY = SumMCY+MCy(i)
       enddo
       AveMCY = SumMCY/dble(Npart)
       DifMCY = VLEBoxlengthY/2.0E0 - AveMCY
    
       IF(abs(DifMCY).GT.(1.0E-3))THEN
           NinSBox = 0
           SeqNinSBox =0
           indNinSBox =0
           NinPropSub = 0
           SeqNinPropSub = 0
           indNinPropSub = 0
           DO i=1, Npart
                MCy(i)=MCy(i)+DifMCY         
                do j = 1, NLJ
                    LJy(i,j) = LJy(i,j)+DifMCY
                enddo      
                do j = 1, NCL
                    CLy(i,j) = CLy(i,j)+DifMCY
                enddo          
               if(MCy(i).gt.VLEBoxlengthY)then
                   MCyold = MCy(i)
                   MCy(i) = MCy(i)-int(MCy(i)/VLEBoxlengthY)*VLEBoxlengthY
                   do j = 1, NLJ
                      LJy(i,j) = LJy(i,j) + (MCy(i) - MCyold)
                   enddo      
                   do j = 1, NCL
                      CLy(i,j) = CLy(i,j) + (MCy(i) - MCyold)
                   enddo 
               elseif(MCy(i).lt.0.0)then
                   MCyold = MCy(i)
                   MCy(i) = MCy(i)-int(MCy(i)/VLEBoxlengthY)*VLEBoxlengthY + VLEBoxlengthY
                   do j = 1, NLJ
                      LJy(i,j) = LJy(i,j)+ (MCy(i) - MCyold)
                   enddo      
                   do j = 1, NCL
                      CLy(i,j) = CLy(i,j)+ (MCy(i) - MCyold)
                   enddo 
               endif     
         
               do iSub=1, SubBox
                  if((MCy(i).GE.LayerYL(iSub)).AND.(MCy(i).LT.LayerYH(iSub)))then
                     NinSBox(iSub)= NinSBox(iSub)+1
                     SeqNinSBox(iSub,NinSBox(iSub))=i
                     indNinSBox(i)=iSub
                     EXIT
                  endif
               enddo   
          
              FLAG=0
              do iPSub = 1, 3
               if((MCy(i).GE.PropLy(iPSub)).AND.(MCy(i).LT.PropHy(iPSub)))then
                 NinPropSub(iPSub)= NinPropSub(iPSub) + 1
                 SeqNinPropSub(iPSub,NinPropSub(iPSub))= i
                 indNinPropSub(i) = iPSub
                 FLAG=1
                 EXIT
               endif
              enddo
              IF(FLAG.EQ.0)indNinPropSub(i) = 0
          ENDDO
     ENDIF
       
     RETURN
     ENDSUBROUTINE LiqPosCheck
    
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE WIDOM
!    Check the position of the fluid, make sure the center is in the middle of 
!    the simulation box
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       SUBROUTINE WIDOM
       USE Constant_M
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       IMPLICIT NONE
    
       integer i, j, iSub, new
    
       real*8 rdn
       real*8 EnergyGhost, CLEnergyGhost, EnergyCGhost, TEnergyGhost
       real*8 EnergyLRC, VirialLRC
       real*8 WGhost
    
      
       new = Npart+1
       DO iSub =1, SubBox
         WGhost = 0.0
         rhoFW(iSub)=dble(NinSBox(iSub))/SubVol(iSub)
         do i=1, NGhost
            call random_number(rdn)
            MCx(new) = rdn*VLEBoxlengthX
            call random_number(rdn)
            MCy(new) =  LayerYL(iSub) + rdn*(LayerYH(iSub)-LayerYL(iSub)) 
            call random_number(rdn)
            MCz(new) = rdn*VLEBoxlengthZ
          
            CALL ROTATE
            do j = 1, NLJ
               LJx(new,j) = Rot(1)*LJx0(j)+ Rot(4)*LJy0(j) + Rot(7)*LJz0(j) + MCx(new) - MCx0
               LJy(new,j) = Rot(2)*LJx0(j)+ Rot(5)*LJy0(j) + Rot(8)*LJz0(j) + MCy(new) - MCy0
               LJz(new,j) = Rot(3)*LJx0(j)+ Rot(6)*LJy0(j) + Rot(9)*LJz0(j) + MCz(new) - MCz0
            enddo

            do j = 1, NCL
               CLx(new,j) = Rot(1)*CLx0(j)+ Rot(4)*CLy0(j) + Rot(7)*CLz0(j) + MCx(new) - MCx0
               CLy(new,j) = Rot(2)*CLx0(j)+ Rot(5)*CLy0(j) + Rot(8)*CLz0(j) + MCy(new) - MCy0
               CLz(new,j) = Rot(3)*CLx0(j)+ Rot(6)*CLy0(j) + Rot(9)*CLz0(j) + MCz(new) - MCz0
            enddo  
          
            call energysingleparticle(new,EnergyGhost, CLEnergyGhost)
          
            call ZLRC(rhoFW(iSub), EnergyLRC, VirialLRC)
            TEnergyGhost = EnergyGhost + CLEnergyGhost + 2.0*EnergyLRC
            WGhost = WGhost + exp(-TEnergyGhost/Temperature)
         enddo
          
          ExcessChemicalPVLE(iSub) = -Temperature*LOG( WGhost/dble(NGhost) )
          IdealChemicalPVLE(iSub) = Temperature*LOG(deBro3*(NinSBox(iSub)+1)/SubVol(iSub))
          ChemicalPVLE(iSub) = IdealChemicalPVLE(iSub) + ExcessChemicalPVLE(iSub)

       ENDDO 
         
       RETURN
       END SUBROUTINE WIDOM
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE FLUCTVLE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
      SUBROUTINE FLUCTVLE(switch)
      USE CONSTANT_M
      USE PHYSICAL_M
      USE MCSETTING_M
      USE SUBBOX_M
      IMPLICIT NONE
     
      integer iSub,switch
     
      real*8, Dimension(:),POINTER:: sum_Ninbin,sum_NNinBin
      real*8, Dimension(:),POINTER:: ave_Ninbin, ave_NNinBin
      real*8, Dimension(:),POINTER:: Fun_NN_bin, LocalNFluctbin
     
        if(switch.eq.0)then
          Allocate(sum_Ninbin(SubBox), sum_NNinBin(SubBox) )
          Allocate(ave_Ninbin(SubBox), ave_NNinBin(SubBox))
          Allocate(Fun_NN_bin(SubBox), LocalNFluctbin(SubBox))
          return
        endif

        if(switch.eq.1)then
          sum_Ninbin         = 0.0e0    !!* 
          sum_NNinBin        = 0.0e0    !!*
          CountN              = 0.0E0    !!*
          return
        endif
     
        IF(switch.eq.2)THEN
           CountN = CountN + 1
           do iSub =1, SubBox
               Sum_Ninbin(iSub)  = Sum_Ninbin(iSub) + dble(NinSBox(iSub))
               Sum_NNinbin(iSub) = Sum_NNinbin(iSub) + dble(NinSBox(iSub))*dble(NinSBox(iSub))
           enddo
           return
        ENDIF
     
        if(switch.eq.3)then   
           do iSub =1, SubBox
               ave_Ninbin(iSub)     = sum_Ninbin(iSub)/dble(CountN) !!*
               ave_NNinBin(iSub)    = sum_NNinBin(iSub)/dble(CountN) !!*
               Fun_NN_bin(iSub)     = ave_NNinBin(iSub) - ave_Ninbin(iSub)*ave_Ninbin(iSub) !!*
               LocalNFluctbin(iSub) = Fun_NN_bin(iSub)/ave_Ninbin(iSub)
           enddo
           write(25,*)'------------Local number fluctuation------------------'
           write(25,'(100E16.6)')(LocalNFluctbin(iSub),iSub=1,SubBox) 
           write(25,'(100E16.6)')(ave_Ninbin(iSub),iSub=1,SubBox) 
           return
         endif   
         
         if(switch.eq.4)then
           DEAllocate(sum_Ninbin,sum_NNinBin)
           DEAllocate(ave_Ninbin, ave_NNinBin)
           DEAllocate(Fun_NN_bin, LocalNFluctbin)
           return
         endif
     
     RETURN
     END SUBROUTINE FLUCTVLE    
    
    
    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE DISTANCE DISTRIBUTION VLE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     SUBROUTINE VLEdistancedistribution(switch)
     USE Constant_M
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE SUBBOX_M
     implicit none
     
     integer i,j,k
     integer bin, switch, NDD, maxbinDVLE
     
     real*8 deltaY, distanceY, ChangeRatio
     real*8, Dimension(:),POINTER::sum_DNbin, ave_DNbin, DNbin
     real*8, Dimension(:),POINTER::dy,dyA,DDensity,DDensityKmolperM3
     
     common/dd/NDD
    
     deltaY         = 0.1e0/scalelength
     maxbinDVLE     = int(VLEboxLengthY0/deltaY)+ 1 
     if(switch.eq.1)then
       Allocate(sum_DNbin(maxbinDVLE), ave_DNbin(maxbinDVLE), DNbin(maxbinDVLE))
       Allocate(dy(maxbinDVLE),dyA(maxbinDVLE),DDensity(maxbinDVLE),DDensityKmolperM3(maxbinDVLE))
       sum_DNbin = 0.0e0
       ave_DNbin = 0.0e0
       DNbin     = 0.0e0
       NDD       = 0
       open(unit=39, file='distanceDVLE',status='unknown')
     endif
     
     if(switch.eq.2)then
        NDD = NDD + 1
        DNbin = 0.0e0  
        ChangeRatio = VLEboxLengthY/VLEboxLengthY0
        DO i = 1, Npart
           DistanceY  = abs(mcy(i))/ChangeRatio
           bin = int(DistanceY/deltaY)+1     
           DNbin(bin) = DNbin(bin) + 1.0e0 
        ENDDO     
        
        do i = 1, maxbinDVLE 
          sum_DNbin(i) = sum_DNbin(i) + DNbin(i)*InterTime
        enddo   
     endif

     if(switch.eq.3)then
        ChangeRatio = VLEboxLengthY/VLEboxLengthY0
        deltaY = deltaY*ChangeRatio
        do i = 1, maxbinDVLE 
          dy(i)  = deltaY*(i-0.5)
          dyA(i) = dy(i)*scalelength
          ave_DNbin(i) = sum_DNbin(i)/TotTime
          DDensity(i)  = ave_DNbin(i)/(VLEBoxlengthX*VLEBoxlengthZ*deltaY)
          DDensityKmolperM3(i) = DDensity(i)/6.0230E23/(Scalelength*1.0e-10)**3/1.0e3
          write(39,'(I4,2E16.6)')Npart, dyA(i),DDensityKmolperM3(i)
        enddo
             
        DEAllocate(sum_DNbin, ave_DNbin, DNbin)
        DEAllocate(dy,dyA,DDensity,DDensityKmolperM3)
     endif
     
   RETURN
   ENDSUBROUTINE VLEdistancedistribution  
