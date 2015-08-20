!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE SHAPECHANGE
!      Adjust the shape of the cross section area (XZ)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     SUBROUTINE ShapeChange(TEnergyOld, EnergyOld, CLEnergyOld,TEnergyNew, EnergyNew, CLEnergyNew)
     USE Constant_M
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE SUBBOX_M
     implicit none
     
     integer i,j,iSub
     integer FLAG1, FLAG2, FLAG3
     
     real*8 rdn1,rdn2, rdn
     real*8 maxChangeL, NewLinX, NewLinY, NewLinZ
     real*8 ChangeRatio1, ChangeRatio2
     real*8 TEnergyOld, EnergyOld, CLEnergyOld
     real*8 TEnergyNew, EnergyNew, CLEnergyNew
     real*8 Energy, CLEnergy,EnergyC
           
         maxChangeL = 0.01  ! Times COLLOSION DIAMETER
         call random_number(rdn1)
         IF(rdn1.LE.0.5)THEN  ! CHANGE THE DIMENTION IN X- DIRECTION
                 call random_number(rdn2)
                 NewLinX = VLEBoxlengthX  + (rdn2-0.5)*maxChangeL
!                 if(abs(NewLinX-VLEBoxlengthX0).gt.(1.0/2.0))CYCLE
                 ChangeRatio1 = NewLinX/VLEBoxlengthX
                 NewLinY = VLEBoxlengthY*VLEBoxlengthX/NewLinX ! Keep the volume constant
                 ChangeRatio2 = NewLinY/VLEBoxlengthY
              
                 DO i=1, Npart
                    MCx(i) = MCx(i)*ChangeRatio1
                    MCy(i) = MCy(i)*ChangeRatio2              
                    do j = 1, NLJ
                       LJx(i,j) = LJx(i,j) + MCx(i) - MCx(i)/ChangeRatio1
                       LJy(i,j) = LJy(i,j) + MCy(i) - MCy(i)/ChangeRatio2
                    enddo
                    do j = 1, NCL
                       CLx(i,j) = CLx(i,j) + MCx(i) - MCx(i)/ChangeRatio1
                       CLy(i,j) = CLy(i,j) + MCy(i) - MCy(i)/ChangeRatio2
                    enddo
                 ENDDO
                 VLEBoxlengthX = NewLinX
                 VLEBoxlengthY = NewLinY
                            
                 call ConfigEnergy(Energy, CLEnergy)
                 TEnergyNew = Energy              
                 EnergyNew = Energy
                 CLEnergyNew = CLEnergy
                 do iSub = 1, SubBox
                    LayerXL(iSub) = LayerXL(iSub)*ChangeRatio1
                    LayerXH(iSub) = LayerXH(iSub)*ChangeRatio1
                    LayerYL(iSub) = LayerYL(iSub)*ChangeRatio2
                    LayerYH(iSub) = LayerYH(iSub)*ChangeRatio2
                    SStepX(iSub)  = SStepX(iSub)*ChangeRatio1
                    SStepY(iSub)  = SStepY(iSub)*ChangeRatio2
                    SubVol(iSub)  = (LayerXH(iSub)-LayerXL(iSub))*(LayerYH(iSub)-LayerYL(iSub))*(LayerZH(iSub)-LayerZL(iSub))
                 enddo   
                 do i=1, 3
                    PropLy(i) = PropLy(i)*ChangeRatio2
                    PropHy(i) = PropHy(i)*ChangeRatio2
                 enddo 
         ELSE  !CHANGE THE DIMENTION IN Z- DIRECTION
                 call random_number(rdn2)
                 NewLinZ = VLEBoxlengthZ  + (rdn2-0.5)*maxChangeL
  !               if(abs(NewLinZ-VLEBoxlengthZ0).gt.(1.0/2.0))CYCLE
                 ChangeRatio1 = NewLinZ/VLEBoxlengthZ
                 NewLinY = VLEBoxlengthY*VLEBoxlengthZ/NewLinZ ! Keep the volume constant
                 ChangeRatio2 = NewLinY/VLEBoxlengthY
              
                 DO i=1, Npart
                    MCz(i) = MCz(i)*ChangeRatio1
                    MCy(i) = MCy(i)*ChangeRatio2
                    do j = 1, NLJ
                       LJz(i,j) = LJz(i,j) + MCz(i) - MCz(i)/ChangeRatio1
                       LJy(i,j) = LJy(i,j) + MCy(i) - MCy(i)/ChangeRatio2
                    enddo
                    do j = 1, NCL
                       CLz(i,j) = CLz(i,j) + MCz(i) - MCz(i)/ChangeRatio1
                       CLy(i,j) = CLy(i,j) + MCy(i) - MCy(i)/ChangeRatio2
                    enddo
                 ENDDO
                 VLEBoxlengthZ = NewLinZ
                 VLEBoxlengthY = NewLinY
                       
                 call ConfigEnergy(Energy, CLEnergy)
                 TEnergyNew = Energy  
                 EnergyNew = Energy
                 CLEnergyNew = CLEnergy
                 do iSub = 1, SubBox
                    LayerZL(iSub) = LayerZL(iSub)*ChangeRatio1
                    LayerZH(iSub) = LayerZH(iSub)*ChangeRatio1
                    LayerYL(iSub) = LayerYL(iSub)*ChangeRatio2
                    LayerYH(iSub) = LayerYH(iSub)*ChangeRatio2
                    SStepZ(iSub)  = SStepZ(iSub)*ChangeRatio1
                    SStepY(iSub)  = SStepY(iSub)*ChangeRatio2
                    SubVol(iSub)  = (LayerXH(iSub)-LayerXL(iSub))*(LayerYH(iSub)-LayerYL(iSub))*(LayerZH(iSub)-LayerZL(iSub))
                 enddo
                 do i=1, 3
                    PropLy(i) = PropLy(i)*ChangeRatio2
                    PropHy(i) = PropHy(i)*ChangeRatio2
                 enddo
         ENDIF
         
      RETURN
      ENDSUBROUTINE ShapeChange
  
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ThermalPropertyChange
!      UPdate the Terms for P, Surface Tension and Evaporation Heat calculation
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
        SUBROUTINE  ThermalPropertyChange(switch, i)
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE SUBBOX_M
        implicit none
      
        integer switch,i,j,iflag
        integer FLAG1, FLAG2, FLAG3
        integer iPSub
        
        
        real*8 EnergyOld, CLEnergyOld
        real*8 EnergyNew, CLEnergyNew
        real*8 TEnergy, TEnergyOld, TEnergyNew
        real*8 Energy,CLEnergy
        real*8 VirialNew, VirialOld, VirialNewF(3), VirialOldF(3)
        real*8 EnergyNewF(3), EnergyOldF(3)
        real*8 FirstDerivEOld,FirstDerivENew
        real*8 FirstDerivEOld_Ale,FirstDerivENew_Ale
        
        common/FLAGSING/FLAG1
        common/SWITCH12/TEnergyOld,EnergyOld,CLEnergyOld,VirialOld,VirialOldF,EnergyOldF
        common/SWITCH13/FirstDerivEOld,FirstDerivEOld_Ale
        
        
          IF(switch.eq.1)THEN
                FLAG1=indNinPropSub(i)
                 call energySingleParticle(i, EnergyOld, CLEnergyOld)
                 TEnergyOld = EnergyOld
                 VirialOld = Virial
                 FirstDerivEOld = DerivE
                 FirstDerivEOld_Ale = DerivE_Ale
                 do iflag=1,3
                    VirialOldF(iflag)=VirialF(iflag)
                    EnergyOldF(iflag)=EnergyF(iflag)
                 enddo
          ENDIF
        
        
          IF(switch.eq.2)THEN
              indNinPropSub(i)=0
              do iPSub = 1, 3
                 if((MCy(i).GE.PropLy(iPSub)).AND.(MCy(i).LT.PropHy(iPSub)))then
                     indNinPropSub(i) = iPSub
                     EXIT
                 endif 
              enddo
              FLAG2=indNinPropSub(i)
              call energysingleparticle(i,EnergyNew, CLEnergyNew)
                    TEnergyNew = EnergyNew
                    Energy   = Energy  + EnergyNew - EnergyOld 
                    TEnergy  = TEnergy + TEnergyNew - TEnergyOld 
                    VirialNew = Virial
                    TVirial = TVirial + VirialNew - VirialOld
                    FirstDerivENew = DerivE
                    FirstDerivENew_Ale = DerivE_Ale
                    FirstDerivE = FirstDerivE + FirstDerivENew - FirstDerivEOld
                    FirstDerivE_Ale = FirstDerivE_Ale + FirstDerivENew_Ale - FirstDerivEOld_Ale
                    do iflag=1,3
                       VirialNewF(iflag)=VirialF(iflag)
                       EnergyNewF(iflag)=EnergyF(iflag)
                    enddo
                    if((FLAG1.GT.0).AND.(FLAG2.EQ.0))then
                          VirialProp(FLAG1)=VirialProp(FLAG1) - VirialOldF(FLAG1)-(VirialOld-VirialOldF(FLAG1))/2.0 + VirialNewF(FLAG1)/2.0
                          if(FLAG1.NE.1)VirialProp(1)=VirialProp(1)-VirialOldF(1)/2.0+VirialNewF(1)/2.0
                          if(FLAG1.NE.2)VirialProp(2)=VirialProp(2)-VirialOldF(2)/2.0+VirialNewF(2)/2.0
                          if(FLAG1.NE.3)VirialProp(3)=VirialProp(3)-VirialOldF(3)/2.0+VirialNewF(3)/2.0
                          
                          EnergyProp(FLAG1)=EnergyProp(FLAG1) - EnergyOldF(FLAG1)-(EnergyOld-EnergyOldF(FLAG1))/2.0 + EnergyNewF(FLAG1)/2.0
                          if(FLAG1.NE.1)EnergyProp(1)=EnergyProp(1)-EnergyOldF(1)/2.0+EnergyNewF(1)/2.0
                          if(FLAG1.NE.2)EnergyProp(2)=EnergyProp(2)-EnergyOldF(2)/2.0+EnergyNewF(2)/2.0
                          if(FLAG1.NE.3)EnergyProp(3)=EnergyProp(3)-EnergyOldF(3)/2.0+EnergyNewF(3)/2.0
                          
                          do j=1, NinPropSub(FLAG1)
                             if(SeqNinPropSub(FLAG1,j).eq.i)then
                                SeqNinPropSub(FLAG1,j)=SeqNinPropSub(FLAG1,NinPropSub(FLAG1))
                                Exit
                             endif
                          enddo
                          NinPropSub(FLAG1)=NinPropSub(FLAG1)-1
                     elseif((FLAG1.GT.0).AND.(FLAG2.GT.0))then
                          IF(FLAG1.NE.FLAG2)THEN
                             VirialProp(FLAG1)=VirialProp(FLAG1) - VirialOldF(FLAG1)-(VirialOld-VirialOldF(FLAG1))/2.0 + VirialNewF(FLAG1)/2.0
                             VirialProp(FLAG2)=VirialProp(FLAG2) - VirialOldF(FLAG2)/2.0 + (VirialNew-VirialNewF(FLAG2))/2.0 + VirialNewF(FLAG2)
                             FLAG3=6 - FLAG1 - FLAG2
                             VirialProp(FLAG3)=VirialProp(FLAG3)- VirialOldF(FLAG3)/2.0+ VirialNewF(FLAG3)/2.0
                             
                             EnergyProp(FLAG1)=EnergyProp(FLAG1) - EnergyOldF(FLAG1)-(EnergyOld-EnergyOldF(FLAG1))/2.0 + EnergyNewF(FLAG1)/2.0
                             EnergyProp(FLAG2)=EnergyProp(FLAG2) - EnergyOldF(FLAG2)/2.0 + (EnergyNew-EnergyNewF(FLAG2))/2.0 + EnergyNewF(FLAG2)
                             EnergyProp(FLAG3)=EnergyProp(FLAG3)- EnergyOldF(FLAG3)/2.0+ EnergyNewF(FLAG3)/2.0
                             
                             do j=1, NinPropSub(FLAG1)
                               if(SeqNinPropSub(FLAG1,j).eq.i)then
                                  SeqNinPropSub(FLAG1,j)=SeqNinPropSub(FLAG1,NinPropSub(FLAG1))
                                  Exit
                                endif
                             enddo
                             NinPropSub(FLAG1)=NinPropSub(FLAG1)-1
                             NinPropSub(FLAG2)=NinPropSub(FLAG2)+1
                             SeqNinPropSub(FLAG2,NinPropSub(FLAG2))=i
                          ELSEIF(FLAG1.EQ.FLAG2)THEN
                             VirialProp(FLAG1)=VirialProp(FLAG1)-VirialOldF(FLAG1)-(VirialOld-VirialOldF(FLAG1))/2.0+VirialNewF(FLAG1)+(VirialNew-VirialNewF(FLAG1))/2.0
                             if(FLAG1.NE.1)VirialProp(1)=VirialProp(1)-VirialOldF(1)/2.0+VirialNewF(1)/2.0
                             if(FLAG1.NE.2)VirialProp(2)=VirialProp(2)-VirialOldF(2)/2.0+VirialNewF(2)/2.0
                             if(FLAG1.NE.3)VirialProp(3)=VirialProp(3)-VirialOldF(3)/2.0+VirialNewF(3)/2.0
                             
                             EnergyProp(FLAG1)=EnergyProp(FLAG1)-EnergyOldF(FLAG1)-(EnergyOld-EnergyOldF(FLAG1))/2.0+EnergyNewF(FLAG1)+(EnergyNew-EnergyNewF(FLAG1))/2.0
                             if(FLAG1.NE.1)EnergyProp(1)=EnergyProp(1)-EnergyOldF(1)/2.0+EnergyNewF(1)/2.0
                             if(FLAG1.NE.2)EnergyProp(2)=EnergyProp(2)-EnergyOldF(2)/2.0+EnergyNewF(2)/2.0
                             if(FLAG1.NE.3)EnergyProp(3)=EnergyProp(3)-EnergyOldF(3)/2.0+EnergyNewF(3)/2.0
                          ENDIF
                     elseif((FLAG1.EQ.0).AND.(FLAG2.GT.0))then
                          VirialProp(FLAG2)=VirialProp(FLAG2) - VirialOldF(FLAG2)/2.0 + VirialNewF(FLAG2) + (VirialNew-VirialNewF(FLAG2))/2.0
                          if(FLAG2.NE.1)VirialProp(1)=VirialProp(1)-VirialOldF(1)/2.0+VirialNewF(1)/2.0
                          if(FLAG2.NE.2)VirialProp(2)=VirialProp(2)-VirialOldF(2)/2.0+VirialNewF(2)/2.0
                          if(FLAG2.NE.3)VirialProp(3)=VirialProp(3)-VirialOldF(3)/2.0+VirialNewF(3)/2.0
                          
                          EnergyProp(FLAG2)=EnergyProp(FLAG2) - EnergyOldF(FLAG2)/2.0 + EnergyNewF(FLAG2) + (EnergyNew-EnergyNewF(FLAG2))/2.0
                          if(FLAG2.NE.1)EnergyProp(1)=EnergyProp(1)-EnergyOldF(1)/2.0+EnergyNewF(1)/2.0
                          if(FLAG2.NE.2)EnergyProp(2)=EnergyProp(2)-EnergyOldF(2)/2.0+EnergyNewF(2)/2.0
                          if(FLAG2.NE.3)EnergyProp(3)=EnergyProp(3)-EnergyOldF(3)/2.0+EnergyNewF(3)/2.0
                          
                          NinPropSub(FLAG2)=NinPropSub(FLAG2)+1
                          SeqNinPropSub(FLAG2,NinPropSub(FLAG2))=i
                     elseif((FLAG1.EQ.0).AND.(FLAG2.EQ.0))then
                          do j=1,3
                             VirialProp(j)=VirialProp(j)-VirialOldF(j)/2.0+VirialNewF(j)/2.0
                             EnergyProp(j)=EnergyProp(j)-EnergyOldF(j)/2.0+EnergyNewF(j)/2.0
                          enddo
                     endif
          ENDIF
        
        RETURN
        ENDSUBROUTINE ThermalPropertyChange
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE TherPropChangeADS
!      UPdate the Terms for Pressure calculation in the bulk phase
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
        SUBROUTINE  TherPropChangeADS(switch, i)
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE SUBBOX_M
        USE PHYSICAL_M
        implicit none
      
        integer switch,i,j,iflag
        integer FLAG1, FLAG2, FLAG3
        integer iPSub
        
        
        real*8 EnergyOld, CLEnergyOld, EnergyCOld
        real*8 EnergyNew, CLEnergyNew, EnergyCNew
        real*8 VirialNew, VirialOld, BulkVirialNewF, BulkVirialOldF
        real*8 EinBulkNewF, EinBulkOldF, sfEninBulkNewF,sfEninBulkOldF

        
        common/FLAGSING/FLAG1
        common/ADSSWITCH12/EnergyOld, VirialOld, BulkVirialOldF, EinBulkOldF, sfEninBulkOldF
        
        
          IF(switch.eq.1)THEN
                 FLAG1=indNinSBox(i)
                 call energySingleParticleADS(i, EnergyOld, CLEnergyOld)
                 VirialOld = Virial
                 BulkVirialOldF = BulkVirialF
                 EinBulkOldF = EinBulkF
                 sfEninBulkOldF = EnergyMatrix(i,i)
          ENDIF
        
        
          IF(switch.eq.2)THEN
              FLAG2=indNinSBox(i)
              call energysingleparticleADS(i,EnergyNew, CLEnergyNew)
                    VirialNew = Virial
                    BulkVirialNewF = BulkVirialF
                    EinBulkNewF = EinBulkF
                    sfEninBulkNewF = EnergyMatrix(i,i)
                    if((FLAG1.EQ.SubFGas).AND.(FLAG2.EQ.SubFGas))then
                          BulkVirial = BulkVirial - BulkVirialOldF-(VirialOld - BulkVirialOldF )/2.0 + BulkVirialNewF + (VirialNew - BulkVirialNewF)/2.0
                          EinBulk = EinBulk - EinBulkOldF - (EnergyOld - EinBulkOldF)/2.0 + EinBulkNewF + (EnergyNew-EinBulkNewF)/2.0 - sfEninBulkOldF + sfEninBulkNewF
                    elseif((FLAG1.EQ.SubFGas).AND.(FLAG2.NE.SubFGas))then
                          BulkVirial = BulkVirial - BulkVirialOldF - (VirialOld - BulkVirialOldF )/2.0 + BulkVirialNewF/2.0
                          EinBulk = EinBulk - EinBulkOldF - (EnergyOld-EinBulkOldF)/2.0 + EinBulkNewF/2.0 - sfEninBulkOldF
                    elseif((FLAG1.NE.SubFGas).AND.(FLAG2.EQ.SubFGas))then
                          BulkVirial = BulkVirial - BulkVirialOldF/2.0 + BulkVirialNewF + (VirialNew - BulkVirialNewF)/2.0
                          EinBulk = EinBulk - EinBulkOldF/2.0 + EinBulkNewF + (EnergyNew - EinBulkNewF)/2.0 + sfEninBulkNewF
                    elseif((FLAG1.NE.SubFGas).AND.(FLAG2.NE.SubFGas))then
                          BulkVirial = BulkVirial - BulkVirialOldF/2.0 + BulkVirialNewF/2.0
                          EinBulk = EinBulk - EinBulkOldF/2.0 + EinBulkNewF/2.0
                    endif  
          ENDIF
        
        RETURN
        ENDSUBROUTINE TherPropChangeADS
        
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
                   
       RETURN
       ENDSUBROUTINE ZLRC   
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE select
!      Select one particle with Ruthenbluth rule
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine select(iS)
        USE SUBBOX_M
        implicit none
        
        integer iS
        real*8 rdn,ws,cumw
        
        call random_number(rdn)
        ws=rdn*TotRate
        cumw = rate(1)
        iS=1
        
        do while(cumw.lt.ws)
           iS=iS+1
           cumw=cumw+rate(iS)
        enddo

        return
        endsubroutine
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE NBOURSelect
!      Select one particle with Ruthenbluth rule WHILE FOR MOVING IN 
!      NEIGHBOUR BINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine NBOURSelect(iS)
        USE SUBBOX_M
        USE PHYSICAL_M
        implicit none
        
        integer i,iS, fSub, tSub1, tSub2
        integer NBiS, ID, NBOURTotN, IDNBOURN(5000)
        integer FLAG1
        
        real*8 rdn,ws,cumw, TotRateBin(5000)
        real*8 NBOURTotRate, NBOURNRate(5000), NBws, NBcumw
        
        FLAG1 = 0
        DO WHILE(FLAG1.EQ.0)
           !-----------------------------------------------
           ! In Scheme 1&3: pick one bin based Rosenbluth
           !-----------------------------------------------
              TotRateBin = 0.0                                 ! Scheme 1&3 - Start 
              fSub = 1
              do i=1, NinSBox(fSub)
                 ID = SeqNinSBox(fSub,i)
                 TotRateBin(fSub) = TotRateBin(fSub) + rate(ID)
              enddo
              call random_number(rdn)
              ws=rdn*TotRate
              cumw = TotRateBin(fSub)
        
              do while(cumw.lt.ws)
                  fSub=fSub+1
                  do i=1, NinSBox(fSub)
                     ID = SeqNinSBox(fSub,i)
                     TotRateBin(fSub) = TotRateBin(fSub) + rate(ID)
                  enddo
                  cumw=cumw+TotRateBin(fSub)
              enddo                                         ! Scheme 1&3 - End
           
           !---------------------------------------
           ! In Scheme 2: pick one bin randomly
           !--------------------------------------- 
           !   call random_number(rdn)                       ! Scheme 2 - Start 
           !   fSub = int(rdn*SubBox) + 1
           !   if(fSub.GT.SubBox)fSub=SubBox                 ! Scheme 2 - END 
        
              IF(fSub.eq.1)THEN
                 tSub1 = fSub + 1
                 tSub2 = -1
              ELSEIF(fSub.eq.SubBox)THEN
                 tSub1 = fSub - 1

                 tSub2 = -1
              ELSE
                 tSub1 = fSub - 1
                 tSub2 = fSub + 1
              ENDIF
           !--------------------------------------------------------------------------------------------
           ! Scheme 1&2: pick one particle based on the Rosenbluth for selected bin AND neighbour bin(s)
           !--------------------------------------------------------------------------------------------   
           ! Scheme 3: pick one particle based on the Rosenbluth for selected bin 
           !-----------------------------------------------------------------------

               NBOURTotRate = 0.0
               NBOURTotN = 0
               IDNBOURN = 0
               DO i = 1, NinSBox(fSub)
                  ID = SeqNinSBox(fSub,i)
                  NBOURTotN = NBOURTotN + 1
                  IDNBOURN(NBOURTotN) = ID
                  NBOURNRate(NBOURTotN) = rate(ID)
                  NBOURTotRate =  NBOURTotRate + rate(ID)
               ENDDO
         !      DO i = 1, NinSBox(tSub1)                            ! Scheme 1&2 - Start 
         !         ID = SeqNinSBox(tSub1,i)
         !         NBOURTotN = NBOURTotN + 1
         !         IDNBOURN(NBOURTotN) = ID
         !         NBOURNRate(NBOURTotN) = rate(ID)
         !         NBOURTotRate =  NBOURTotRate + rate(ID)
         !      ENDDO
         !      IF(tSub2.GT.0)THEN
         !         DO i = 1, NinSBox(tSub2)
         !            ID = SeqNinSBox(tSub2,i)
         !            NBOURTotN = NBOURTotN + 1
         !            IDNBOURN(NBOURTotN) = ID
         !            NBOURNRate(NBOURTotN) = rate(ID)
         !            NBOURTotRate =  NBOURTotRate + rate(ID)     
         !         ENDDO
         !      ENDIF                                             ! Scheme 1&2 - End
               IF(NBOURTotN.GT.0)THEN
                  call random_number(rdn)
                  NBws=rdn*NBOURTotRate
                  NBcumw = NBOURNRate(1)
                  NBiS=1
                  iS = IDNBOURN(NBiS)
                  do while(NBcumw.lt.NBws)
                     NBiS=NBiS+1
                     NBcumw=NBcumw+NBOURNRate(NBiS)
                     iS = IDNBOURN(NBiS)
                  enddo
                  FLAG1 = 1
               ENDIF          
        ENDDO
           
        return
        endsubroutine NBOURSelect
        
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      SUBROUTINE POSITIONCHECK
!      CHECK THE DISTANCE BETWEEN TWO PARTICLES IF < 0.01*SIGAMA
!      SHOULD REJECT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE  PositionCheck(i,FLAG1)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       implicit none
       
       integer i, j,k,m, FLAG1
       
       real*8 xmclj, ymclj, zmclj
       real*8 xlj, ylj, zlj, ljd
       real*8 xcl, ycl, zcl, cld
       real*8 sigma
       
         FLAG1 = 0

         DO j = 1, Npart
            IF(i.NE.j)THEN
                xmclj = abs(mcx(i) - mcx(j))
                ymclj = abs(mcy(i) - mcy(j))
                zmclj = abs(mcz(i) - mcz(j))        
               !!------------
               !! LJ PART
               !!------------                   
                do k=1, NLJ                      ! the LJ sites of particle i
                    do m = 1, NLJ                 ! the LJ sites of particle j 
                       sigma = (sigmaFF(k)+sigmaFF(m))/2.0e0
                       xlj = abs(ljx(j,m) - ljx(i,k))
                       ylj = abs(ljy(j,m) - ljy(i,k))
                       zlj = abs(ljz(j,m) - ljz(i,k))
                       if(pbcx)then
                            if(xmclj > VLEBoxlengthX/2.0D0) xlj = xlj - VLEBoxlengthX
                       endif
                       if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ylj = ylj - VLEBoxlengthY
                       endif
                       if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zlj = zlj - VLEBoxlengthZ
                       endif 
                       ljd = sqrt(xlj*xlj + ylj*ylj + zlj*zlj)
                       if(ljd .LT. 0.01*sigma)then
                          FLAG1 = 1
                          RETURN
                          EXIT
                       endif   
                    enddo
                enddo

               !!-----------------
               !! COULOMB PART
               !!-----------------                   
                  do k = 1, NCL
                     do m = 1, NCL
                        sigma = (sigmaFF(k)+sigmaFF(m))/2.0e0
                        xcl =  abs(CLx(j,m) - CLx(i,k))
                        ycl =  abs(CLy(j,m) - CLy(i,k))
                        zcl =  abs(CLz(j,m) - CLz(i,k))
                        if(pbcx)then      
                           if(xmclj > VLEBoxlengthX/2.0D0) xcl = xcl - VLEBoxlengthX
                        endif
                        if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ycl = ycl - VLEBoxlengthY
                        endif
                        if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zcl = zcl - VLEBoxlengthZ
                        endif
                        cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                        if(cld .LT. 0.01*sigma)then
                           FLAG1 = 1
                           RETURN
                           EXIT
                        endif
                     enddo
                  enddo                
            ENDIF  
         ENDDO

       RETURN
       ENDSUBROUTINE PositionCheck
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE PairECHECK
!      CHECK THE INTERACTION TWO PARTICLES IF >75
!      SHOULD REJECT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE  PairECheck(i,FLAG2)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       implicit none
       
       integer i, j,k,m, FLAG2
       
       real*8 xmclj, ymclj, zmclj
       real*8 xmcn, ymcn,zmcn
       real*8 mcd2, sigma, welldepth
       real*8 xlj, ylj, zlj, ljd2, U,V,U1st
       real*8 xcl, ycl, zcl, cld, CLU,CLV,CLU1st
       real*8 PairE
       
         FLAG2 = 0

         DO j = 1, Npart
            IF(i.NE.j)THEN
                PairE = 0.0
                xmclj = abs(mcx(i) - mcx(j))
                ymclj = abs(mcy(i) - mcy(j)) 
                zmclj = abs(mcz(i) - mcz(j))
                xmcn = xmclj
                ymcn = ymclj
                zmcn = zmclj
                 if(pbcx)then
                       if(xmclj > VLEBoxlengthX/2.0D0) xmcn = xmclj - VLEBoxlengthX
                 endif
                 if(pbcy)then
                       if(ymclj > VLEBoxlengthY/2.0D0) ymcn = ymclj - VLEBoxlengthY
                 endif
                 if(pbcz)then
                      if(zmclj > VLEBoxlengthZ/2.0D0) zmcn = zmclj - VLEBoxlengthZ
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
!                ----------------------------------------------------
!                ENSURE THAT WE DEAL WITH THE NEAREST PERIODIC IMAGES
!                ----------------------------------------------------
                       if(pbcx)then
                            if(xmclj > VLEBoxlengthX/2.0D0) xlj = xlj - VLEBoxlengthX
                       endif
                       if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ylj = ylj - VLEBoxlengthY
                       endif
                       if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zlj = zlj - VLEBoxlengthZ
                       endif 
!                -----------------------
!                INTER-PARTICLE DISTANCE
!                -----------------------
                       ljd2 = xlj*xlj + ylj*ylj + zlj*zlj
!                ----------------
!                POTENTIAL ENERGY
!                ----------------
                      call potentialEnergy(ljd2,sigma,welldepth,U,V,U1st)
                       PairE = PairE + U
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
                        if(pbcx)then      
                           if(xmclj > VLEBoxlengthX/2.0D0) xcl = xcl - VLEBoxlengthX
                        endif
                        if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ycl = ycl - VLEBoxlengthY
                        endif
                        if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zcl = zcl - VLEBoxlengthZ
                        endif
                        cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                        call coulombforce(cld, charges(k),charges(m),CLU,CLV,CLU1st)
                        PairE = PairE + CLU
                     enddo
                  enddo   
                  IF(PairE/Temperature.GT.75.0)THEN
                      FLAG2 = 1
                      RETURN
                  ENDIF
               endif   
            ENDIF  
         ENDDO

       RETURN
       ENDSUBROUTINE PairECheck      
       
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ENSEMBLESUM
!      SUM UP THE PROPERTIES IN SAMPLING STAGE
!      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE ENSEMBLESUM
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       implicit none
       
       integer i
       
       real*8 rdn, TEnergy
       real*8 numberOfSampling, EnergyTotal,FirstDerivETotal, FirstDerivETotal_Ale
       real*8 totNinPropSub(3), VirialPropTotal(3), EnergyPropTotal(3)
       
       common/SAMPSUM1/TEnergy,numberOfSampling, EnergyTotal,FirstDerivETotal, FirstDerivETotal_Ale
       common/SAMPSUM2/totNinPropSub, VirialPropTotal, EnergyPropTotal
       
        
            numberOfSampling  = numberOfSampling + 1.0
            call random_number(rdn)
            InterTime = log(1.0/rdn)/TotRate   
            TotTime = TotTime + InterTime
            EnergyTotal = EnergyTotal + TEnergy*InterTime          
            do i = 1, SubBox
                totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
            enddo 
            do i = 1, 3
                 VirialPropTotal(i) = VirialPropTotal(i) + VirialProp(i)*InterTime
                 EnergyPropTotal(i) = EnergyPropTotal(i) + EnergyProp(i)*InterTime
                 totNinPropSub(i)   = totNinPropSub(i) + NinPropSub(i)*InterTime
             enddo
             FirstDerivETotal = FirstDerivETotal + FirstDerivE*InterTime
             FirstDerivETotal_Ale = FirstDerivETotal_Ale + FirstDerivE_Ale*InterTime
             
             if(WIDOM_VLE)then
                 call WIDOM
                 do i = 1,SubBox
                   ExChemPVLETotal(i) = ExChemPVLETotal(i) + ExcessChemicalPVLE(i)*InterTime
                   IdChemPVLETotal(i) = IdChemPVLETotal(i) + IdealChemicalPVLE(i)*InterTime  
                   ChemPVLETotal(i)   = ChemPVLETotal(i) + ChemicalPVLE(i)*InterTime
                 enddo
              endif
              
       RETURN
       ENDSUBROUTINE ENSEMBLESUM
       
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      SUBROUTINE ConfPOSITIONCHECK
!      CHECK THE DISTANCE BETWEEN TWO PARTICLES IF < 0.01*SIGAMA
!      SHOULD REJECT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE  ConfPositionCheck(FLAG1)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       implicit none
       
       integer i, j,k,m, FLAG1
       
       real*8 xmclj, ymclj, zmclj
       real*8 xlj, ylj, zlj, ljd
       real*8 xcl, ycl, zcl, cld
       real*8 sigma
       
         FLAG1 = 0
         DO i = 1, Npart-1
           DO j = i+1, Npart
                xmclj = abs(mcx(i) - mcx(j))
                ymclj = abs(mcy(i) - mcy(j))
                zmclj = abs(mcz(i) - mcz(j))        
               !!------------
               !! LJ PART
               !!------------                   
                do k=1, NLJ                      ! the LJ sites of particle i
                    do m = 1, NLJ                 ! the LJ sites of particle j 
                       sigma = (sigmaFF(k)+sigmaFF(m))/2.0e0
                       xlj = abs(ljx(j,m) - ljx(i,k))
                       ylj = abs(ljy(j,m) - ljy(i,k))
                       zlj = abs(ljz(j,m) - ljz(i,k))
                       if(pbcx)then
                            if(xmclj > VLEBoxlengthX/2.0D0) xlj = xlj - VLEBoxlengthX
                       endif
                       if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ylj = ylj - VLEBoxlengthY
                       endif
                       if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zlj = zlj - VLEBoxlengthZ
                       endif 
                       ljd = sqrt(xlj*xlj + ylj*ylj + zlj*zlj)
                       if(ljd .LT. 0.01*sigma)then
                          FLAG1 = 1
                          RETURN
                          EXIT
                       endif   
                    enddo
                enddo

               !!-----------------
               !! COULOMB PART
               !!-----------------                   
                  do k = 1, NCL
                     do m = 1, NCL
                        sigma = (sigmaFF(k)+sigmaFF(m))/2.0e0
                        xcl =  abs(CLx(j,m) - CLx(i,k))
                        ycl =  abs(CLy(j,m) - CLy(i,k))
                        zcl =  abs(CLz(j,m) - CLz(i,k))
                        if(pbcx)then      
                           if(xmclj > VLEBoxlengthX/2.0D0) xcl = xcl - VLEBoxlengthX
                        endif
                        if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ycl = ycl - VLEBoxlengthY
                        endif
                        if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zcl = zcl - VLEBoxlengthZ
                        endif
                        cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                        if(cld .LT. 0.01*sigma)then
                           FLAG1 = 1
                           RETURN
                           EXIT
                        endif
                     enddo
                  enddo                
           ENDDO
         ENDDO
       
       RETURN
       ENDSUBROUTINE ConfPositionCheck
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ConfPairECHECK
!      CHECK THE INTERACTION TWO PARTICLES IF >75
!      SHOULD REJECT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE  ConfPairECheck(FLAG2)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       implicit none
       
       integer i, j,k,m, FLAG2
       
       real*8 xmclj, ymclj, zmclj
       real*8 xmcn, ymcn,zmcn
       real*8 mcd2, sigma, welldepth
       real*8 xlj, ylj, zlj, ljd2, U,V,U1st
       real*8 xcl, ycl, zcl, cld, CLU,CLV,CLU1st
       real*8 PairE
       
         FLAG2 = 0

         DO i = 1, Npart-1
            DO j = i+1, Npart
                PairE = 0.0
                xmclj = abs(mcx(i) - mcx(j))
                ymclj = abs(mcy(i) - mcy(j)) 
                zmclj = abs(mcz(i) - mcz(j))
                xmcn = xmclj
                ymcn = ymclj
                zmcn = zmclj
                 if(pbcx)then
                       if(xmclj > VLEBoxlengthX/2.0D0) xmcn = xmclj - VLEBoxlengthX
                 endif
                 if(pbcy)then
                       if(ymclj > VLEBoxlengthY/2.0D0) ymcn = ymclj - VLEBoxlengthY
                 endif
                 if(pbcz)then
                      if(zmclj > VLEBoxlengthZ/2.0D0) zmcn = zmclj - VLEBoxlengthZ
                 endif 
                 mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
                IF(mcd2.gt.r2cutoff)THEN
                   CYCLE
                ELSE
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
!                ----------------------------------------------------
!                ENSURE THAT WE DEAL WITH THE NEAREST PERIODIC IMAGES
!                ----------------------------------------------------
                       if(pbcx)then
                            if(xmclj > VLEBoxlengthX/2.0D0) xlj = xlj - VLEBoxlengthX
                       endif
                       if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ylj = ylj - VLEBoxlengthY
                       endif
                       if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zlj = zlj - VLEBoxlengthZ
                       endif 
!                -----------------------
!                INTER-PARTICLE DISTANCE
!                -----------------------
                       ljd2 = xlj*xlj + ylj*ylj + zlj*zlj
!                ----------------
!                POTENTIAL ENERGY
!                ----------------
                       call potentialEnergy(ljd2,sigma,welldepth,U,V,U1st)
                       PairE = PairE + U
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
                        if(pbcx)then      
                           if(xmclj > VLEBoxlengthX/2.0D0) xcl = xcl - VLEBoxlengthX
                        endif
                        if(pbcy)then
                            if(ymclj > VLEBoxlengthY/2.0D0) ycl = ycl - VLEBoxlengthY
                        endif
                        if(pbcz)then
                            if(zmclj > VLEBoxlengthZ/2.0D0) zcl = zcl - VLEBoxlengthZ
                        endif
                        cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                        call coulombforce(cld, charges(k),charges(m),CLU,CLV,CLU1st)
                        PairE = PairE + CLU
                     enddo
                  enddo 
                  IF(PairE/Temperature.GT.75.0)THEN
                     FLAG2 = 1
                     RETURN
                  ENDIF
               ENDIF   
            ENDDO
         ENDDO
       
       RETURN
       ENDSUBROUTINE ConfPairECheck      

       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE STEELE - POTENTIAL
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE steelepotential(i, Energy)
       USE Constant_M
       USE FLUID_SOLID_M
       USE MCSETTING_M
       USE PHYSICAL_M
       implicit none
    
       integer i,j, k, miniL 
       
       real*8 Energy, iEnergy, rhosperM3, rho_s
       real*8 mcdistanceZ,distanceZ
       real*8 SigmaSF, WellDepthSF
       
              
       rhosperM3 = 114.0E27  
       Energy    = 0.0e0
       rho_s     = rhosperM3*((SCALELENGTH*1.0E-10)**3)
       
       DO k = 1, Nsteele
         miniL = 0
         do j=1,NLJ
            sigmaSF     = (sigmaFF(j) + sigmaSS)/2.0e0
            welldepthSF = sqrt(welldepthFF(j)*welldepthSS) 
            distanceZ = abs(LJz(i,j)-Psteele(k))
            IF(distanceZ .LT. Dist_Limit*sigmaSF)THEN
               Energy = PairEsf_Limit*Temperature
               MiniL = 1
               EXIT
            ELSE
               iEnergy = 4.0E0*pi*WellDepthSF*rho_s*(SigmaSF**2)*GrapheneLayer*((SigmaSF/distanceZ)**10/5.0E0 &
&                   - (SigmaSF/distanceZ)**4/2.0E0 - SigmaSF**4/(6.0E0*GrapheneLayer*(distanceZ+0.61e0*GrapheneLayer)**3))
               Energy = Energy + iEnergy
            ENDIF
            IF(MiniL.EQ.1)EXIT
         enddo    
       ENDDO
       
       IF(Energy/Temperature .GT. PairEsf_Limit)Energy = PairEsf_Limit*Temperature
       
       RETURN
       ENDSUBROUTINE steelepotential

 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ACCESSIBLE VOLUME
!      THIS SUBROUTINE IS USED TO GET THE TOTAL ACCESSIBLE VOLUME 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE accessiblevolume
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       implicit none
       
       integer i,j,m,iInsert, successinsert,iCycle, iiCycle
       integer bflag
       integer Nrotate, Maxrotate
       
       real*8 rdn
       real*8 imaxerror(100000), iAccVolume(100000),error, ierror
       real*8 x1,y1,z1, Energy, EnergyC
       real*8 CLEnergy  
       real*8 LJxnew(20),LJynew(20),LJznew(20) 
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 AccLengthX, AccLengthY, AccLengthZ, AccBulkVolume
       real*8 AccX0, AccY0, AccZ0

       iInsert       = 0
       successinsert = 0
       iCycle        = 0
       error         = 1.0e0
       iiCycle       = 0
       open(unit=22,file='accparticle',status='unknown' )
!       open(unit=20,file='reasonableinsert',status='unknown' )
!       open(unit=19,file='allinsert',status='unknown' )
       IF(SURFACE)THEN
          AccX0 = 0.0E0
          AccY0 = 0.0E0
          AccZ0 = 0.0E0
          AccLengthX = BoxLengthX
          AccLengthY = BoxLengthY
          AccLengthZ = AccVZH
          AccBulkVolume = AccLengthX*AccLengthY*AccLengthZ
       ENDIF
       
       IF(PORE)THEN
          AccX0 = 0.0E0
          AccLengthX = BoxLengthX
          AccY0 = AccVYL
          AccLengthY = AccVYH - AccVYL
          IF(NPartinY.EQ.1)THEN
             AccZ0 = 0.0E0
             AccLengthZ = BoxLengthZ
          ELSE
             AccZ0 = (BoxLengthZ - Part2Lz)/2.0E0
             AccLengthZ = Part2Lz
          ENDIF
          AccBulkVolume = AccLengthX*AccLengthY*AccLengthZ
       ENDIF


       do while (error .gt. 1.0e-3) 
            iCycle = iCycle + 1
            imaxerror(iCycle) = 0.0E0
          do i=1,1000                
             iInsert = iInsert + 1 
           
             call random_number(rdn)
             x1 = rdn * AccLengthX + AccX0 
             call random_number(rdn)
             y1 = rdn * AccLengthY + AccY0
             call random_number(rdn)
             z1 = rdn * AccLengthZ + AccZ0

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
                 
                MCx(1) = x1
                MCy(1) = y1
                MCz(1) = z1
                do m = 1, NLJ
                   LJx(1,m) = LJxnew(m) + MCx(1) - MCx0
                   LJy(1,m) = LJynew(m) + MCy(1) - MCy0
                   LJz(1,m) = LJznew(m) + MCz(1) - MCz0
                enddo
         
                do m =1, NCL
                   CLx(1,m) = CLxnew(m) + MCx(1) - MCx0
                   CLy(1,m) = CLynew(m) + MCy(1) - MCy0
                   CLz(1,m) = CLznew(m) + MCz(1) - MCz0
                enddo               
                call SFEnergySingle(1, EnergyC )
                    if (EnergyC .le. 0.0)then 
                        successinsert = successinsert + 1
                        EXIT
                    endif
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
       
       write(62,*)'iInsert, successinsert',iInsert, successinsert
!       write(22,*)'iInsert, successinsert',iInsert, successinsert
!       write(22,*)'friction',dble(successinsert)/dble(iInsert)
       
       AccVolume = dble(successinsert)/dble(iInsert)*AccBulkVolume
!       write(22,*)'BulkVolume, AccVolume',BulkVolume, AccVolume
     
       RETURN
       ENDSUBROUTINE ACCESSIBLEVOLUME
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE InserSurf
!    INSERT THE NEW ADDED PARTICLES IN THE SIMULATION BOX FOR THE CASE OF
!    SURFACE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE InsertSurf(LastN,AimN)   
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer LastN, AimN
      integer i,j, flag,k
      real*8 rdn, distance2,distance
      
        i=1       
        do while(i .le. AimN)
           call random_number(rdn)
           MCx(LastN+i) = rdn*BoxLengthX
           call random_number(rdn)
           MCy(LastN+i) = rdn*BoxLengthY
           call random_number(rdn)
           IF(MassTSFER)THEN
              !MCz(LastN+i) = rdn*(LayerZH(SubBox)-LayerZL(SubBox)) + LayerZL(SubBox)
              MCz(LastN+i) = rdn*(InsertZH-InsertZL) + InsertZL
           ELSE
              MCz(LastN+i) = rdn*BoxLengthZ
           ENDIF
          
           flag=0 
           
           IF(Steele)THEN
              do k=1,NSteele
                 if( abs(MCz(LastN+i)-PSteele(k)) .LT. 0.75e0 )then
                   flag=1
                   EXIT
                 endif
              enddo
           ENDIF
           if(flag.eq.1)CYCLE

            
           do j=1,LastN+i-1
              distance2 = (MCx(LastN+i)-MCx(j))**2.0 + (MCy(LastN+i)-MCy(j))**2.0 + (MCz(LastN+i)-MCz(j))**2.0
              distance = sqrt(distance2)
              if(distance.LT.0.5E0)then
                 flag=1
                 EXIT
              endif
          enddo
          if(flag.eq.1)CYCLE
                     
              flag=0
             do j= 1, SubBox
                if(MCz(LastN+i).GE. LayerZL(j) .AND. MCz(LastN+i).LT. LayerZH(j) )then
                   NinSBox(j)= NinSBox(j)+1
                   SeqNinSBox(j,NinSBox(j))=LastN+i
                   indNinSBox(LastN+i)=j
                   flag=1
                endif
                if(flag.eq.1)EXIT
             enddo

            do j = 1, NLJ
             LJx(LastN+i,j) = LJx0(j) + MCx(LastN+i) - MCx0
             LJy(LastN+i,j) = LJy0(j) + MCy(LastN+i) - MCy0
             LJz(LastN+i,j) = LJz0(j) + MCz(LastN+i) - MCz0
            enddo
         
           do j =1, NCL
             CLx(LastN+i,j) = CLx0(j) + MCx(LastN+i) - MCx0
             CLy(LastN+i,j) = CLy0(j) + MCy(LastN+i) - MCy0
             CLz(LastN+i,j) = CLz0(j) + MCz(LastN+i) - MCz0
           enddo
           
           i=i+1

      enddo
      
     return
     END SUBROUTINE InsertSurf
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE InserPore
!    INSERT THE NEW ADDED PARTICLES IN THE SIMULATION BOX FOR THE CASE OF
!    PORE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE InsertPore(LastN,AimN)   
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer LastN, AimN
      integer i,j, flag,k
      real*8 rdn, distance2,distance
      
        i=1       
        do while(i .le. AimN)
           call random_number(rdn)
           MCx(LastN+i) = rdn*BoxLengthX
           call random_number(rdn)
           MCz(LastN+i) = rdn*BoxLengthZ
           call random_number(rdn)
           IF(MassTSFER .OR. MeanFPath)THEN
              MCy(LastN+i) = rdn*Part1Ly
           ELSEIF( (NPartinY.EQ.3) .AND. (i.GT.int(AimN/2)) )THEN
              MCy(LastN+i) = rdn*Part3Ly + Part1Ly + Part2Ly
           ELSE
              MCy(LastN+i) = rdn*Part1Ly
           ENDIF
          
           flag=0 
           do j=1,LastN+i-1
              distance2 = (MCx(LastN+i)-MCx(j))**2.0 + (MCy(LastN+i)-MCy(j))**2.0 + (MCz(LastN+i)-MCz(j))**2.0
              distance = sqrt(distance2)
              if(distance.LT.0.5E0)then
                 flag=1
                 EXIT
              endif
           enddo
           if(flag.eq.1)CYCLE
                     
              flag=0
             do j= 1, SubBox
                if(MCy(LastN+i).GE. LayerYL(j) .AND. MCy(LastN+i).LT. LayerYH(j) )then
                   NinSBox(j)= NinSBox(j)+1
                   SeqNinSBox(j,NinSBox(j))=LastN+i
                   indNinSBox(LastN+i)=j
                   flag=1
                endif
                if(flag.eq.1)EXIT
             enddo

            do j = 1, NLJ
             LJx(LastN+i,j) = LJx0(j) + MCx(LastN+i) - MCx0
             LJy(LastN+i,j) = LJy0(j) + MCy(LastN+i) - MCy0
             LJz(LastN+i,j) = LJz0(j) + MCz(LastN+i) - MCz0
            enddo
         
           do j =1, NCL
             CLx(LastN+i,j) = CLx0(j) + MCx(LastN+i) - MCx0
             CLy(LastN+i,j) = CLy0(j) + MCy(LastN+i) - MCy0
             CLz(LastN+i,j) = CLz0(j) + MCz(LastN+i) - MCz0
           enddo
           
           i=i+1

      enddo
      
     return
     END SUBROUTINE InsertPore     
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE SFENERGYSINGLE
!      THIS SUBROUTINE CALCULATES THE ENERGY OF ONE PARTICLE WITH THE SOLID
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SFEnergySingle(i, EnergyC)
      USE FLUID_SOLID_M
      USE MCSETTING_M
      USE PHYSICAL_M
      IMPLICIT NONE
      
      integer i, j
      
      real*8  EnergyC, UC
      real*8  BU, PEnergyC
      
         EnergyC  = 0.0E0
         if(steele)then
            call steelepotential(i,UC)
            EnergyC = EnergyC + UC 
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
          endif

               
      RETURN
      ENDSUBROUTINE SFEnergySingle
      
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Bojan Potential
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     SUBROUTINE BojanPotential(bi,i,BU) 
     USE Constant_M
     USE FLUID_SOLID_M
     USE MCSETTING_M
     USE PHYSICAL_M
     USE  POREFIGURE_M
     implicit none
     
     integer bi, i, j, k, MiniL
     
     real*8 y,z,BU,Cz(20),CY
     real*8 totforce
     real*8 BWellDepthSF
     real*8 layerTotal
     real*8 Inctotforce, DisP_L, DisP_CZ, DisCY
     real*8 Clo1totforce,Clo2totforce


       BU = 0.0E0    
       Cz = 0.0e0
       
       IF(bi.le.NS)THEN
          do j = 1, StriplayerN(bi)
             layerTotal = StriplayerN(bi)
             
             if(TopStrip(bi))then
                 Cz(j) = StripZ0(bi)+dble(j-1)*Stripgap(bi)
             else
                 Cz(j) = StripZ0(bi)-dble(j-1)*Stripgap(bi)
             endif
             
             if(BojanSlit)then
                layerTotal = StriplayerN(bi)*2
                if(StripZ0(bi).gt. 0.0e0)then
                   Cz(j+StriplayerN(bi)) = Cz(j)+BoxlengthZ-2.0e0*StripZ0(bi)+dble(StriplayerN(bi)-1)*Stripgap(bi)
                else
                   Cz(j+StriplayerN(bi)) = Cz(j)+BoxlengthZ+dble(StriplayerN(bi)-1)*Stripgap(bi)
                endif
             endif
           enddo
        ENDIF
              
       DO k = 1, NLJ
           IF(bi.le.NS)THEN
              totforce   = 0.0e0   
              BsigmaSF = (sigmaFF(k)+BsigmaSS(bi))/2.0e0
              BwelldepthSF = sqrt(welldepthFF(k)*BwelldepthSS(bi))
              Cy =  LJy(i,k) - StripC(bi)          
              Py =  StripW(bi)/2.0e0 - Cy
              Ny = -StripW(bi)/2.0e0 - Cy
          
              do j = 1, layerTotal
                 MiniL = 0
                 BdeltaZ = abs(Cz(j)-LJz(i,k))
                 IF(BdeltaZ.LT. Dist_Limit*BsigmaSF)THEN
                     BU = PairEsf_Limit*Temperature
                     MiniL = 1
                     EXIT
                 ELSE
                    call BojanCalculation          
                    totforce = totforce + subforce
                 ENDIF
              enddo
       
              IF(MiniL.EQ.1)THEN
                 RETURN
              ELSE
                 BU = BU + 2.0e0*pi*rhosperM2(bi)*(BsigmaSF**2.0e0)*BwelldepthSF*totforce
              ENDIF
           ENDIF  
           
          IF(bi.eq.NS+2)then   ! For Close1End .OR. Close2Ends
             IF(Close1End)THEN
                 Clo1totforce = 0.0e0
                 IF(Clo1Hard)Then
                    Clo1totforce = 0.0e0
                 ELSE
                    BsigmaSF = (sigmaFF(k)+Clo1sigmaSS)/2.0e0
                    BwelldepthSF = sqrt(welldepthFF(k)*Clo1welldepthSS)   
                    Cy =  LJZ(i,k)- Clo1Mid_z
                    Py =  Clo1LengthZ/2.0e0 - Cy
                    Ny = -Clo1LengthZ/2.0e0 - Cy
                    do j = 1, Clo1LayerN
                       BdeltaZ = LJy(i,k)+(j-1)*Clo1gap
                       call BojanCalculation
                       Clo1totforce = Clo1totforce + subforce
                    enddo
                 ENDIF
                 BU = BU+2.0e0*pi*Clo1rhosperM2*(BsigmaSF**2.0e0)*BwelldepthSF*Clo1totforce  
             ENDIF
             IF(Close2Ends)THEN
                 Clo2totforce = 0.0e0
                 BsigmaSF = (sigmaFF(k)+Clo2sigmaSS)/2.0e0
                 BwelldepthSF = sqrt(welldepthFF(k)*Clo2welldepthSS)   
                 Cy =  LJZ(i,k)- Clo2Mid_z
                 Py =  Clo2LengthZ/2.0e0 - Cy
                 Ny = -Clo2LengthZ/2.0e0 - Cy
                 do j = 1, Clo2LayerN
                    BdeltaZ = BoxlengthY-LJy(i,k)+(j-1)*Clo2gap
                    call BojanCalculation
                    Clo2totforce = Clo2totforce + subforce
                 enddo
                 BU = BU + 2.0e0*pi*Clo2rhosperM2*(BsigmaSF**2.0e0)*BwelldepthSF*Clo2totforce
             ENDIF
         ENDIF
                   
       ENDDO
              
       IF(BU/Temperature .GT. PairEsf_Limit)BU = PairEsf_Limit*Temperature
 
    RETURN
    ENDSUBROUTINE BojanPotential     
    
!+++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Bojan Calculation
!
!+++++++++++++++++++++++++++++++++++++              
       SUBROUTINE BojanCalculation
       USE FLUID_SOLID_M
       implicit none
                
       real*8 Prepulsion,Pattration,Nrepulsion,Nattration
       real*8 subforce1
       real*8 sub1,sub2,sub3,sub4,sub5,sigmadNy,sigmadPy, zdsigma2
     
        IF( (BdeltaZ .lt. 0.1e0) .and. (Py*Ny .gt. 0.0e0) )THEN
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
          ELSE           
             call Brepatt(Py,BdeltaZ,BsigmaSF,Prepulsion,Pattration)
             call Brepatt(Ny,BdeltaZ,BsigmaSF,Nrepulsion,Nattration)
             subforce  = (Prepulsion - Nrepulsion)- (Pattration - Nattration)
          ENDIF
               
        RETURN
        ENDSUBROUTINE BojanCalculation
    
!+++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Bojan Repulsion and attraction
!
!+++++++++++++++++++++++++++++++++++++++++++++++++             
       SUBROUTINE Brepatt (y, z, BsigmaSF, Repulsion, Attraction)
       implicit none
       
       real*8 y,z,z2,sumy2z2
       real*8 Repulsion, Attraction
       real*8 BsigmaSF
       
       
            z2 = z*z
            sumy2z2 = y*y + z2
            
            Repulsion = y*(BsigmaSF**10.0e0)* ( 0.2e0/(z2**4.0e0) + 0.1e0/(z2**3.0e0)/sumy2z2    &
 &                      + 3.0e0/4.0e1/((z2*sumy2z2)**2.0e0) + 1.0e0/1.6e1/z2/(sumy2z2**3.0e0)  &
 &                      + 7.0e0/1.28e2/(sumy2z2**4.0e0) )/Z2/sqrt(sumy2z2)
            
            Attraction = y*(BsigmaSF**4.0e0)*(0.5e0/z2 + 0.25e0/sumy2z2)/Z2/sqrt(sumy2z2)
            
            
       RETURN     
       ENDSUBROUTINE Brepatt
   
   
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE MC_PORE
!    FOR THE Kinetic MC for adsorpiton in PORE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       subroutine MC_PORE
       USE Constant_M
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE SUBBOX_M
       USE ROTATE_M
       USE POREFIGURE_M
       implicit none
    
       integer NumberOfMove, NumberOfAcceptanceMove
       integer iEqui, iMove, iCycle
       integer i, j, fSub, cutoffopt
       integer numberofequilibrium
       integer iSub, iS, iRot,iSRot
       integer FLAG1, FLAG2, FLAG3
       integer Nbranch, ibranch, iSOld, OverZ
    
       real*8 secnds, time,rdn
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
       real*8 Nb_rho, Nb_rhoMolPerM3, Ex_Nb_rho, Ex_Nb_rhoMolPerM3 
       real*8 Ave_MTtime, TotMTtime,InterTimeOld,PresentTime
                    
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

          BulkVolume = 0.0E0
          DO i = 1, SubBox
             BulkVolume  = BulkVolume + SubVol(i)
          ENDDO
          
          IF(NpartinY.EQ.2)THEN
            RatioV1 = (BoxLengthX*BoxLengthZ*Part1Ly)/BulkVolume
          ELSEIF(NpartinY.EQ.3)THEN
            RatioV1 = (BoxLengthX*BoxLengthZ*Part1Ly)/BulkVolume
            RatioV2 = RatioV1 + (BoxLengthX*Part2Lz*Part2Ly)/BulkVolume
          ENDIF 

          IF(ExtraBox)THEN
             RatioV1 = SubVol(SubBox)/BulkVolume
          ENDIF

          IF( (NPartinY.EQ.1).AND.(Steele.OR.Bojan) )THEN
             TotSurfaceArea = 2.0*BoxlengthX*BoxlengthY
          ELSE
             TotSurfaceArea = 2.0*BoxlengthX*Part2Ly
          ENDIF
         
         if(Close1End.AND.(.NOT. Clo1Hard))then
              TotSurfaceArea = TotSurfaceArea+opClo1LengthZ*BoxlengthX
         endif
         
         if(Close2Ends)then
             TotSurfaceArea = TotSurfaceArea+opClo2lengthZ*BoxlengthX
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
       

        open(unit=10, file='PoreIsotherm.txt',status='unknown')
        rewind(10)
        write(10,*)'------------------------------------------------------------------------------------------------'
        write(10,'(A4,10A16)') 'no.','rho','P','Npart','Abs_rhoinPore','Exc_rhoinPore','SurfE','N-Nbulk','U-Ubulk','Ave_NinBin'
        write(10,'(1x,10A16)')'(-)','(mol/m3)','(pa)','(-)','(mol/m3)','(mol/m3)','(umol/m2)','(-)','(Joule)','(-)'
        write(10,*)'-------------------------------------------------------------------------------------------------'
!       ------------------------
!       STARTING CONFIGURATION
!       ------------------------
        ALLOCATE(rate(5000), AccuR(5000), Energyi(5000))
        ALLOCATE(EnergyMatrix(5000,5000),VirialMatrix(5000,5000))
        Npart = Npart  !+ int( SEinput(1)*TotSurfaceArea)+1
        NinSBox = 0
        SeqNinSBox = 0
        IndNinSBox = 0
        
        CALL InsertPore(0,Npart)
        
        IF(LocalD2D) call DDistribution_2D(0)
        Nbranch = 1
        if(NVTDes)Nbranch = 2
     DO ibranch =1, Nbranch
        DO iden = 1, NumOfrho
            write(62,*)'Inserting new particles'
            if(ibranch.eq.1)then
               CALL InsertPore(Npart,AddN)
               Npart = Npart + AddN 
            elseif(ibranch.eq.2)then
               CALL DelPore
            endif   
            write(62,*)'Npart=',Npart
            write(51,*)'Npart=',Npart
            IF(MassTSFER)THEN
              UpFlux = 0.0
              DownFlux = 0.0
              TotMTtime  = 0.0
              NetFlux = 0.0
            ENDIF
            EnergyMatrix = 0.0E0
            VirialMatrix = 0.0E0
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
            write(62,*)'Calculating ConfigEnergy'
            call ConfigEnergyADS(Energy,CLEnergy, EnergyC)
            TEnergy = Energy + EnergyC
            write(62,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
            write(51,*)'TEnergy, Energy, CLEnergy, EnergyC', TEnergy, Energy, CLEnergy, EnergyC
            call random_number(rdn)
            InterTimeOld = log(1.0/rdn)/TotRate 

            iSOld = 0
            OverZ = -1
            DO iEqui = 1, Ncycle
                if (mod(iEqui,100000).eq.0)then
                   write(62,*)'====================kMCPore EQULIBRATION==================='
                   write(51,*)'====================kMCPore EQULIBRATION==================='
                   write(62,646)iEqui, EnergyTotal/TotTime/dble(Npart),Npart
                   write(51,646)iEqui, EnergyTotal/TotTime/dble(Npart),Npart
      646          format(1x,'CYCLES:',I12, '  AveEnergy:',E15.6, '   Npart:',I5)
                endif
                
!               -----------------------------
!               ADJUST THE SIMULATION BOX
!               -----------------------------         
                if(mod(iEqui,100).eq.0)then
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
!               ---------------------------
!               Randomly pick one particle
!               ---------------------------
                if((iSOld.NE.0).and. (OverZ.EQ.1))then
                       iS = iSOld
                   else
                       call select(iS)
                endif

                CALL TherPropChangeADS(1,iS)
                        
!               -------------------------
!               Pick a position randomly 
!               -------------------------
                IF(MeanFPath)THEN
                   call PoreMFPMCpick(iS,OverZ)
                   IF(OverZ.EQ.1)iSOld = iS
                   IF(OverZ.EQ.-1)iSOld = 0
                ELSE
                   call MassCpick(iS)    
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
                        CALL PairEnergyADS(iS)
                        TEnergy = TEnergy + Energyi(iS) 
                        CALL TherPropChangeADS(2,iS)
            
                call random_number(rdn)
                InterTime = log(1.0/rdn)/TotRate   
                TotTime = TotTime + InterTime
                EnergyTotal = EnergyTotal + TEnergy*InterTime
                IF(MassTSFER)THEN
                      IF(fSub.NE.iSub)CALL FLUXCALCULATE(fSub, iSub)
                      Ave_MTtime = (InterTimeOld+InterTime)/2.0
                      TotMTtime = TotMTtime + Ave_MTtime
                      InterTimeOld = InterTime
                ENDIF       
                do i = 1, SubBox
                     totNinSBox(i) = totNinSBox(i) + NinSBox(i)*InterTime
                enddo 
                    IF(MassTSFER .AND. mod(iEqui,100).EQ.0)THEN
                       do i =1,NFluxSurf
                          IF(TotMTtime.ne.0)NetFlux(i) = (UpFlux(i)-DownFlux(i))/TotMTtime
                       enddo
                       do i = 1, SubBox
                          ave_NinSBox(i) = totNinSBox(i)/TotTime
                       enddo 
                       PresentTime = PresentTime + TotTime/2.0e0
                       write(11,'(I10,100E16.6)')iEqui,PresentTime,(ave_NinSBox(i)/SubVol(i),i=1,SubBox) !,(NetFlux(i),i=1,NFluxSurf) 
                       PresentTime = PresentTime + TotTime/2.0e0
                       DownFlux = 0.0
                       UpFlux   = 0.0
                       TotMTtime  = 0e0
                       totNinSBox = 0.0
                       TotTime = 0.0
                       NetFlux = 0.0
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

     
            EnergyTotal            = 0.0E0
            numberOfSampling       = 0.0E0    
            BulkVirialTotal        = 0.0E0
            BulkEnergyTotal        = 0.0E0
            ExChemPVLETotal        = 0.0
            IdChemPVLETotal        = 0.0 
            ChemPVLETotal          = 0.0
            TotTime                = 0.0
  
             IF(LocalD2D) call DDistribution_2D(1)
        
             totNinSBox     = 0.0
             Nmove=SubBox*(SubBox-1)/2
         
             iSOld = 0
             OverZ = -1
             do iCycle = 1, Ncycles
                  if (mod(iCycle,100000).eq.0)then
                     write(62,*)'================kMCPore SAMPLING==================='
                     write(51,*)'================kMCPore SAMPLING==================='
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
                     ENDIF
                  endif           
!                 ----------------------------------
!                 Pick one particle with Rosenbluth
!                 ----------------------------------
                  if((iSOld.NE.0).and. (OverZ.EQ.1))then
                         iS = iSOld
                     else
                         call select(iS)
                  endif
                  CALL TherPropChangeADS(1,iS)
              
!                   -------------------------
!                   Pick a position randomly 
!                   -------------------------
                    IF(MeanFPath)THEN
                       call PoreMFPMCpick(iS,OverZ)
                       IF(OverZ.EQ.1)iSOld = iS
                       IF(OverZ.EQ.-1)iSOld = 0
                    ELSE
                       call MassCpick(iS)    
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
                           rateTSBox(i)  = rateTSBox(i) + rateSBox(i)*InterTime
                           timeSBox(i)   = timeSBox(i) + InterTime
                        enddo 
                        
 !                 IF(VLEradius)call VLEradiusdistribution(2)
 !                 IF(VLEdistance)call VLEdistancedistribution(2)   
                   IF(mod(iCycle,1000).eq.0)THEN
                      IF(LocalD2D) call DDistribution_2D(2)
                   ENDIF
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
       
           IdChemPkVLE = Temperature*LOG(deBro3/(BoxlengthX*BoxlengthY*BoxlengthZ))
           ExChemPkVLE = Temperature*LOG(numberOfSampling/TotTime)
           ChemPkVLE = IdChemPkVLE + ExChemPkVLE
           ChemPkVLEJPerMole = ChemPkVLE * SCALEENERGY*Rg
           do i=1, SubBox
              ChemPSBox(i) = Temperature*LOG(deBro3*rateTSBox(i)/timeSBox(i)/SubVol(i)) * SCALEENERGY*Rg
           enddo

           
           rho(iden) = rho_SBox(SubFGas)
           P(iden) = rho(iden)*Temperature + BulkVirialAve/SubVol(SubFGas)
           CNBulk = rho(iden)*AccVolume + (NPartinY-1)*ave_NinSBox(SubFGas)
           CEBulk = EinBulkAve*AccVolume/SubVol(SubFGas) + (NPartinY-1)*EinBulkAve
           CEBulkJoule = CEBulk*SCALEENERGY*kB
           call ZLRC(rho(iden), EnergyLRC, VirialLRC)
           PLRC = rho(iden)*VirialLRC
           P(iden) = P(iden) + PLRC
           
           SEinput(iden) = (Npart-CNBulk)/totSurfaceArea
           IF(NPartinY.EQ.1)THEN
              if(ExtraBox)then
                 Nb_rho = dble(Npart-ave_NinSBox(SubFGas))/AccVolume
                 Ex_Nb_rho = Nb_rho - rho(iden)
              else
                 SEinput(iden) = dble(Npart)/totSurfaceArea
                 Nb_rho = dble(Npart)/AccVolume
              endif
           ELSEIF(NPartinY.EQ.2)THEN
              SEinput(iden) = (Npart-ave_NinSBox(SubFGas))/totSurfaceArea
              Nb_rho = (Npart-ave_NinSBox(SubFGas))/AccVolume
              Ex_Nb_rho = Nb_rho - rho(iden)
           ELSEIF(NPartinY.EQ.3)THEN
              SEinput(iden) = (Npart-ave_NinSBox(SubFGas)-ave_NinSBox(SubBox))/totSurfaceArea
              Nb_rho = (Npart-ave_NinSBox(SubFGas)-ave_NinSBox(SubBox))/AccVolume
              Ex_Nb_rho = Nb_rho  - rho(iden)
           ENDIF   
        
!          ----------------------
!          DIMENSIONAL QUANTITIES
!          ----------------------
           rhoMolPerM3              = rho(iden)/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
           SEinputuMolPerM2         = SEinput(iden)/(SCALELENGTH*1.0E-10)**2/AvogadroNumber*1.0e6
           PressurePa               = P(iden)*SCALEENERGY*kB/(SCALELENGTH*1E-10)**3
           EnergyAverageJoule       = EnergyAverage*SCALEENERGY*kB 
           Nb_rhoMolPerM3           = Nb_rho/(SCALELENGTH*1.0E-10)**3/AvogadroNumber  
           Ex_Nb_rhoMolPerM3        = Ex_Nb_rho/(SCALELENGTH*1.0E-10)**3/AvogadroNumber
     
 !          IF(VLEradius)call VLEradiusdistribution(3)
 !          IF(VLEdistance)call VLEdistancedistribution(3)
            IF(LocalD2D)call DDistribution_2D(3) 
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
           &      1x, 'Ns_Density (Mol/m2) ------------------------ : ', e20.10,/, &
           &      1x, 'Pressure (Pa) ------------------------------ : ', e20.10,/, &
           &      1x, 'AccessVolume (m3) -------------------------- : ', e20.10,/, &
           &      1x, 'Average Energy per Particle (J) ------------ : ', e20.10,/, &
           &      1x, 'Time (sec) --------------------------------- : ', f20.10/)
!            -------------------------------
!            STORE THE POSITION OF PARTICLES
!            -------------------------------
             write(52,*) 'Npart=', Npart
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
           write(10,716)iden,rhoMolPerM3,PressurePa, Npart, Nb_rhoMolPerM3, Ex_Nb_rhoMolPerM3,SEinputuMolPerM2,Npart-CNBulk, EnergyAverageJoule-CEBulkJoule, (ave_NinSBox(i), i=1,SubBox), (ChemPSBox(i), i=1,SubBox)
716        FORMAT(1x,I4,2E18.8,I10,200E18.8)  
           
flush(51)
flush(52)
flush(10)

           ENDDO !iden
     ENDDO !ibranch
             
        DEALLOCATE(rate, AccuR, Energyi)
        DEALLOCATE(EnergyMatrix,VirialMatrix)
        IF(LocalD2D)call DDistribution_2D(4) 
        
     RETURN
     END SUBROUTINE MC_PORE          
     
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE MassCpick
!    Pick the mass center for a particle in PORE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
     SUBROUTINE MassCpick(iS)
     USE FLUID_SOLID_M
     USE MCSETTING_M
     IMPLICIT NONE
     
     integer iS
     real*8 rdn, rdn1
                   
                   call random_number(rdn)
                   MCx(iS) = rdn*BoxlengthX 
                   call random_number(rdn)
                   IF(NPartinY.EQ.1)THEN
                       MCy(iS) = rdn*BoxlengthY
                       IF(ExtraBox)THEN
                          call random_number(rdn)
                          if(rdn.lt.RatioV1)then
                              call random_number(rdn1)
                              MCz(iS) = rdn1*ExtraBoxH + BoxlengthZ 
                           else
                              call random_number(rdn1)
                              MCz(iS) = rdn1*BoxlengthZ 
                           endif                           
                       ELSE
                          call random_number(rdn1)
                          MCz(iS) = rdn1*BoxlengthZ 
                       ENDIF 
                   ELSEIF(NPartinY.EQ.2)THEN
                       if(rdn.lt.RatioV1)then
                           call random_number(rdn1)
                           MCy(iS) = rdn1*Part1Ly
                           call random_number(rdn1)
                           MCz(iS) = rdn1*BoxlengthZ 
                       else
                           call random_number(rdn1)
                           MCy(iS) = rdn1*Part2Ly + Part1Ly
                           call random_number(rdn1)
                           MCz(iS) = rdn1*Part2Lz + (BoxLengthZ-Part2Lz)/2.0E0
                       endif
                   ELSEIF(NPartinY.EQ.3)THEN 
                       if(rdn.lt.RatioV1)then
                           call random_number(rdn1)
                           MCy(iS) = rdn1*Part1Ly
                           call random_number(rdn1)
                           MCz(iS) = rdn1*BoxlengthZ 
                       elseif(rdn.ge.RatioV1 .AND. rdn.le.RatioV2)then
                           call random_number(rdn1)
                           MCy(iS) = rdn1*Part2Ly + Part1Ly
                           call random_number(rdn1)
                           MCz(iS) = rdn1*Part2Lz + (BoxLengthZ-Part2Lz)/2.0E0
                       else
                           call random_number(rdn1)
                           MCy(iS) = rdn1*Part3Ly + Part1Ly + Part2Ly
                           call random_number(rdn1)
                           MCz(iS) = rdn1*BoxLengthZ
                       endif
                   ENDIF    
     RETURN
     END SUBROUTINE MassCPick     
     
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE SurfMCpick
!    Pick the mass center for a particle in surface
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
     SUBROUTINE SurfMCpick(iS)
     USE FLUID_SOLID_M
     USE MCSETTING_M
     USE SUBBOX_M
     IMPLICIT NONE
     
     integer iS, fSub, tSub, tSub1, tSub2
     
     real*8 rdn, rdn1, totSubVol
     real*8 SRatioV1, SRatioV2
     
       IF(NBOURBIN)THEN
             fSub = indNinSBox(iS)
             IF(fSub.eq.1)THEN
                tSub1 = fSub + 1
                totSubVol = SubVol(fSub) + SubVol(tSub1) 
             ELSEIF(fSub.eq.SubBox)THEN
                tSub1 = fSub - 1
                totSubVol = SubVol(fSub) + SubVol(tSub1)
             ELSE
                tSub1 = fSub - 1
                tSub2 = fSub + 1
                totSubVol = SubVol(fSub) + SubVol(tSub1) + SubVol(tSub2)
             ENDIF
             
             SRatioV1 = SubVol(fSub)/totSubVol
             SRatioV2 = (SubVol(fSub) + SubVol(tSub1))/totSubVol
             call random_number(rdn)
             MCx(iS) = rdn*BoxlengthX
             call random_number(rdn)
             MCy(iS) = rdn*BoxlengthY
             call random_number(rdn)
             IF(rdn.le.SRatioV1)THEN
                tSub = fSub
             ELSEIF(rdn.gt.SRatioV1 .AND. rdn.le.SRatioV2)THEN
                tSub = tSub1
             ELSE
                tSub = tSub2
             ENDIF
             call random_number(rdn)
             MCz(iS) = rdn*(LayerZH(tSub)-LayerZL(tSub)) + LayerZL(tSub) 
        ELSE
                   call random_number(rdn)
                   MCx(iS) = rdn*BoxlengthX
                   call random_number(rdn)
                   MCy(iS) = rdn*BoxlengthY
                   call random_number(rdn)
                   MCz(iS) = rdn*BoxlengthZ 
       ENDIF  
     RETURN
     END SUBROUTINE SurfMCPick     
     
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE DelPore
!    Delete PARTICLES FROM SIMULATION BOX - DESORPTION FOR PORE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE DelPore   
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer idel,DelN
      integer DelN1, DelN2, DelBox
      
      real*8 rdn
      
      DelN = addN
      IF(NpartinY.EQ.2)THEN
          DelBox = 1
          if(NinSBox(DelBox).LT.addN)DelN=NinSBox(DelBox)
          if(DelN.LT.1)then
             write(10,*)'Oops, the number of partics in gas bin is 0!'
             STOP
          endif 
          CALL DelParticle(DelBox, DelN)    
      ENDIF
      
      IF(NpartinY.EQ.3)THEN
         if( (NinSBox(1)+NinSBox(SubBox)).lt. addN)DelN = NinSBox(1)+NinSBox(SubBox)
         if(DelN.LT.1)then
            write(10,*)'Oops, the number of partics in gas bin is 0!'
            STOP
         endif 
         if(DelN.LT.addN)then
            DelN1 = NinSBox(1)
            DelN2 = NinSBox(SubBox)
         else
            if(NinSBox(1).ge.DelN)then
               DelN1 = DelN
               DelN2 = 0
            else
               DelN1 = NinSBox(1)
               DelN2 = DelN - DelN1 
            endif
         endif
         
           CALL DelParticle(1,DelN1)
           CALL DelParticle(SubBox,DelN2)  
      ENDIF
      
     return
     END SUBROUTINE DelPore
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE DelParticle
!    Delete PARTICLES FROM SELECTED BINS
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE DelParticle(DelBox,DelN)   
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer idel,i,j,k, DelN,ID
      integer aimSub, DelBox
      
      real*8 rdn
      
         idel=1       
         do while(idel .le. DelN)
             call random_number(rdn)
             ID  = int(rdn*NinSBox(DelBox)) + 1
             if(ID.gt.NinSBox(DelBox)) ID = NinSBox(DelBox)
             i=SeqNinSBox(DelBox,ID)        
             SeqNinSBox(DelBox,ID)= SeqNinSBox(DelBox,NinSBox(DelBox))
             NinSBox(DelBox) = NinSBox(DelBox)-1  
             aimSub = IndNinSBox(Npart)
             do j=1, NinSBox(aimSub)
                if(SeqNinSBox(aimSub,j).eq.Npart)then
                   SeqNinSBox(aimSub,j)=i
                   EXIT
                endif
              enddo
              IndNinSBox(i)=aimSub
              mcx(i)=mcx(Npart)
              mcy(i)=mcy(Npart)
              mcz(i)=mcz(Npart)
              do j = 1, NLJ
                 LJx(i,j) = LJx(Npart,j)
                 LJy(i,j) = LJy(Npart,j)
                 LJz(i,j) = LJz(Npart,j)
              enddo
              do j = 1, NCL
                 CLx(i,j) = CLx(Npart,j)
                 CLy(i,j) = CLy(Npart,j)
                 CLz(i,j) = CLz(Npart,j)
              enddo            
              Npart=Npart-1          
              idel=idel+1
         enddo
      
     return
     END SUBROUTINE DelParticle
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE FLUXCALCULATE
!    Delete PARTICLES FROM SELECTED BINS
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE FLUXCALCULATE(fSub, tSub)
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
       integer i,fSub, tSub
       
       IF(fSub .LT. tSub)THEN
          do i = fSub, tSub-1
             UpFlux(i)= UpFlux(i)+1
          enddo
       ELSEIF(fSub .GT. tSub)THEN
          do i = fSub-1, tSub,-1
             DownFlux(i) = DownFlux(i)+1
          enddo       
       ENDIF

     
     
      RETURN
      ENDSUBROUTINE FLUXCALCULATE
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE EntropyCalculation
!    For calculation the entropy of the system
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE EntropyCalculation(switch,TEnergy)
      USE Constant_M
      USE PHYSICAL_M
      USE MCSETTING_M
      USE SUBBOX_M
      implicit none
      
      integer i,switch
      
      real*8  TEnergy, deltaUBin, EntroTotTime
      
      common/Etime/EntroTotTime
      
         IF(switch.eq.1)THEN
           deltaUBin = -0.001E0
           NumUBin = int(abs(LowUsys-HighUsys)/abs(deltaUBin))+1
 !          deltaUBin = (HighUsys-LowUsys)/dble(NumUBin)
           ALLOCATE(UBinL(NumUBin), UBinH(NumUBin), UBinTime(NumUBin),ProbaUBin(NumUBin))
           UBinL(1) = LowUsys
           UBinH(1) = UBinL(1) + deltaUBin
           do i =2, NumUBin
              UBinL(i) = UBinH(i-1)
              UBinH(i) = UBinL(i) + deltaUBin
           enddo
           UBinTime = 0.0
           EntroTotTime = 0.0
         ENDIF
         
         IF(switch.eq.2)THEN
           do i =1, NumUBin
              IF(UBinL(i).GE.TEnergy .AND. UBinH(i).LT.TEnergy )THEN
                 UBinTime(i) = UBinTime(i) + InterTime*exp( -(UBinL(i)+UBinH(i))/2.0/Temperature )
                 EntroTotTime = EntroTotTime + InterTime*exp(-(UBinL(i)+UBinH(i))/2.0/Temperature )
                 EXIT
              ENDIF
           enddo
         ENDIF
         
         IF(switch.eq.3)THEN
           EntropySys = 0.0e0
           do i =1, NumUBin
              ProbaUBin(i) = UBinTime(i)/EntroTotTime
              IF(ProbaUBin(i).GT.0.0)EntropySys = EntropySys+ProbaUBin(i)*log(ProbaUBin(i))
!              write(*,*)'5940',EntropySys,ProbaUBin(i), log(ProbaUBin(i))
!              read(*,*)
           enddo
           EntropySys = -kB*EntropySys
           DEALLOCATE(UBinL, UBinH, UBinTime,ProbaUBin)
         ENDIF
      
      RETURN
      ENDSUBROUTINE EntropyCalculation

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE MFPMCpick
!    Pick a position for Mass Center based on the Mean Free Path
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE MFPMCpick(iS,OverZ)
      USE Constant_M
      USE FLUID_SOLID_M
      USE PHYSICAL_M
      USE MCSETTING_M
      IMPLICIT NONE 
      
      integer i,iS, FLAG, IteT,IteN, OverZ
      
      real*8 Ave_rho, PathL0, PathL, rdn
      real*8 iSX, iSY, iSZ
      real*8 newiSX, newiSY, newiSZ
      real*8 ZStopRatio,newiSX1, newiSY1, newiSZ1
      real*8 ExpTravelD, Tot_TravelD, TravelD, IniTravelD
      real*8 FromZ,ToZ,LocalDN, Local_rho, IteRatio
      real*8 IteiSX(10),IteiSY(10),IteiSZ(10)
      
      logical IteratePick
      
         TotPickN = TotPickN + 1
         Ave_rho = dble(Npart)/BulkVolume
         PathL0  = 1.0/(sqrt(2.0)*pi*(SigmaFF(1)**2.0)*Ave_rho)
         call random_number(rdn)
         PathL   = PathL0*log(1.0/rdn)
         ExpTravelD = PathL
         Tot_TravelD = 0.0
         FLAG = 0
     !    WRITE(*,*)'6075',Ave_rho,PathL0,PathL
         
!         DO WHILE(FLAG .EQ. 0)
            iSX = 0.0
            iSY = ExpTravelD
            iSZ = 0.0
            CALL ROTATE
         
            newiSX = Rot(1)*iSX + Rot(4)*iSY + Rot(7)*iSZ + MCx(iS)  
            newiSY = Rot(2)*iSX + Rot(5)*iSY + Rot(8)*iSZ + MCy(iS)  
            newiSZ = Rot(3)*iSX + Rot(6)*iSY + Rot(9)*iSZ + MCz(iS)
         
            !-----------------------------------------------------------------------
            ! For the molecules already at the top or bottom,
            ! they only sampling the half sphere that towards the inside of the box
            !----------------------------------------------------------------------
             IF(MeanFPath .AND. (.NOT. pbcz))THEN
  !              if((newiSZ.GT.StopZH).OR.(newiSZ.LT.StopZL))then
  !                  OutZPickN = OutZPickN + 1
  !                 call random_number(rdn)
  !                  newiSX = rdn*BoxLengthX
  !                  call random_number(rdn)
  !                  newiSY = rdn*BoxLengthY
  !                  call random_number(rdn)
  !                  newiSZ = rdn*BoxLengthZ
  !               endif
                
 !               IF((MCz(iS).EQ.StopZH) .AND. (newiSZ.GT.StopZH))THEN  !Only sampling the half sphere inside the box
 !                  newiSZ = 2.0*StopZH - newiSZ
 !               ENDIF
 !               IF((MCz(iS).EQ.StopZL) .AND. (newiSZ.LT.StopZL))THEN
 !                  newiSZ = 2.0*StopZL - newiSZ
 !               ENDIF
 !               newiSZ1 = newiSZ
 !               IF(newiSZ.GT.StopZH) newiSZ1 =  StopZH 
 !               IF(newiSZ.LT.StopZL) newiSZ1 =  StopZL 
 !                  if(newiSZ.NE.MCz(iS).and. (newiSZ1.NE.newiSZ))then
 !                    ZStopRatio = (newiSZ1-MCz(iS))/(newiSZ-MCz(iS))
 !                    newiSX1 = ZStopRatio*(newiSX-MCx(iS))+MCx(iS) 
 !                    newiSY1 = ZStopRatio*(newiSY-MCy(iS))+MCy(iS) 
 !                    newiSX = newiSX1
 !                    newiSY = newiSY1
 !                    newiSZ = newiSZ1
 !                  else
 !                    newiSZ = newiSZ1
 !                  endif  
                  if((newiSZ.GT.StopZH).OR.(newiSZ.LT.StopZL))then
 !                      call random_number(rdn)
 !                      PathL   = PathL0*log(1.0/rdn)
 !                      ExpTravelD = PathL
 !                      CYCLE
                       OverZ = 1
                       return
                  else
                       OverZ = -1
 !                      FLAG = 1
                  endif
             ENDIF
 !            TravelD = sqrt((MCx(iS)-newiSX)**2.0 + (MCy(iS)-newiSY)**2.0 + (MCz(iS)-newiSZ)**2.0)


             IF(pbcZ)THEN
                IF(newiSZ.GT.BoxLengthZ)THEN
                    newiSZ = newiSZ - int(newiSZ/BoxLengthZ)*BoxLengthZ
                ELSEIF(newiSZ.LT.0.0)THEN
                    newiSZ = newiSZ - int(newiSZ/BoxLengthZ)*BoxLengthZ + BoxLengthZ
                ENDIF
             ENDIF
         
             IF(newiSX.GT.BoxLengthX)THEN
                newiSX = newiSX - int(newiSX/BoxLengthX)*BoxLengthX
             ELSEIF(newiSX.LT.0.0)THEN
                newiSX = newiSX - int(newiSX/BoxLengthX)*BoxLengthX + BoxLengthX
             ENDIF
         
             IF(newiSY.GT.BoxLengthY)THEN
                newiSY = newiSY - int(newiSY/BoxLengthY)*BoxLengthY
             ELSEIF(newiSY.LT.0.0)THEN
                newiSY = newiSY - int(newiSY/BoxLengthY)*BoxLengthY + BoxLengthY
             ENDIF
             IniTravelD = sqrt((MCx(iS)-newiSX)**2.0 + (MCy(iS)-newiSY)**2.0 + (MCz(iS)-newiSZ)**2.0)
                        
             IteratePick = .FALSE.
             IteN = 1
             IF(IteratePick)THEN
                DO IteT = 1,IteN
                   if(IteT.eq.1)then
                      IF(newiSZ.GT.MCz(iS))THEN
                         FromZ = MCz(iS)
                         ToZ   = newiSZ
                      ELSE
                         FromZ = newiSZ
                         ToZ   = MCz(iS)
                      ENDIF
                   elseif(IteT.gt.1)then
                      IF(IteiSZ(IteT-1).GT.MCz(iS))THEN
                         FromZ = MCz(iS)
                         ToZ   = IteiSZ(IteT-1)
                      ELSE
                         FromZ = IteiSZ(IteT-1)
                         ToZ   = MCz(iS)
                      ENDIF
                   endif
             
                   LocalDN = 0.0
                   DO i =1, Npart
                      if(MCz(i).GE.FromZ .AND. MCz(i).LE.ToZ)then
                        LocalDN = LocalDN + 1.0
                      endif
                   ENDDO
                   Local_rho = LocalDN/(ToZ-FromZ)/BoxLengthX/BoxLengthY
                   PathL0  = 1.0/(sqrt(2.0)*pi*(SigmaFF(1)**2.0)*Local_rho)
                   call random_number(rdn)
                   PathL   = PathL0*log(1.0/rdn)
                   IteRatio = PathL /IniTravelD 
                   IteiSX(IteT) = MCx(iS)+(newiSX-MCx(iS))*IteRatio
                   IteiSY(IteT) = MCy(iS)+(newiSY-MCy(iS))*IteRatio
                   IteiSZ(IteT) = MCz(iS)+(newiSZ-MCz(iS))*IteRatio
                   newiSZ1 = IteiSZ(IteT)
                   IF(IteiSZ(IteT).GT.StopZH) newiSZ1 =  StopZH 
                   IF(IteiSZ(IteT).LT.StopZL) newiSZ1 =  StopZL 
                   if(newiSZ1.NE.IteiSZ(IteT))then
                     ZStopRatio = (newiSZ1-MCz(iS))/(IteiSZ(IteT)-MCz(iS))
                     newiSX1 = ZStopRatio*(IteiSX(IteT)-MCx(iS))+MCx(iS) 
                     newiSY1 = ZStopRatio*(IteiSY(IteT)-MCy(iS))+MCy(iS) 
                     IteiSX(IteT) = newiSX1
                     IteiSY(IteT) = newiSY1
                     IteiSZ(IteT) = newiSZ1
                   endif  
                 !  WRITE(*,*)'6190',FromZ,ToZ,MCz(iS)
                 !  WRITE(*,*)'6191',Local_rho,PathL0,PathL
                 !  READ(*,*)
                ENDDO
                   IF(IteiSX(IteN).GT.BoxLengthX)THEN
                      IteiSX(IteN) = IteiSX(IteN) - int(IteiSX(IteN)/BoxLengthX)*BoxLengthX
                   ELSEIF(IteiSX(IteN).LT.0.0)THEN
                      IteiSX(IteN) = IteiSX(IteN) - int(IteiSX(IteN)/BoxLengthX)*BoxLengthX + BoxLengthX
                   ENDIF
         
                   IF(IteiSY(IteN).GT.BoxLengthY)THEN
                      IteiSY(IteN) = IteiSY(IteN) - int(IteiSY(IteN)/BoxLengthY)*BoxLengthY
                   ELSEIF(IteiSY(IteN).LT.0.0)THEN
                      IteiSY(IteN) = IteiSY(IteN) - int(IteiSY(IteN)/BoxLengthY)*BoxLengthY + BoxLengthY
                   ENDIF

                   MCx(iS) = IteiSX(IteN)
                   MCy(iS) = IteiSY(IteN)
                   MCz(iS) = IteiSZ(IteN)
              ELSE
                   MCx(iS) = newiSX
                   MCy(iS) = newiSY
                   MCz(iS) = newiSZ
             ENDIF
       !      Tot_TravelD = Tot_TravelD + TravelD
       !      IF(Tot_TravelD.GE.(PathL-1.0E-10))THEN
       !         FLAG = 1
       !         EXIT
       !      ELSEIF(Tot_TravelD.LT.PathL)THEN
       !         ExpTravelD = PathL - Tot_TravelD
       !      ENDIF
!      ENDDO
           
      
      RETURN
      END SUBROUTINE MFPMCpick
      
      

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!    SUBROUTINE StopPosition
!    Calcute the position where the replusion energy equlas to the 
!    kinetic energy of the molecule
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE StopPosition
      USE Constant_M
      USE FLUID_SOLID_M
      USE PHYSICAL_M
      USE MCSETTING_M
      IMPLICIT NONE 
      
      integer i, j, FLAG
      
      real*8 BelowZ, AboveZ, decreaseby, MidZ
      real*8 UC, StopRatio
      
         
           
           IF( (.NOT. Steele).OR.(SigmaSS.EQ.0.0))THEN
              write(62,*)'THERE IS NO SOLID FOR STOP THE PARTICLE'
              StopZL = 0.0E0
              StopZH = BoxLengthZ
           ENDIF
           
           SingleKineticE = 4.0E0/pi*Temperature

           i = Npart+1
           MCx(i) = MCx0
           MCy(i) = MCy0
           MCz(i) = MCz0
           
           do j=1,NLJ
              LJx(i,j) = LJx0(j) + MCx(i)
              LJy(i,j) = LJy0(j) + MCy(i)
              LJz(i,j) = LJz0(j) + MCz(i)
           enddo
           
           BelowZ = 1.0E0
           decreaseby = 0.1E0 
           AboveZ = BelowZ - decreaseby
           FLAG = 0
           DO WHILE(FLAG.eq.0)
              MCz(i) =  AboveZ
              do j=1,NLJ
                 LJz(i,j) = LJz0(j) + MCz(i)
              enddo
              CALL steelepotential(i,UC)
              if(UC.GE.SingleKineticE)then
                 AboveZ = MCz(i)
                 FLAG = 1
              else
                 BelowZ = AboveZ
                 AboveZ = BelowZ - decreaseby
              endif
           ENDDO
           
           FLAG = 0
           DO WHILE(FLAG.eq.0)
              MidZ = (AboveZ+BelowZ)/2.0
              MCz(i) = MidZ
              do j=1,NLJ
                 LJz(i,j) = LJz0(j) + MCz(i)
              enddo
              CALL steelepotential(i,UC)
              StopRatio = (UC-SingleKineticE)/SingleKineticE
              if(abs(StopRatio).LE.0.001E0)then
                 FLAG=1
                 StopZL = MCz(i)
              else
                 if(StopRatio.lt.0.0)BelowZ = MidZ
                 if(StopRatio.gt.0.0)AboveZ = MidZ
              endif      
           ENDDO
           
           StopZH = BoxLengthZ - 1.0E-5    ! minus a very small number to help the program to find the last bin
 !          WRITE(*,*)'6201',StopZL,StopZH
 !          READ(*,*)
 !             StopZL = 0.0E0
 !             StopZH = BoxLengthZ
      
      RETURN
      ENDSUBROUTINE StopPosition
      
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE DISTANCE DISTRIBUTION
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
       SUBROUTINE distancedistribution(switch)
         USE Constant_M
         USE FLUID_SOLID_M
         USE PHYSICAL_M
         USE MCSETTING_M
         USE SUBBOX_M
       implicit none
     
       integer i,j,k
       integer bin, switch,  rhoBinH
       integer AngleBin, Dmaxbin
     
       real*8 distanceZ, DdeltaZ
       real*8 distance, Cy, Cz
       real*8 binL, binH, minBinD, SumTime
       real*8 maxAngle, Angle, deltaZfAngle, dAngle
       real*8, Dimension(:),POINTER::sum_DNbin, ave_DNbin, DNbin
       real*8, Dimension(:),POINTER::dz,dzA,DDensity,DDensityKmolperM3
       real*8, Dimension(:,:),POINTER::sum_DNbinA, ave_DNbinA, DNbinA
       real*8, Dimension(:,:),POINTER::DDensityA,DDensityKmolperM3A

       real*8, Dimension(:,:),POINTER::Sep_Sum_DNbin,Sep_ave_DNbin,Sep_DNbin
       real*8, Dimension(:,:),POINTER::Sep_DDensity,Sep_DDensityKmolperM3
       real*8 sumrho
     
       common/dd/Dmaxbin
       common/dd2/DdeltaZ, SumTime
    
     
       IF(switch.eq.0)THEN
          DdeltaZ         = 0.1e0/scalelength
          Dmaxbin         = int(BoxLengthZ/DdeltaZ)+ 1 
          !IF(NVT)Dmaxbin = 100
          Allocate(sum_DNbin(Dmaxbin), ave_DNbin(Dmaxbin), DNbin(Dmaxbin))
          Allocate(dz(Dmaxbin),dzA(Dmaxbin),DDensity(Dmaxbin),DDensityKmolperM3(Dmaxbin))
          open(unit=31, file='distanceD',status='unknown')
          if(OrienD)then
             deltaAngle = 0.05e0
             if(SymVector)maxAngle = pi/2.0e0
             if(.NOT. SymVector)maxAngle = pi
             maxAngleBin = int(maxAngle/deltaAngle)+1
             Allocate(sum_DNbinA(Dmaxbin,maxAngleBin), ave_DNbinA(Dmaxbin,maxAngleBin), DNbinA(Dmaxbin,maxAngleBin))
             Allocate(DDensityA(Dmaxbin,maxAngleBin),DDensityKmolperM3A(Dmaxbin,maxAngleBin))
          endif
       ENDIF
     
       IF(switch.eq.1)THEN
           sum_DNbin = 0.0e0
           ave_DNbin = 0.0e0
           DNbin     = 0.0e0
           SumTime   = 0.0e0
          write(31,*)'Pressure(pa)   distance(A)   LocalD(Kmol/m3)'
          write(31,*)'==========================================='
          IF(OrienD)THEN
             sum_DNbinA = 0.0e0
             ave_DNbinA = 0.0e0
             DNbinA     = 0.0e0
          ENDIF
       ENDIF
     
       IF(switch.eq.2)THEN
          SumTime = SumTime + InterTime
          DNbin = 0.0e0  
          IF(OrienD)then
           DNbinA     = 0.0e0
          ENDIF 
          DO i = 1, Npart
             DistanceZ  = abs(mcz(i))
             bin = int(DistanceZ/DdeltaZ)+1     
             if(bin.le.Dmaxbin)then
                DNbin(bin) = DNbin(bin) + 1.0e0                
                IF(OrienD)THEN
                   deltaZfAngle = LJz(i,fSite)- LJz(i,tSite)
                   if(deltaZfAngle.ne.0.0)then
                      if(SymVector)then
                        Angle = ACOS(abs(deltaZfAngle)/VectorL)
                      else
                        Angle = ACOS(deltaZfAngle/VectorL)
                      endif
                   else
                      Angle = pi/2.0
                   endif
                   AngleBin = int(Angle/deltaAngle)+1
                   DNbinA(bin, AngleBin) = DNbinA(bin, AngleBin) + 1.0
                ENDIF 
             endif        
          ENDDO     
        
          do i = 1, Dmaxbin
            sum_DNbin(i) = sum_DNbin(i) + DNbin(i)*InterTime
            IF(OrienD)THEN
              do j = 1, maxAngleBin
                sum_DNbinA(i,j) = sum_DNbinA(i,j) + DNbinA(i,j)*InterTime
              enddo
            ENDIF
          enddo
       ENDIF

       IF(switch.eq.3)THEN
          do i = 1, Dmaxbin
            dz(i)  = DdeltaZ*(i-0.5)
            dzA(i) = dz(i)*scalelength
            ave_DNbin(i) = sum_DNbin(i)/SumTime
            DDensity(i)  = ave_DNbin(i)/(BoxLengthX*BoxLengthY*DdeltaZ)
            DDensityKmolperM3(i) = DDensity(i)/6.0230E23/(Scalelength*1.0e-10)**3/1.0e3
            IF(NVT)write(31,'(I4,2E16.6)')Npart, dzA(i),DDensityKmolperM3(i)
            IF(OrienD)THEN
               do j =1, maxAngleBin
                  dAngle = deltaAngle*(j-0.5)
                  ave_DNbinA(i,j) = sum_DNbinA(i,j)/SumTime
                  if(SymVector)DDensityA(i,j) = ave_DNbinA(i,j)/(BoxLengthX*BoxLengthY*DdeltaZ*SIN(dAngle)*deltaAngle)
                  if(.NOT. SymVector)DDensityA(i,j) = 2.0*ave_DNbinA(i,j)/(BoxLengthX*BoxLengthY*DdeltaZ*SIN(dAngle)*deltaAngle)
                  DDensityKmolperM3A(i,j) = DDensityA(i,j)/6.0230E23/(Scalelength*1.0e-10)**3/1.0e3
                  if(NVT)write(36,'(I8,4E16.6)')Npart, dzA(i),dAngle,DDensityA(i,j),DDensityKmolperM3A(i,j)
               enddo   
            ENDIF
          enddo 
       ENDIF
     
       IF(switch.eq.4)THEN
         DEALLOCATE(sum_DNbin, ave_DNbin, DNbin)
         DEALLOCATE(dz,dzA,DDensity,DDensityKmolperM3)
         if(OrienD)then
           DEALLOCATE(sum_DNbinA, ave_DNbinA, DNbinA)
           DEALLOCATE(DDensityA,DDensityKmolperM3A)
         endif
       ENDIF
     
     RETURN
     ENDSUBROUTINE distancedistribution    
      

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE PoreMFPMCpick
!    Pick a position for Mass Center based on the Mean Free Path for "PORE"
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE PoreMFPMCpick(iS,OverZ)
      USE Constant_M
      USE FLUID_SOLID_M
      USE PHYSICAL_M
      USE MCSETTING_M
      IMPLICIT NONE 
      
      integer i,iS, FLAG, IteT,IteN, OverZ
      
      real*8 Ave_rho, PathL0, PathL, rdn
      real*8 iSX, iSY, iSZ
      real*8 newiSX, newiSY, newiSZ
      real*8 ZStopRatio,newiSX1, newiSY1, newiSZ1
      real*8 ExpTravelD, Tot_TravelD, TravelD
      real*8 FromZ,ToZ,LocalDN, Local_rho, IteRatio
      
      
         TotPickN = TotPickN + 1
         Ave_rho = dble(Npart)/BulkVolume
         PathL0  = 1.0/(sqrt(2.0)*pi*(SigmaFF(1)**2.0)*Ave_rho)
         call random_number(rdn)
         PathL   = PathL0*log(1.0/rdn)
         ExpTravelD = PathL
         FLAG = 0
         
            iSX = 0.0
            iSY = ExpTravelD
            iSZ = 0.0
            CALL ROTATE
         
            newiSX = Rot(1)*iSX + Rot(4)*iSY + Rot(7)*iSZ + MCx(iS)  
            newiSY = Rot(2)*iSX + Rot(5)*iSY + Rot(8)*iSZ + MCy(iS)  
            newiSZ = Rot(3)*iSX + Rot(6)*iSY + Rot(9)*iSZ + MCz(iS)
         
            !-----------------------------------------------------------------------
            ! For the molecules already at the top or bottom,
            ! they only sampling the half sphere that towards the inside of the box
            !----------------------------------------------------------------------
             IF(MeanFPath .AND. (.NOT. pbcz))THEN
                 if( ((newiSY.GE.AccVYL).AND.(newiSY.LE.AccVYH))  & 
             &      .AND. ((newiSZ.GT.AccVZH).OR.(newiSZ.LT.AccVZL)) )then
                       OverZ = 1
                       return
                  elseif(newiSY.LT.0.0E0 .OR. newiSY.GT.BoxLengthY)then
                       OverZ = 1
                       return
                  else
                       OverZ = -1
                  endif
             ENDIF

             
            IF(newiSZ.GT.BoxLengthZ)THEN
                newiSZ = newiSZ - int(newiSZ/BoxLengthZ)*BoxLengthZ
            ELSEIF(newiSZ.LT.0.0)THEN
                newiSZ = newiSZ - int(newiSZ/BoxLengthZ)*BoxLengthZ + BoxLengthZ
            ENDIF
             
         
             IF(newiSX.GT.BoxLengthX)THEN
                newiSX = newiSX - int(newiSX/BoxLengthX)*BoxLengthX
             ELSEIF(newiSX.LT.0.0)THEN
                newiSX = newiSX - int(newiSX/BoxLengthX)*BoxLengthX + BoxLengthX
             ENDIF
         
             IF(pbcy)THEN
                 IF(newiSY.GT.BoxLengthY)THEN
                    newiSY = newiSY - int(newiSY/BoxLengthY)*BoxLengthY
                 ELSEIF(newiSY.LT.0.0)THEN
                    newiSY = newiSY - int(newiSY/BoxLengthY)*BoxLengthY + BoxLengthY
                 ENDIF
             ENDIF
                                     

             MCx(iS) = newiSX
             MCy(iS) = newiSY
             MCz(iS) = newiSZ
                        
      
      RETURN
      END SUBROUTINE PoreMFPMCpick
      
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE DDistribution_2D
!      FOR CALCULATING THE LOCAL DENSITY DISTRIBUTION IN 2D       
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     subroutine DDistribution_2D(switch)
     USE Constant_M
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE F2DLOCALD_M
     USE SUBBOX_M
     implicit none
     
     integer i,j,k
     integer switch
     integer i1, i2, i3, iCount
     integer iBin, inX, inY, inZ, iSec 
     integer MaxBinfAve
     
     real*8 Mindelta, NDD 
     
     common/d2d/NDD    
     
          
     IF(switch.eq.0)THEN
          unitNo = 200
          open(unit=25, file='Ini2DLocalD.txt', status='unknown')
          rewind(25)
          read(25,*)Sec2D
          read(25,*)(SecYL2D(i),i=1,Sec2D)
          read(25,*)(SecYH2D(i),i=1,Sec2D)
          read(25,*)(SecZL2D(i),i=1,Sec2D)
          read(25,*)(SecZH2D(i),i=1,Sec2D)
          read(25,*)(deltaY2D(i),i=1,Sec2D) ! Input in reduced unit
          read(25,*)(deltaZ2D(i),i=1,Sec2D) ! Input in reduced unit
          read(25,*) AveRadius              ! Input in reduced unit
          DO i = 1, Sec2D
             SecYL2D(i) = SecYL2D(i)/ScaleLength
             SecYH2D(i) = SecYH2D(i)/ScaleLength
             SecZL2D(i) = SecZL2D(i)/ScaleLength
             SecZH2D(i) = SecZH2D(i)/ScaleLength
          ENDDO
          totBin2D = 0
          IF(deltaY2D(1).LT.deltaZ2D(1))THEN
              Mindelta = deltaY2D(1)
          ELSE
              Mindelta = deltaZ2D(1)
          ENDIF
          DO i = 1, Sec2D
             IF(mod((SecYH2D(i)-SecYL2D(i)),deltaY2D(1)) .eq. 0.0)THEN
                 NinYSec2D(i) = int((SecYH2D(i)-SecYL2D(i))/deltaY2D(1))
             ELSE
                 NinYSec2D(i) = int((SecYH2D(i)-SecYL2D(i))/deltaY2D(1)) + 1
                 deltaY2D(i)  = (SecYH2D(i)-SecYL2D(i))/NinYSec2D(i)
                 IF(deltaY2D(i).LT.Mindelta)Mindelta = deltaY2D(i)
             ENDIF
             IF(mod((SecZH2D(i)-SecZL2D(i)),deltaZ2D(1)) .eq. 0.0)THEN
                 NinZSec2D(i) = int((SecZH2D(i)-SecZL2D(i))/deltaZ2D(1))
             ELSE 
                 NinZSec2D(i) = int((SecZH2D(i)-SecZL2D(i))/deltaZ2D(1)) + 1
                 deltaZ2D(i)  = (SecZH2D(i)-SecZL2D(i))/NinZSec2D(i)
                 IF(deltaZ2D(i).LT.Mindelta)Mindelta = deltaZ2D(i)
             ENDIF
             Bin2DSec(i)  =  NinYSec2D(i) * NinZSec2D(i)            
             BinVolSec(i) = BoxLengthX*deltaY2D(i)*deltaZ2D(i)
             totBin2D = totBin2D + Bin2DSec(i)
           ENDDO
           
           MaxBinfAve = (int(2.0*AveRadius/Mindelta)+2)**2.0
           
           ALLOCATE(Nin2DBin(totBin2D), Sum_Nin2DBin(totBin2D), ave_Nin2DBin(totBin2D))
           ALLOCATE(Bin2DY(totBin2D),Bin2DZ(totBin2D))
           ALLOCATE(DDensity2D(totBin2D),Smooth_DDensity2D(totBin2D),DDensity2DKmolperM3(totBin2D))
           ALLOCATE(SumD_Bin2D(totBin2D), Bin2DVol(totBin2D), SumV_Bin2D(totBin2D))
           ALLOCATE(NBinfAve(totBin2D), SeqBinfAve(totBin2D,MaxBinfAve))
           
           iCount = 0
           DO i = 1, Sec2D
             if(i.gt.1)iCount = iCount + Bin2DSec(i-1)
                do i1 = 1, NinYSec2D(i)
                   do i2 = 1, NinZSec2D(i)
                      iBin = iCount + (i1-1)*NinZSec2D(i) + i2
                      Bin2DY(iBin) = (i1-1)*deltaY2D(i) + deltaY2D(i)/2.0 + SecYL2D(i)
                      Bin2DZ(iBin) = (i2-1)*deltaZ2D(i) + deltaZ2D(i)/2.0 + SecZL2D(i)
                      Bin2DVol(iBin) = BinVolSec(i)
                   enddo
                enddo
           ENDDO  
          
          CALL FindBinfAve
          Nin2DBin = 0
     ENDIF
     
     IF(switch.eq.1)THEN
         sum_Nin2DBin = 0.0e0
         ave_Nin2DBin = 0.0e0
         NDD          = 0.0e0
     ENDIF
     
     IF(switch.eq.2)THEN
        NDD = NDD + InterTime
        Nin2DBin  = 0 
        DO i = 1, Npart
           iSec = 0
           do i1=1,Sec2D
              if( MCy(i).GE.SecYL2D(i1) .AND. MCy(i).LT.SecYH2D(i1) )then
                 iSec = i1
                 EXIT
              endif
           enddo
           IF(iSec.gt.0)THEN
              inY = int( (MCy(i)-SecYL2D(iSec))/deltaY2D(iSec) ) + 1
              inZ = int( (MCz(i)-SecZL2D(iSec))/deltaZ2D(iSec) ) + 1
              iBin = 0
              do i1 = 1, iSec - 1
                 iBin = iBin + Bin2DSec(i1)
              enddo
              iBin = iBin + (inY-1)*NinZSec2D(iSec) + inZ
              Nin2DBin(iBin) = Nin2DBin(iBin) + 1
           ENDIF
        ENDDO     
        
        DO i = 1, totBin2D
           sum_Nin2DBin(i) = sum_Nin2DBin(i) + Nin2DBin(i)*InterTime
        ENDDO
     ENDIF

     IF(switch.eq.3)THEN    
          do i = 1, totBin2D
             ave_Nin2DBin(i) = sum_Nin2DBin(i)/NDD
         enddo 
         CALL LocalD2DSmooth 
     ENDIF
     
     IF(switch.eq.4)THEN
        DEALLOCATE(Nin2DBin, Sum_Nin2DBin, ave_Nin2DBin)
        DEALLOCATE(Bin2DY,Bin2DZ)
        DEALLOCATE(DDensity2D,Smooth_DDensity2D,DDensity2DKmolperM3)
        DEALLOCATE(SumD_Bin2D, Bin2DVol, SumV_Bin2D)
        DEALLOCATE(NBinfAve, SeqBinfAve)
     ENDIF
     
   RETURN
   ENDSUBROUTINE DDistribution_2D
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE FindBinfAve
!      Find the bins in the range for smoothing the 2D Local density distribution 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     subroutine FindBinfAve
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE F2DLOCALD_M
     implicit none
     
     integer i, j
     
     real*8 distance2Bin
     
         
     NBinfAve   = 0
     SeqBinfAve = 0
     SumV_Bin2D = 0.0E0
       
    
     DO i = 1, totBin2D-1
        DO j = i+1, totBin2D
  !         distance2Bin = ((Bin2DY(i)-Bin2DY(j))**2.0e0 + (Bin2DZ(i)-Bin2DZ(j))**2.0e0)**0.5e0
           IF( abs(Bin2DY(i)-Bin2DY(j)).LE. AveRadius .AND.   & 
       &       abs(Bin2DZ(i)-Bin2DZ(j)).LE. AveRadius   )  THEN
  !          IF(distance2Bin .LE.AveRadius)THEN
              NBinfAve(i) = NBinfAve(i) + 1
              SeqBinfAve(i, NBinfAve(i)) = j
              SumV_Bin2D(i) = SumV_Bin2D(i) + Bin2DVol(j)
              NBinfAve(j) = NBinfAve(j) + 1
              SeqBinfAve(j, NBinfAve(j)) = i
              SumV_Bin2D(j) = SumV_Bin2D(j) + Bin2DVol(i)
           ENDIF    
       ENDDO
     ENDDO

     RETURN
     ENDSUBROUTINE FindBinfAve   
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE LocalD2DSmooth
!      Smooth the 2D Local density distribution by averaging one bin with its      
!      neighbour bins
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     subroutine LocalD2DSmooth
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE F2DLOCALD_M
     implicit none
     
     integer i, j, fixUnit
     
     character*200 filename, filename1, filename2
     
         
     SumD_Bin2D  = 0.0e0
           
     DO i = 1, totBin2D
        DO j = 1, NBinfAve(i)
            SumD_Bin2D(i)  = SumD_Bin2D(i) + ave_Nin2DBin(SeqBinfAve(i,j))
       ENDDO
     ENDDO
     
     unitNo = unitNo + 1
     fixUnit = unitNo
     write(filename1,'(I8)') Npart
     write(filename2,'(I5)') iden
     filename = trim(filename1)//trim(filename2)//'Pa2D.txt'
     open(unit=fixUnit, file=filename,status='unknown')
     rewind(fixUnit)
     write(fixUnit,'(I4,I15)') Npart,totBin2D
     
     Do i = 1, totBin2D
        Smooth_DDensity2D(i) = ( SumD_Bin2D(i) + ave_Nin2DBin(i) )/(SumV_Bin2D(i) + Bin2DVol(i))
        DDensity2DKmolperM3(i) = Smooth_DDensity2D(i)/6.0230E23/(Scalelength*1.0e-10)**3/1.0e3
        write(fixUnit,'(I6,5E18.8)')Npart, Bin2DY(i),Bin2DZ(i),DDensity2DKmolperM3(i), ave_Nin2DBin(i), Bin2DVol(i)
     ENDDO

     RETURN
     ENDSUBROUTINE LocalD2DSmooth    