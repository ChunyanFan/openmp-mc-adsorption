!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE BinJudge
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     subroutine BinJudge(x1,y1,z1,bin)
     USE MCSETTING_M
     USE FLUID_SOLID_M
     USE FLUCTUATION_M
     implicit none
     
     integer i, j, k, bin, switch, NDD
     
     real*8 distanceZ 
     real*8 distance, Cy, Cz
     real*8 x1,y1,z1
      
         
!      deltaZ  = 0.5e0/scalelength
      DistanceZ = abs(Z1)
           
           if(Bojan)then
              do j = 1, NS
                if( y1.ge. (StripC(j)-StripW(j)/2.0e0).and. &
                  & y1.le. (StripC(j)+StripW(j)/2.0e0))then
                       distance = abs(z1- StripZ0(j))                                            
                       if(distance .lt. DistanceZ)DistanceZ = distance
                endif
                 if(y1.lt. (StripC(j)-StripW(j)/2.0e0))then
                       Cy = abs(y1-(StripC(j)-StripW(j)/2.0e0))
                 else
                         Cy = abs(y1-(StripC(j)+StripW(j)/2.0e0))
                      do k = 1,StriplayerN(j)
                         Cz = abs(z1-(StripZ0(j)-(k-1)*Stripgap(j)))
                         distance = sqrt(Cy**2.0e0 + Cz**2.0e0)
                         if(distance .lt. DistanceZ)DistanceZ = distance
                      enddo
                endif         
              enddo
           endif
           
           if(NC3 .gt. 0)then
              do j = 1, NC3
                distance = sqrt((x1-xc(j))**2.0e0 + (y1-yc(j))**2.0e0 + (z1-zc(j))**2.0e0)
                if(distance .lt. DistanceZ)DistanceZ = distance
              enddo
           endif
           
          bin = int(DistanceZ/deltaZ)+1     
       
   return
   end           
            
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Particleposition
!      Check if the particle is in the pore or in the gas phase
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     subroutine Particleposition(i,flag)
     USE MCSETTING_M
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE SUBBOX_M
     implicit none
     
     integer i, flag

        if(MCy(i).ge. AccVYL .and. MCy(i).le. AccVYH)then
           if(ExtraBin)then
             if(MCz(i).le. BoxLengthZ)then
                if(flag.eq.1)Ninpore = Ninpore + 1
                if(flag.eq.-1)Ninpore = Ninpore - 1
             endif
           else
                if(flag.eq.1)Ninpore = Ninpore + 1
                if(flag.eq.-1)Ninpore = Ninpore - 1
           endif
        endif
       
       IF(DecayESF)THEN
          if(MCy(i).ge. CollectLY .and. MCy(i).le. CollectHY)then
            if(flag.eq.1)CollectN = CollectN + 1
            if(flag.eq.-1)CollectN = CollectN - 1
          endif
       ENDIF
       
     return
     end subroutine particleposition         
    
     
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE CLMove
!    CHANGE THE LOCATION OF CHARGES ON ADSORBATE MOLECULES
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       SUBROUTINE CLMove(TEnergy,Energy,CLEnergy,EnergyC)
        USE FLUID_SOLID_M
        USE PHYSICAL_M
        USE MCSETTING_M
        USE FLUCTUATION_M
        USE POREFIGURE_M
        USE  FLEXICLINFO_M
       implicit none
       
       integer N, i,j, new
       
       real*8 TEnergy, rdn, DeltaEnergy
       real*8 RdmCLStep, OldD, OldD2
       real*8 OldIntraE, NewIntraE, CoulombOld, CoulombNew
       real*8 TEnergyOld, TEnergyNew
       real*8 Energy, CLEnergy, EnergyC
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       
       if(Npart.eq.0)return
       N    = Npart
       new = N+1
!      ---------------------------
!      RANDOMLY SELECT A PARTICLE 
!      ---------------------------
       call random_number(rdn)
       i  = int(rdn*N) + 1
!      ------------------------------------------------------------------------------
!      THIS IS TO ENSURE THAT i SHOULD NOT EXCEED N, AS THIS CAN OCCURS WHEN rand()=1
!      ------------------------------------------------------------------------------
       if(i.gt.N) i = N 
       ! CALCULATE THE IntraE of i first
       OldIntraE = IntraE(i)
!       CALL SingleCoulomb(i, CoulombOld)
       call energysingleparticle(1,i, EnergyOld, CLEnergyOld, EnergyCOld)
       TEnergyOld = EnergyOld+ CLEnergyOld+ EnergyCOld
!      -------------------------------------------------
!      INSTORE THE OLD POSITION OF picked CHARGES FIRST
!      -------------------------------------------------
       call random_number(rdn)
       j = int(rdn*NpickedCL) + 1
       if (j.gt.NpickedCL) j = NpickedCL

!       do j =1, NpickedCL   
          CLx(new,indexPickedCL(j))  = CLx(i,indexPickedCL(j))
          CLy(new,indexPickedCL(j))  = CLy(i,indexPickedCL(j))
          CLz(new,indexPickedCL(j))  = CLz(i,indexPickedCL(j))

!       enddo
!      ------------------------------------------------
!      SET THE STEP FOR MOVEING THE LOCATION OF CHARGES
!      ------------------------------------------------
       call random_number(rdn)
       RdmCLStep = (rdn-0.5)*CLStepRange
       Lequi(i) = Lequi(i) + RdmCLStep
!      ----------------------------------------------------
!      NOW MOVE THE CHARGES 
!      ------------------------------------------------------------
!      1ST: ASSUMEED MOVE THE MASS CENTER OF PARTICLE i TO ORIGINAL
!      2ND: MOVE THE CHARGES
!      3RD: MOVED THE MASS CENTER BACK TO ITS REAL POSITION
!      ------------------------------------------------------------
!       do j = 1, NpickedCL
          ECLx(i,indexPickedCL(j)) = ECLx(i,indexPickedCL(j))- mcx(i) ! NOTE: the initial mass center of particle (MCX0,MCY0,MCZ0) has to be origianl point
          ECLy(i,indexPickedCL(j)) = ECLy(i,indexPickedCL(j))- mcy(i)
          ECLz(i,indexPickedCL(j)) = ECLz(i,indexPickedCL(j))- mcz(i)
          OldD2 = ECLx(i,indexPickedCL(j))**2 + ECLy(i,indexPickedCL(j))**2 + ECLz(i,indexPickedCL(j))**2
          OldD = sqrt(OldD2)
          CLx(i,indexPickedCL(j)) = ECLx(i,indexPickedCL(j))*(1.0 + RdmCLStep/OldD) + mcx(i)
          CLy(i,indexPickedCL(j)) = ECLy(i,indexPickedCL(j))*(1.0 + RdmCLStep/OldD) + mcy(i)
          CLz(i,indexPickedCL(j)) = ECLz(i,indexPickedCL(j))*(1.0 + RdmCLStep/OldD) + mcz(i)
          LJx(i,indexPickedLJ(j)) = CLx(i,indexPickedCL(j))
          LJy(i,indexPickedLJ(j)) = CLy(i,indexPickedCL(j))
          LJz(i,indexPickedLJ(j)) = CLz(i,indexPickedCL(j))
          ECLx(i,indexPickedCL(j)) = ECLx(i,indexPickedCL(j))+ mcx(i) ! MOVE BACK
          ECLy(i,indexPickedCL(j)) = ECLy(i,indexPickedCL(j))+ mcy(i)
          ECLz(i,indexPickedCL(j)) = ECLz(i,indexPickedCL(j))+ mcz(i)
!       enddo
       NewIntraE = MaxPenalty*(1-Lequi(i)/Lequi0)**2.0e0
!       CALL SingleCoulomb(i, CoulombNew)
       call energysingleparticle(1,i, EnergyNew, CLEnergyNew, EnergyCNew)
       TEnergyNew = EnergyNew+ CLEnergyNew+ EnergyCNew
       DeltaEnergy = NewIntraE - OldIntraE + TEnergyNew-TEnergyOld
       call random_number(rdn)
       IF(rdn .LT. exp(-DeltaEnergy))THEN  !Accept the move of charges
          TEnergy = TEnergy + DeltaEnergy
          Energy  = Energy - EnergyOld + EnergyNew
          CLEnergy  = CLEnergy - CLEnergyOld + CLEnergyNew
          EnergyC  = EnergyC - EnergyCOld + EnergyCNew
          IntraE(i) = NewIntraE
       ELSE
 !         do j = 1, NpickedCL
             CLx(i,indexPickedCL(j)) = CLx(new,indexPickedCL(j))
             CLy(i,indexPickedCL(j)) = CLy(new,indexPickedCL(j))
             CLz(i,indexPickedCL(j)) = CLz(new,indexPickedCL(j))
             LJx(i,indexPickedLJ(j)) = CLx(i,indexPickedCL(j))
             LJy(i,indexPickedLJ(j)) = CLy(i,indexPickedCL(j))
             LJz(i,indexPickedLJ(j)) = CLz(i,indexPickedCL(j))
             Lequi(i) = Lequi(i) - RdmCLStep
 !         enddo
      ENDIF 
       
      return
      END SUBROUTINE CLMove
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE SingleCoulomb
!    Calculte the Coulomb force of single particle with all the others
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE SingleCoulomb(i, SCoulomb)
        USE FLUID_SOLID_M
        USE PHYSICAL_M
        USE MCSETTING_M
        USE FLUCTUATION_M
        USE POREFIGURE_M
        USE  FLEXICLINFO_M
      implicit none
      
        integer i, j, k, m
        
        real*8 SCoulomb
        real*8 xcl,ycl,zcl
        real*8 cld, CLU, CLV
        real*8 xmclj, ymclj, zmclj
        
        SCoulomb = 0.0e0
        do j = 1,Npart
           if(j .NE. i)then
                xmclj = abs(mcx(j)-mcx(i))
                ymclj = abs(mcy(j)-mcy(i))
                zmclj = abs(mcz(j)-mcz(i))
             do k = 1, NCL
                 do m = 1, NCL
                    xcl =  abs(CLx(j,m) - CLx(i,k))
                    ycl =  abs(CLy(j,m) - CLy(i,k))
                    zcl =  abs(CLz(j,m) - CLz(i,k))
                    if(pbcx)then                
                       if(xmclj > CarbonLengthx1/2.0D0) xcl = xcl - CarbonLengthx1
                    endif
                    if(pbcy)then
                       if(ymclj > CarbonLengthy1/2.0D0) ycl = ycl - CarbonLengthy1
                    endif
                    if(pbcz)then
                       if(zmclj > BoxLengthZ/2.0D0) zcl = zcl - BoxLengthZ
                    endif
                    cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                    call coulombforce(cld, charges(k),charges(m),CLU,CLV)
                    SCoulomb = SCoulomb + CLU
                 enddo
              enddo   
            endif
         enddo
         
      return
      END SUBROUTINE SingleCoulomb
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE InserNVT
!    INSERT THE NEW ADDED PARTICLES IN THE SIMULATION BOX
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE InsertNVT   
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer i,j, flag,k
      real*8 rdn, distance2,distance
      
        i=1       
        do while(i .le. addN)
           call random_number(rdn)
           MCx(Npart+i) = rdn*Carbonlengthx1
           call random_number(rdn)
           MCy(Npart+i) = rdn*Carbonlengthy1
           call random_number(rdn)
           MCz(Npart+i) = rdn*Boxlengthz
          
           flag=0 
           
           do k=1,NSteele
              if( abs(MCz(Npart+i)-PSteele(k)) .LT. 0.75e0 )then
                flag=1
                EXIT
              endif
           enddo
           if(flag.eq.1)CYCLE

            
           do j=1,Npart+i-1
              distance2 = (MCx(Npart+i)-MCx(j))**2.0 + (MCy(Npart+i)-MCy(j))**2.0 + (MCz(Npart+i)-MCz(j))**2.0
              distance = sqrt(distance2)
              if(distance.LT.0.75E0)then
                 flag=1
                 EXIT
              endif
          enddo
          if(flag.eq.1)CYCLE
                     
           if(SDisplace)then
              flag=0
             do j= 1, SubBox
                if(MCz(Npart+i).GE. LayerZL(j) .AND. MCz(Npart+i).LT. LayerZH(j) )then
                   NinSBox(j)= NinSBox(j)+1
                   SeqNinSBox(j,int(NinSBox(j)))=Npart+i
                   indNinSBox(Npart+i)=j
                   flag=1
                endif
                if(flag.eq.1)EXIT
             enddo
           endif

            do j = 1, NLJ
             LJx(Npart+i,j) = LJx0(j) + MCx(Npart+i) - MCx0
             LJy(Npart+i,j) = LJy0(j) + MCy(Npart+i) - MCy0
             LJz(Npart+i,j) = LJz0(j) + MCz(Npart+i) - MCz0
            enddo
         
           do j =1, NCL
             CLx(Npart+i,j) = CLx0(j) + MCx(Npart+i) - MCx0
             CLy(Npart+i,j) = CLy0(j) + MCy(Npart+i) - MCy0
             CLz(Npart+i,j) = CLz0(j) + MCz(Npart+i) - MCz0
           enddo
           
           i=i+1

      enddo
      
     return
     END SUBROUTINE InsertNVT

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE DelNVT
!    Delete PARTICLES FROM SIMULATION BOX - DESORPTION FOR CMC
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE DelNVT   
        USE FLUID_SOLID_M
        USE MCSETTING_M
        USE PHYSICAL_M
        USE SUBBOX_M
      implicit none
      
      integer idel,i,j,k, DelN,ID
      integer aimSub
      
      real*8 rdn
      
      DelN = addN
      IF(NinSBox(SubBox).LT.addN)DelN=NinSBox(SubBox)
      IF(DelN.LT.1)THEN
        write(15,*)'Oops, the number of partics in gas bin are 0!'
        STOP
      ENDIF      
      idel=1       
      do while(idel .le. DelN)
          call random_number(rdn)
          ID  = int(rdn*NinSBox(SubBox)) + 1
          if(ID.gt.NinSBox(SubBox)) ID = NinSBox(SubBox)
          i=SeqNinSBox(SubBox,ID)        
          SeqNinSBox(SubBox,ID)= SeqNinSBox(SubBox,int(NinSBox(SubBox)))
          NinSBox(SubBox) = NinSBox(SubBox)-1  
          aimSub = IndNinSBox(Npart)
          do j=1, NinSBox(aimSub)
             if(SeqNinSBox(aimSub,j).eq.Npart)then
                SeqNinSBox(aimSub,j)=i
                EXIT
             endif
           enddo
           IndNinSBox(i)=aimSub
           SFenergy(i) = SFenergy(Npart)
           call particleposition(i,-1)
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
     END SUBROUTINE DelNVT


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE NBulkCheck
!    COUNT THE NUMBER OF PARTICLES IN THE BULK PHASE FOR NVT ENSEMBLE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE NBulkCheck(Ncount)   
      USE MCSETTING_M
      USE FLUID_SOLID_M
      USE PHYSICAL_M
      implicit none
      
      integer i, Ncount
      
      Ncount = 0
      
      do i=1, Npart
         if((MCz(i).LE.AccVZH))Ncount = Ncount+1
      enddo  
        
      return
      END SUBROUTINE NBulkCheck
      
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE SecMCMOVE
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine SecMove(TEnergy,Energy,CLEnergy,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer N, i,J,new,ID, flag, repeatN
       integer iSub, fSub,tSub
       
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       real*8 TEnergy, TEnergyOld, TEnergyNew
       real*8 Energy,EnergyC, CLEnergy
       real*8 deltaEnergy, rdn, rdn1
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 ForLog
       
       logical inMove, Movein, Moveout
        
       N = Npart
       if(N.EQ.0)return      
!      -------------------------
!      PICK ONE SUBBOX RANDOMLY
!      -------------------------
       call random_number(rdn)
       iSub = int(rdn*SubBox) + 1
       if(iSub.gt.SubBox) iSub = SubBox      
!      ---------------------------------------------------------------------------
!      Equal possibility for move in, move out or displacement in the same SubBox
!      ---------------------------------------------------------------------------
       call random_number(rdn1)
       IF(rdn1.le.0.75e0)THEN         !Pick one particle and move out or displace in the same SubBox
          if(NinSBox(iSub).eq.0)return 
          call random_number(rdn)
          ID  = int(rdn*NinSBox(iSub)) + 1
          if(ID.gt.NinSBox(iSub)) ID = NinSBox(iSub)
          i=SeqNinSBox(iSub,ID)  
          call random_number(rdn)
          if(rdn.le.0.5e0)then  !displace in the same SubBox
            inMove  = .TRUE.
            Moveout = .FALSE.  
            Movein  = .FALSE.
            AtpinMove(iSub) = AtpinMove(iSub) + 1
          else     ! Move out from this SubBox
            inMove  = .FALSE.
            Moveout = .TRUE. 
            Movein  = .FALSE. 
          endif
       ELSE                    !Pick one particle from outside to this bin
            inMove  = .FALSE.
            Moveout = .FALSE. 
            Movein  = .TRUE. 
            flag=0
            do while(flag.eq.0)
              call random_number(rdn)
              fSub = int(rdn*SubBox) + 1
              if(fSub.gt.SubBox) fSub = SubBox      
              if(fSub .NE. iSub)then
                flag=1
                EXIT
              endif
           enddo
           if(NinSBox(fSub).eq.0)return
           call random_number(rdn)
           ID  = int(rdn*NinSBox(fSub)) + 1
           if(ID.gt.NinSBox(fSub)) ID = NinSBox(fSub)
           i=SeqNinSBox(fSub,ID)  
       ENDIF
!      --------------------------------------
!      OLD ENERGY OF PARTICLE i IN ADSORPTION
!      --------------------------------------
          call energySingleParticle(1,i, EnergyOld, CLEnergyOld, EnergyCOld )
          TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld
!         -------------------------------------------------
!         Store the information for i before any movement
!         ------------------------------------------------- 
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
           enddo

         IF(inMove)THEN
             call random_number(rdn)
             mcx(i) = mcx(i) + (rdn - 0.5E0)*SStepX(iSub)
             if(.NOT. pbcx)then
               if((mcx(i) .lt. 0.E0) .or. (mcx(i) .gt. CarbonLengthx1))then
                  mcx(i)=mcx(new)
                  return
               endif 
             endif
       
             call random_number(rdn)
             mcy(i) = mcy(i) + (rdn - 0.5E0)*SStepY(iSub)  
             if(.NOT. pbcy)then    
               if((mcy(i) .lt. 0.E0) .or. (mcy(i) .gt. CarbonLengthy1))then
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  return
               endif
             endif
        
             call random_number(rdn)
             mcz(i) = mcz(i) + (rdn - 0.5E0)*SStepZ(iSub)
             if(mcz(i).lt.LayerZL(iSub) .OR. mcz(i).ge.LayerZH(iSub))then
                mcx(i) = mcx(new)
                mcy(i) = mcy(new)
                mcz(i) = mcz(new)
                return
             endif
             if(.NOT. pbcz)then    
               if((mcz(i) .lt. 0.E0) .or. (mcz(i) .gt. BoxLengthZ))then
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  mcz(i) = mcz(new)
                  return
               endif
              endif
!                if (mcz(i).gt.LayerZH(iSub)) then                     
!                   mcz(i) = mcz(i) - int((mcz(i)-LayerZL(iSub))/(LayerZH(iSub)-LayerZH(iSub)))*(LayerZH(iSub)-LayerZH(iSub))
!                elseif (mcz(i).lt.LayerZL(iSub)) then
!                   mcz(i) = mcz(i) - int((mcz(i)-LayerZL(iSub))/(LayerZH(iSub)-LayerZH(iSub)))*(LayerZH(iSub)-LayerZH(iSub)) + (LayerZH(iSub)-LayerZH(iSub))
!                endif
          ELSEIF(Moveout)THEN
                fSub = iSub
                call random_number(rdn)
                mcx(i) = rdn*CarbonLengthx1
                call random_number(rdn)
                mcy(i) = rdn*CarbonLengthy1
                flag=0
                do while(flag.eq.0)
                   call random_number(rdn)
                   mcz(i)=rdn*BoxLengthZ
                   if(mcz(i).lt.LayerZL(iSub) .OR. mcz(i).ge.LayerZH(iSub))then
                     flag=1
                     EXIT
                   endif
                enddo
                do j=1,SubBox
                   if(mcz(i).ge. LayerZL(j) .and. mcz(i).lt. LayerZH(j))then
                      tSub = j
                       EXIT
                   endif
                enddo
          ELSEIF(Movein)then
              tSub = iSub
              call random_number(rdn)
              mcx(i) = rdn*CarbonLengthx1
              call random_number(rdn)
              mcy(i) = rdn*CarbonLengthy1
              call random_number(rdn)
              mcz(i) = LayerZL(iSub) + rdn*(LayerZH(iSub)-LayerZL(iSub)) 
          ENDIF
!           if(Bojan)then
!                    do j = 1, NS
!                       if( StripZ0(j) .gt. 0.0e0 )then
!                         if(   mcz(i) .lt. StripZ0(j) .and. &
!                           &   mcy(i) .gt. (StripC(j)- StripW(j)/2.0e0)  .and. &
!                           &   mcy(i) .lt. (StripC(j)+ StripW(j)/2.0e0) ) then
!                           mcx(i) = mcx(new)
!                           mcy(i) = mcy(new)
!                           mcz(i) = mcz(new)
!                          return  
!                         endif
!                       endif
!                    enddo
!           endif 
            
       
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
          enddo     

!          -------------------------
!          NEW ENERGY OF PARTICLE I
!          -------------------------
           call  energySingleParticle(1,i, EnergyNew, CLEnergyNew, EnergyCNew)
           TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
       
!          ------------------------------------------------
!          CHECK TO SEE WHETHER THE MOVE SHOULD BE ACCEPTED
!          ------------------------------------------------
           deltaEnergy =  (TEnergyNew - TEnergyOld)/Temperature
           if(Movein .OR. Moveout)then
               ForLog = SubVol(fSub)*(NinSBox(tSub)+1)/SubVol(tSub)/NinSBox(fSub)
!               ForLog = (BulkVolume-SubVol(tSub))*(NinSBox(tSub)+1)/SubVol(tSub)/(dble(Npart)-NinSBox(tSub))/dble(SubBox-2)
!               if(fSub.eq.SubBox)then
!                  ForLog = SubVol(fSub)*(dble(Npart)-NinSBox(fSub)+1)/(BulkVolume-SubVol(fSub))/NinSBox(fSub)
!               elseif(tSub.eq.SubBox)then
!                  ForLog = (BulkVolume-SubVol(tSub))*(NinSBox(tSub)+1)/SubVol(tSub)/(dble(Npart)-NinSBox(tSub))
!               else
!                  ForLog = SubVol(fSub)*(NinSBox(tSub)+1)/SubVol(tSub)/NinSBox(fSub)
!               endif      
               deltaEnergy = deltaEnergy + log(ForLog)       
           endif
!           ------------------------------------------
!           ACCEPT THE MOVE IF deltaEnergy IS NEGATIVE
!           ------------------------------------------
            call random_number(rdn)
             if(rdn.lt.exp(-deltaEnergy))then  
                   SFenergy(i) = EnergyCNew
                   Energy   = Energy  + EnergyNew - EnergyOld 
                   CLEnergy = CLEnergy + CLEnergyNew - CLEnergyOld
                   EnergyC  = EnergyC + EnergyCNew - EnergyCOld               
                   TEnergy  = TEnergy + TEnergyNew - TEnergyOld
                    if(inMove)then
                         SucinMove(iSub) = SucinMove(iSub)+1
                    elseif(Moveout)then
                         SucMoveout(iSub) = SucMoveout(iSub)+1
                         SucMovein(tSub) =  SucMovein(tSub)+1
                         SeqNinSBox(iSub,ID)= SeqNinSBox(iSub,int(NinSBox(iSub)))
                         NinSBox(iSub) = NinSBox(iSub)-1
                         NinSBox(tSub) = NinSBox(tSub)+1
                         SeqNinSBox(tSub,int(NinSBox(tSub))) = i
                   elseif(Movein)then
                         SucMovein(iSub) = SucMovein(iSub)+1
                         SucMoveout(fSub) = SucMoveout(fSub) +1
                         SeqNinSBox(fSub,ID)= SeqNinSBox(fSub,int(NinSBox(fSub)))
                         NinSBox(fSub) = NinSBox(fSub)-1
                         NinSBox(iSub) = NinSBox(iSub)+1
                         SeqNinSBox(iSub,int(NinSBox(iSub))) = i
                   endif
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
                  enddo
             endif
             
         return
      end subroutine SecMove
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       SUBROUTINE SECADJUSTMENT
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine secAdjustment
        USE PHYSICAL_M
        USE MCSETTING_M
        USE SUBBOX_M
        implicit none
        
        integer i
        
        real*8 AcceptanceRatio
        
!       --------------------------------------------------------------
!       CHECK THE ACCEPTANCE RATIO & THEN ADJUST THE DISPLACEMENT STEP
!       --------------------------------------------------------------
        do i=1, SubBox
           if(AtpinMove(i).gt.10)then
                AcceptanceRatio = real(SucinMove(i))/real(AtpinMove(i))
           else
                !CYCLE
                AcceptanceRatio = 1.0e0
           endif
       
           if(AcceptanceRatio.lt.AccRatio)then
            SStepX(i) = SStepX(i)*0.95E0
            SStepY(i) = SStepY(i)*0.95E0
            SStepZ(i) = SStepZ(i)*0.95E0
           endif
           if(AcceptanceRatio.gt.AccRatio)then
            SStepX(i) = SStepX(i)*1.05E0
            SStepY(i) = SStepY(i)*1.05E0
            SStepZ(i) = SStepZ(i)*1.05E0
           endif
!          ----------------------------------------
!          INSURE THAT THE deltaR HAS A UPPER BOUND
!          ----------------------------------------
           if((SStepX(i) .gt. IniStepX(i)) .OR. (SStepX(i) .lt. 0.001E0))SStepX(i) = IniStepX(i)
           if((SStepY(i) .gt. IniStepY(i)) .OR. (SStepY(i) .lt. 0.001E0))SStepY(i) = IniStepY(i)
           if((SStepZ(i) .gt. IniStepZ(i)) .OR. (SStepZ(i) .lt. 0.001E0))SStepZ(i) = IniStepZ(i)
        enddo
!       ------------------
!       RESET THE COUNTERS
!       ------------------
        AtpinMove = 0
        SucinMove = 0
        
        return
        end SUBROUTINE secAdjustment
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Smcexc
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine Smcexc(iSub,TEnergy,Energy,CLEnergy,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer i,J,new,ID
       integer iSub
       
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       real*8 TEnergy, TEnergyOld, TEnergyNew
       real*8 Energy,EnergyC, CLEnergy
       real*8 deltaEnergy, rdn, rdn1
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 ForLog, arg, aimSub
       
       logical inMove, Movein, Moveout
        
!      -----------------------------------------------------------------------------
!      Equal possibility for insertion, deletion or displacement in the same SubBox
!      -----------------------------------------------------------------------------
          call random_number(rdn1)
          IF(rdn1.le.0.75e0)THEN         !Pick one particle and do the displacement or deletion
             if(NinSBox(iSub).eq.0)return 
             call random_number(rdn)
             ID  = int(rdn*NinSBox(iSub)) + 1
             if(ID.gt.NinSBox(iSub)) ID = NinSBox(iSub)
             i=SeqNinSBox(iSub,ID)  
             call random_number(rdn)
             if(rdn.le.0.5e0)then  !displace in the same SubBox
               inMove  = .TRUE.
               Moveout = .FALSE.  
               Movein  = .FALSE.
               AtpinMove(iSub) = AtpinMove(iSub) + 1
             else     ! Move out from this SubBox
               inMove  = .FALSE.
               Moveout = .TRUE. 
               Movein  = .FALSE. 
             endif
          ELSE                    !Insert one particle into this bin
               inMove  = .FALSE.
               Moveout = .FALSE. 
               Movein  = .TRUE. 
          ENDIF
!         --------------------------------------
!         OLD ENERGY OF PARTICLE i IN ADSORPTION
!         --------------------------------------
          if(inMove.OR.Moveout)then
             call energySingleParticle(1,i, EnergyOld, CLEnergyOld, EnergyCOld)
             TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld
          endif
          
          if(Moveout)then  !delete
             arg = dble(NinSBox(iSub))*exp(TEnergyOld/Temperature)/(zactivity(iPC)*SubVol(iSub))
             call random_number(rdn)
             if(rdn .lt. arg)then     ! DELETE  
               Energy   = Energy - EnergyOld 
               CLEnergy = CLEnergy - CLEnergyOld
               EnergyC  = EnergyC - EnergyCOld 
               TEnergy  = TEnergy - TEnergyOld
               SucMoveout(iSub) = SucMoveout(iSub)+1
               SeqNinSBox(iSub,ID)= SeqNinSBox(iSub,int(NinSBox(iSub)))
               NinSBox(iSub) = NinSBox(iSub)-1  
               aimSub = IndNinSBox(Npart)
               do j=1, NinSBox(int(aimSub))
                  if(SeqNinSBox(int(aimSub),j).eq.Npart)then
                     SeqNinSBox(int(aimSub),j)=i
                     EXIT
                   endif
               enddo
               IndNinSBox(i)=aimSub
               SFenergy(i) = SFenergy(Npart)
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
             endif
            return
          endif
          
!         -------------------------------------------------
!         Store the information for i before any movement
!         ------------------------------------------------- 
          if(inMove)then
             new = Npart + 1
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
             enddo
             
             call random_number(rdn)
             mcx(i) = mcx(i) + (rdn - 0.5E0)*SStepX(iSub)
             if(.NOT. pbcx)then
               if((mcx(i) .lt. 0.E0) .or. (mcx(i) .gt. CarbonLengthx1))then
                 mcx(i)=mcx(new)
                 return
                endif 
             endif
             
             call random_number(rdn)
             mcy(i) = mcy(i) + (rdn - 0.5E0)*SStepY(iSub)  
             if(.NOT. pbcy)then    
                if((mcy(i) .lt. 0.E0) .or. (mcy(i) .gt. CarbonLengthy1))then
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  return
                endif
             endif
             
             call random_number(rdn)
             mcz(i) = mcz(i) + (rdn - 0.5E0)*SStepZ(iSub)
             if(mcz(i).lt.LayerZL(iSub) .OR. mcz(i).ge.LayerZH(iSub))then
                mcx(i) = mcx(new)
                mcy(i) = mcy(new)
                mcz(i) = mcz(new)
                return
             endif
             
             if(.NOT. pbcz)then    
               if((mcz(i) .lt. 0.E0) .or. (mcz(i) .gt. BoxLengthZ))then
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  mcz(i) = mcz(new)
                  return
                endif
              endif
!            
              call positioncheck(mcx(i),mcy(i),mcz(i))  
!             ---------------------------------------
!             FOR ORIENTATIONAL MOVES    
!             ---------------------------------------
!             Put the mass center on original point
!             ---------------------------------------
               do j = 1, NLJ
                  LJxnew(j) = LJx(i,j) - MCx(new)
                  LJynew(j) = LJy(i,j) - MCy(new)
                  LJznew(j) = LJz(i,j) - MCz(new)
               enddo
               do j = 1, NCL
                  CLxnew(j) = CLx(i,j) - MCx(new)
                  CLynew(j) = CLy(i,j) - MCy(new)
                  CLznew(j) = CLz(i,j) - MCz(new)
                enddo
!             ----------
!             Rotate
!             ----------
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
              enddo     

!             -------------------------
!             NEW ENERGY OF PARTICLE I
!             -------------------------
              call  energySingleParticle(1,i, EnergyNew, CLEnergyNew, EnergyCNew)
              TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
              deltaEnergy =  (TEnergyNew - TEnergyOld)/Temperature
              call random_number(rdn)
              if(rdn.lt.exp(-deltaEnergy))then  
                 SFenergy(i) = EnergyCNew
                 Energy   = Energy  + EnergyNew - EnergyOld 
                 CLEnergy = CLEnergy + CLEnergyNew - CLEnergyOld
                 EnergyC  = EnergyC + EnergyCNew - EnergyCOld  
                 TEnergy  = TEnergy + TEnergyNew - TEnergyOld
                 SucinMove(iSub)=SucinMove(iSub)+1
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
                 enddo
              endif 
              return
            endif
         
            if(Movein)then   !Insertion
              new=Npart+1
              call random_number(rdn)
              mcx(new) = rdn*CarbonLengthx1
              call random_number(rdn)
              mcy(new) = rdn*CarbonLengthy1
              call random_number(rdn)
              mcz(new) = LayerZL(iSub) + rdn*(LayerZH(iSub)-LayerZL(iSub)) 
           
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
                    
              call energysingleparticle(1,new,EnergyNew, CLEnergyNew, EnergyCNew)
              TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
              arg = zactivity(iPC)*SubVol(iSub)*exp(-TEnergyNew/Temperature)/dble(NinSBox(iSub)+1)     
              call random_number(rdn)
              if(rdn .lt. arg)then
                SFenergy(new) = EnergyCNew
                Npart         = Npart + 1
                Energy        = Energy + EnergyNew
                CLEnergy      = CLEnergy + CLEnergyNew
                EnergyC       = EnergyC + EnergyCNew
                TEnergy       = TEnergy + TEnergyNew
                SucMovein(iSub) = SucMovein(iSub)+ 1
                NinSBox(iSub) = NinSBox(iSub)+1
                SeqNinSBox(iSub,int(NinSBox(iSub))) = Npart
                IndNinSBox(Npart) = iSub
              endif
              return
            endif
         return
      end subroutine Smcexc
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Smcexc1
!      Insertion
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine Smcexc1(iSub,TEnergy,Energy,CLEnergy,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer i,J,new,ID
       integer iSub
       
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       real*8 TEnergy,TEnergyNew
       real*8 Energy,EnergyC, CLEnergy
       real*8 rdn
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 arg, ExtTEnergyNew
       
              new=Npart+1
              call random_number(rdn)
              mcx(new) = rdn*CarbonLengthx1
              call random_number(rdn)
              mcy(new) = LayerYL(iSub) + rdn*(LayerYH(iSub)-LayerYL(iSub)) 
              call random_number(rdn)
              mcz(new) = LayerZL(iSub) + rdn*(LayerZH(iSub)-LayerZL(iSub)) 
           
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
               
              indNinSBox(new)=iSub      
              call energysingleparticle(1,new,EnergyNew, CLEnergyNew, EnergyCNew)
              TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
              arg = zactivity(iPC)*SubVol(iSub)*exp(-TEnergyNew/Temperature)/dble(NinSBox(iSub)+1) 
              if(ExtraBin.and.(iSub.eq.SubBox))then
                 ExtTEnergyNew = ExtEnergy1+ExtCLEnergy1 
                 arg = zactivity(iPC)*SubVol(iSub)*exp(-ExtTEnergyNew/Temperature)/dble(NinSBox(iSub)+1) 
              endif
              call random_number(rdn)
              if(rdn .lt. arg)then
                IF(ExtraBin.and.(iSub.eq.SubBox))THEN
                  SFenergy(new) = 0.0e0
                  ExtEnergy = ExtEnergy + ExtEnergy1
                  ExtCLEnergy = ExtCLEnergy + ExtCLEnergy1
                  ExtTEnergy = ExtTEnergy + ExtEnergy1 + ExtCLEnergy1               
                ELSE
                  SFenergy(new) = EnergyCNew
                  Energy        = Energy + EnergyNew
                  CLEnergy      = CLEnergy + CLEnergyNew
                  EnergyC       = EnergyC + EnergyCNew
                  TEnergy       = TEnergy + TEnergyNew
                ENDIF
                Npart  = Npart + 1
                SucMovein(iSub) = SucMovein(iSub)+ 1
                NinSBox(iSub) = NinSBox(iSub)+1
                SeqNinSBox(iSub,int(NinSBox(iSub))) = Npart
                IndNinSBox(Npart) = iSub
                call particleposition(new,1)
              endif
      return
      end subroutine Smcexc1
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Smcexc2
!      DELETION
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine Smcexc2(iSub,TEnergy,Energy,CLEnergy,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer i,J,ID
       integer iSub
       
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 TEnergy, TEnergyOld, ExtTEnergyOld
       real*8 Energy,EnergyC, CLEnergy
       real*8 rdn
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 arg, aimSub
       
       logical inMove, Movein, Moveout
        
             if(NinSBox(iSub).eq.0)return 
             call random_number(rdn)
             ID  = int(rdn*NinSBox(iSub)) + 1
             if(ID.gt.NinSBox(iSub)) ID = NinSBox(iSub)
             i=SeqNinSBox(iSub,ID)  
             
             call energySingleParticle(1,i, EnergyOld, CLEnergyOld, EnergyCOld)
             TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld
             arg = dble(NinSBox(iSub))*exp(TEnergyOld/Temperature)/(zactivity(iPC)*SubVol(iSub))
             if(ExtraBin.and.(iSub.eq.SubBox))then
                ExtTEnergyOld = ExtEnergy1 +  ExtCLEnergy1
                arg = dble(NinSBox(iSub))*exp(ExtTEnergyOld/Temperature)/(zactivity(iPC)*SubVol(iSub))
             endif
             call random_number(rdn)
             if(rdn .lt. arg)then     ! DELETE 
               if(ExtraBin.and.(isub.eq.SubBox))then
                  ExtEnergy = ExtEnergy - ExtEnergy1
                  ExtCLEnergy = ExtCLEnergy - ExtCLEnergy1
                  ExtTEnergy = ExtTEnergy - ExtEnergy1- ExtCLEnergy1
               else 
                  Energy   = Energy - EnergyOld 
                  CLEnergy = CLEnergy - CLEnergyOld
                  EnergyC  = EnergyC - EnergyCOld 
                  TEnergy  = TEnergy - TEnergyOld
               endif
               SucMoveout(iSub) = SucMoveout(iSub)+1
               SeqNinSBox(iSub,ID)= SeqNinSBox(iSub,int(NinSBox(iSub)))
               NinSBox(iSub) = NinSBox(iSub)-1  
               aimSub = IndNinSBox(Npart)
               do j=1, NinSBox(int(aimSub))
                  if(SeqNinSBox(int(aimSub),j).eq.Npart)then
                     SeqNinSBox(int(aimSub),j)=i
                     EXIT
                   endif
               enddo
               IndNinSBox(i)=aimSub
               SFenergy(i) = SFenergy(Npart)
               call particleposition(i,-1)
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
             endif
          
      return
      end subroutine Smcexc2

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE Smcexc3
!      DISPLACEMENT
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine Smcexc3(iSub,TEnergy,Energy,CLEnergy,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer i,J,new,ID
       integer iSub
       
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       real*8 TEnergy, TEnergyOld, TEnergyNew
       real*8 Energy,EnergyC, CLEnergy
       real*8 deltaEnergy, rdn, rdn1
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 ForLog, arg, aimSub
       real*8 ExtEnergyOld, ExtCLEnergyOld,ExtEnergyNew, ExtCLEnergyNew
       real*8 ExtTEnergyOld, ExtTEnergyNew
       
       logical inMove, Movein, Moveout
        
             if(NinSBox(iSub).eq.0)return
             AtpinMove(iSub) = AtpinMove(iSub) + 1
             call random_number(rdn)
             ID  = int(rdn*NinSBox(iSub)) + 1
             if(ID.gt.NinSBox(iSub)) ID = NinSBox(iSub)
             i=SeqNinSBox(iSub,ID)  

             call energySingleParticle(1,i, EnergyOld, CLEnergyOld, EnergyCOld)

             TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld
             if(ExtraBin.and.(iSub.eq.SubBox))then
               ExtTEnergyOld = ExtEnergy1+ExtCLEnergy1
               ExtEnergyOld = ExtEnergy1
               ExtCLEnergyOld = ExtCLEnergy1
             endif
             new = Npart + 1
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
             enddo
             
             call random_number(rdn)
             mcx(i) = mcx(i) + (rdn - 0.5E0)*SStepX(iSub)
             if(.NOT. pbcx)then
               if((mcx(i) .lt. 0.E0) .or. (mcx(i) .gt. CarbonLengthx1))then
                 mcx(i)=mcx(new)
                 return
                endif 
             endif
             
             call random_number(rdn)
             mcy(i) = mcy(i) + (rdn - 0.5E0)*SStepY(iSub)  
             if(mcy(i).lt.LayerYL(iSub) .OR. mcy(i).ge.LayerYH(iSub))then
                mcx(i) = mcx(new)
                mcy(i) = mcy(new)
                return
             endif
             if(.NOT. pbcy)then    
                if((mcy(i) .lt. 0.E0) .or. (mcy(i) .gt. CarbonLengthy1))then
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  return
                endif
             endif
             
             call random_number(rdn)
             mcz(i) = mcz(i) + (rdn - 0.5E0)*SStepZ(iSub)
             if(mcz(i).lt.LayerZL(iSub) .OR. mcz(i).ge.LayerZH(iSub))then
                mcx(i) = mcx(new)
                mcy(i) = mcy(new)
                mcz(i) = mcz(new)
                return
             endif
             if(.NOT. pbcz)then    
               if((mcz(i) .lt. 0.E0) .or. (mcz(i) .gt. BoxLengthZ))then
                  mcx(i) = mcx(new)
                  mcy(i) = mcy(new)
                  mcz(i) = mcz(new)
                  return
                endif
              endif
            
              call positioncheck(mcx(i),mcy(i),mcz(i))  
!             ---------------------------------------
!             FOR ORIENTATIONAL MOVES    
!             ---------------------------------------
!             Put the mass center on original point
!             ---------------------------------------
               do j = 1, NLJ
                  LJxnew(j) = LJx(i,j) - MCx(new)
                  LJynew(j) = LJy(i,j) - MCy(new)
                  LJznew(j) = LJz(i,j) - MCz(new)
               enddo
               do j = 1, NCL
                  CLxnew(j) = CLx(i,j) - MCx(new)
                  CLynew(j) = CLy(i,j) - MCy(new)
                  CLznew(j) = CLz(i,j) - MCz(new)
                enddo
!             ----------
!             Rotate
!             ----------
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
              enddo     

!             -------------------------
!             NEW ENERGY OF PARTICLE I
!             -------------------------             
              call  energySingleParticle(1,i, EnergyNew, CLEnergyNew, EnergyCNew)
              TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
              deltaEnergy =  (TEnergyNew - TEnergyOld)/Temperature
              if(ExtraBin.and.(iSub.eq.SubBox))then
                ExtTEnergyNew = ExtEnergy1+ExtCLEnergy1
                ExtEnergyNew = ExtEnergy1
                ExtCLEnergyNew = ExtCLEnergy1
                deltaEnergy =  (ExtTEnergyNew - ExtTEnergyOld)/Temperature
              endif

              call random_number(rdn)
              if(rdn.lt.exp(-deltaEnergy))then  
                 if(ExtraBin.and.(iSub.eq.SubBox))then
                    SFenergy(i) = 0.0e0
                    ExtEnergy   = ExtEnergy  + ExtEnergyNew - ExtEnergyOld 
                    ExtCLEnergy = ExtCLEnergy + ExtCLEnergyNew - ExtCLEnergyOld
                    ExtTEnergy  = ExtTEnergy + ExtTEnergyNew - ExtTEnergyOld 
                 else
                    SFenergy(i) = EnergyCNew
                    Energy   = Energy  + EnergyNew - EnergyOld 
                    CLEnergy = CLEnergy + CLEnergyNew - CLEnergyOld
                    EnergyC  = EnergyC + EnergyCNew - EnergyCOld  
                    TEnergy  = TEnergy + TEnergyNew - TEnergyOld
                 endif   
                 SucinMove(iSub)=SucinMove(iSub)+1
                 call particleposition(new,-1)
                 call particleposition(i,1)
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
                 enddo
              endif 
         
    return
    end subroutine Smcexc3      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE SecMOVE1
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine SecMove1(fSub,tSub,TEnergy,Energy,CLEnergy,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer i,J,new,ID, flag, repeatN, new1
       integer iSub1,iSub2, fSub,tSub,dSub
       integer iSwap, iDisp, k, FLAG1, FLAG2
       
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       real*8 TEnergy, TEnergyOld, TEnergyNew
       real*8 Energy,EnergyC, CLEnergy
       real*8 deltaEnergy, rdn, rdn1
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 ForLog, Crit, VirialOld, VirialNew
       real*8 ExtTEnergyOld,ExtEnergyOld, ExtCLEnergyOld
       real*8 ExtTEnergyNew,ExtEnergyNew, ExtCLEnergyNew 
       real*8 VirialNewF, VirialOldF
          

          if(NinSBox(fSub).eq.0)return
          call random_number(rdn) 
          ID  = int(rdn*NinSBox(fSub)) + 1
          if(ID.gt.NinSBox(fSub)) ID = NinSBox(fSub)
          i=SeqNinSBox(fSub,ID)  
          
          call energySingleParticle(1,i, EnergyOld, CLEnergyOld, EnergyCOld )
          TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld
          VirialOld = Virial + VirialCL
          VirialOldF = VirialF
          if(ExtraBin .and. (fSub.eq.SubBox))then
             ExtTEnergyOld=ExtEnergy1+ExtCLEnergy1
             ExtEnergyOld = ExtEnergy1
             ExtCLEnergyOld = ExtCLEnergy1
          endif
          new=Npart+1
          indNinSBox(new)=tSub
          call random_number(rdn)
          MCx(new) = rdn*CarbonLengthx1
          call random_number(rdn)
          MCy(new) = rdn*CarbonLengthy1
          call random_number(rdn)
          MCz(new) = LayerZL(tSub) + rdn*(LayerZH(tSub)-LayerZL(tSub)) 
 
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
                    
          call energysingleparticle(1,new,EnergyNew, CLEnergyNew, EnergyCNew)
          TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
          Crit = dble(NinSBox(fSub))*SubVol(tSub)/SubVol(fSub)/dble(NinSBox(tSub)+1)*exp(-(TEnergyNew-TEnergyOld)/Temperature)
          if(ExtraBin .and. (tSub.eq.SubBox))then
             ExtTEnergyNew=ExtEnergy1+ExtCLEnergy1
             ExtEnergyNew = ExtEnergy1
             ExtCLEnergyNew = ExtCLEnergy1
             Crit = dble(NinSBox(fSub))*SubVol(tSub)/SubVol(fSub)/dble(NinSBox(tSub)+1)*exp(-(ExtTEnergyNew-TEnergyOld)/Temperature)
!             write(62,*)'7573 tSub=SubBox,tSub, Crit',tSub, Crit
!             read(62,*)
          endif
          if(ExtraBin .and. (fSub.eq.SubBox))then
             Crit = dble(NinSBox(fSub))*SubVol(tSub)/SubVol(fSub)/dble(NinSBox(tSub)+1)*exp(-(TEnergyNew-ExtTEnergyOld)/Temperature)
          endif
              call random_number(rdn)
              if(rdn.lt.Crit)then
                 SucMoveout(fSub) = SucMoveout(fSub)+1
                 SucMovein(tSub) =  SucMovein(tSub)+1
                 SeqNinSBox(fSub,ID)= SeqNinSBox(fSub,int(NinSBox(fSub)))
                 NinSBox(fSub) = NinSBox(fSub)-1
                 NinSBox(tSub) = NinSBox(tSub)+1
                 SeqNinSBox(tSub,int(NinSBox(tSub))) = i
                 indNinSBox(i)=tSub
                 call particleposition(i,-1)
                 call particleposition(new,1)
                 IF(.NOT. ExtraBin)THEN
                    if(fSub.eq.SubBox)then
                       call EBulkChange(1,i,new)
                    elseif(tSub.eq.SubBox)then
                       call EBulkChange(2,i,new)
                    else
                       call EBulkChange(4,i,new)
                    endif
                 ENDIF
                 MCx(i) = MCx(new)
                 MCy(i) = MCy(new)
                 MCz(i) = MCz(new)
                 do j = 1, NLJ
                    ljx(i,j) = ljx(new,j)
                    ljy(i,j) = ljy(new,j)
                    ljz(i,j) = ljz(new,j)
                 enddo
                 do j = 1, NCL
                    CLx(i,j) = CLx(new,j)
                    CLy(i,j) = CLy(new,j)
                    CLz(i,j) = CLz(new,j)
                 enddo
                 if(ExtraBin .and. (fSub.eq.SubBox))then
                    SFenergy(i)  = EnergyCNew
                    Energy       = Energy  + EnergyNew 
                    CLEnergy     = CLEnergy + CLEnergyNew
                    EnergyC      = EnergyC + EnergyCNew               
                    TEnergy      = TEnergy + TEnergyNew  
                    ExtEnergy    = ExtEnergy  - ExtEnergyOld 
                    ExtCLEnergy  = ExtCLEnergy - ExtCLEnergyOld
                    ExtTEnergy   = ExtTEnergy - ExtTEnergyOld  
                 elseif(ExtraBin .and. (tSub.eq.SubBox))then 
                    SFenergy(i)  = 0.0E0
                    Energy       = Energy  - EnergyOld 
                    CLEnergy     = CLEnergy - CLEnergyOld
                    EnergyC      = EnergyC - EnergyCOld               
                    TEnergy      = TEnergy - TEnergyOld  
                    ExtEnergy    = ExtEnergy  + ExtEnergyNew 
                    ExtCLEnergy  = ExtCLEnergy + ExtCLEnergyNew
                    ExtTEnergy   = ExtTEnergy + ExtTEnergyNew
                 else
                    call energysingleparticle(1,i,EnergyNew, CLEnergyNew, EnergyCNew)
                    TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
                    SFenergy(i) = EnergyCNew
                    Energy   = Energy  + EnergyNew - EnergyOld 
                    CLEnergy = CLEnergy + CLEnergyNew - CLEnergyOld
                    EnergyC  = EnergyC + EnergyCNew - EnergyCOld               
                    TEnergy  = TEnergy + TEnergyNew - TEnergyOld 
                    VirialNew = Virial + VirialCL
                    TVirial = TVirial + VirialNew - VirialOld
                    VirialNewF = VirialF
                      IF( (fSub.eq.SubBox) .AND. (tSub.ne.SubBox) )THEN
                         BulkVirial = BulkVirial - VirialOldF - (VirialOld-VirialOldF)/2.0 + VirialNewF/2.0
                      ELSEIF( (fSub.eq.SubBox) .AND. (tSub.eq.SubBox) )THEN
                         BulkVirial = BulkVirial - VirialOldF-(VirialOld-VirialOldF)/2.0+VirialNewF+(VirialNew-VirialNewF)/2.0
                      ELSEIF( (fSub.ne.SubBox) .AND. (tSub.eq.SubBox) )THEN
                         BulkVirial = BulkVirial - VirialOldF/2.0 + VirialNewF + (VirialNew-VirialNewF)/2.0   
                      ELSEIF( (fSub.ne.SubBox) .AND. (tSub.ne.SubBox) )THEN
                         BulkVirial = BulkVirial - VirialOldF/2.0 + VirialNewF/2.0
                      ENDIF
                 endif          
              endif

      return
      end subroutine SecMove1

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE SecMCMOVE2
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine SecMove2(dSub,TEnergy,Energy,CLEnergy,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE FLUCTUATION_M
       USE POREFIGURE_M
       USE SUBBOX_M
       implicit none
       
       integer i,J,new,ID 
       integer dSub
       
       real*8 EnergyOld, CLEnergyOld, EnergyCOld
       real*8 EnergyNew, CLEnergyNew, EnergyCNew
       real*8 TEnergy, TEnergyOld, TEnergyNew
       real*8 Energy,EnergyC, CLEnergy
       real*8 deltaEnergy, rdn
       real*8 LJxnew(20),LJynew(20),LJznew(20)
       real*8 CLxnew(20),CLynew(20),CLznew(20)
       real*8 ExtTEnergyOld, ExtEnergyOld, ExtCLEnergyOld
       real*8 ExtTEnergyNew, ExtEnergyNew, ExtCLEnergyNew
       real*8 VirialNew, VirialOld, VirialNewF, VirialOldF

              
                 AtpinMove(dSub)=AtpinMove(dSub)+1
                 call random_number(rdn) 
                 ID  = int(rdn*NinSBox(dSub)) + 1
                 if(ID.gt.NinSBox(dSub)) ID = NinSBox(dSub)
                 i=SeqNinSBox(dSub,ID)  
                 call energySingleParticle(1,i, EnergyOld, CLEnergyOld, EnergyCOld )
                 TEnergyOld = EnergyOld + CLEnergyOld + EnergyCOld
                 VirialOld  = Virial + VirialCL
                 VirialOldF = VirialF
                 if(ExtraBin .and. (dSub.eq.SubBox))then
                   ExtTEnergyOld = ExtEnergy1 + ExtCLEnergy1
                   ExtEnergyOld = ExtEnergy1
                   ExtCLEnergyOld = ExtCLEnergy1
                 endif

                 new = Npart + 1
                 MCx(new) = MCx(i)
                 MCy(new) = MCy(i) 
                 MCz(new) = MCz(i)
                 do j = 1, NLJ
                    LJx(new,j) = LJx(i,j)
                    LJy(new,j) = LJy(i,j)
                    LJz(new,j) = LJz(i,j)
                 enddo
                 do j = 1, NCL
                    CLx(new,j) = CLx(i,j)
                    CLy(new,j) = CLy(i,j)
                    CLz(new,j) = CLz(i,j)
                 enddo
                 
                 call random_number(rdn)
                 MCx(i) = MCx(i) + (rdn - 0.5E0)*SStepX(dSub)
                 if(.NOT. pbcx)then
                   if((MCx(i) .lt. 0.E0) .or. (MCx(i) .gt. CarbonLengthx1))then
                      MCx(i)=MCx(new)
                      return
                   endif 
                 endif
       
                 call random_number(rdn)
                 MCy(i) = MCy(i) + (rdn - 0.5E0)*SStepY(dSub)  
                 if(.NOT. pbcy)then    
                   if((MCy(i) .lt. 0.E0) .or. (MCy(i) .gt. CarbonLengthy1))then
                      MCx(i) = MCx(new)
                      MCy(i) = MCy(new)
                      return
                   endif
                 endif 
                  
                 call random_number(rdn)
                 MCz(i) = MCz(i) + (rdn - 0.5E0)*SStepZ(dSub)
                 IF(ExtraBin .AND. (dSub.eq.SubBox))THEN
                      if ( MCz(i).gt.LayerZH(dSub)) then                     
                        MCz(i) = MCz(i) - int((MCz(i)-LayerZL(dSub))/(LayerZH(dSub)-LayerZL(dSub)))*(LayerZH(dSub)-LayerZL(dSub))
                      elseif (MCz(i).lt.LayerZL(dSub)) then
                        MCz(i) = MCz(i) - int((MCz(i)-LayerZL(dSub))/(LayerZH(dSub)-LayerZL(dSub)))*(LayerZH(dSub)-LayerZL(dSub)) + (LayerZH(dSub)-LayerZL(dSub))
                      endif
                 ELSE
                    if(MCz(i).lt.LayerZL(dSub) .OR. MCz(i).ge.LayerZH(dSub))then
                       MCx(i) = MCx(new)
                       MCy(i) = MCy(new)
                       MCz(i) = MCz(new)
                       return
                    endif
                    if(.NOT. pbcz)then    
                      if((MCz(i) .lt. 0.E0) .or. (MCz(i) .gt. BoxLengthZ))then
                         MCx(i) = MCx(new)
                         MCy(i) = MCy(new)
                         MCz(i) = MCz(new)
                         return
                      endif
                    endif
                 ENDIF
                 call positioncheck(MCx(i),MCy(i),MCz(i))  
!             ---------------------------------------
!              FOR ORIENTATIONAL MOVES    
!             ---------------------------------------
!              Put the mass center on original point
!             ---------------------------------------
                do j = 1, NLJ
                   LJxnew(j) = LJx(i,j) - MCx(new)
                   LJynew(j) = LJy(i,j) - MCy(new)
                   LJznew(j) = LJz(i,j) - MCz(new)
                enddo
                do j = 1, NCL
                   CLxnew(j) = CLx(i,j) - MCx(new)
                   CLynew(j) = CLy(i,j) - MCy(new)
                   CLznew(j) = CLz(i,j) - MCz(new)
                enddo
!             ----------
!              Rotate
!             ----------
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
                enddo     
               call  energySingleParticle(1,i, EnergyNew, CLEnergyNew, EnergyCNew)
               TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
               deltaEnergy =  (TEnergyNew - TEnergyOld)/Temperature
               VirialNew = Virial + VirialCL
               VirialNewF = VirialF
               if(ExtraBin .and. (dSub.eq.SubBox))then
                   ExtTEnergyNew = ExtEnergy1 + ExtCLEnergy1
                   ExtEnergyNew = ExtEnergy1
                   ExtCLEnergyNew = ExtCLEnergy1
                   deltaEnergy =  (ExtTEnergyNew - ExtTEnergyOld)/Temperature
               endif

               call random_number(rdn)
               if(rdn.lt.exp(-deltaEnergy))then 
                   if(ExtraBin .and. (dSub.eq.SubBox))then
                      ExtEnergy   = ExtEnergy  + ExtEnergyNew - ExtEnergyOld 
                      ExtCLEnergy = ExtCLEnergy + ExtCLEnergyNew - ExtCLEnergyOld
                      ExtTEnergy  = ExtTEnergy + ExtTEnergyNew - ExtTEnergyOld
                   else 
                      SFenergy(i) = EnergyCNew
                      Energy   = Energy  + EnergyNew - EnergyOld 
                      CLEnergy = CLEnergy + CLEnergyNew - CLEnergyOld
                      EnergyC  = EnergyC + EnergyCNew - EnergyCOld               
                      TEnergy  = TEnergy + TEnergyNew - TEnergyOld 
                      TVirial = TVirial + VirialNew - VirialOld
                      IF(dSub.eq.SubBox)THEN   
                         call EBulkChange(3,new,i)
                      ELSE
                         call EBulkChange(4,new,i)
                      ENDIF  
                      IF(dSub.eq.SubBox)THEN
                         BulkVirial = BulkVirial - VirialOldF-(VirialOld-VirialOldF)/2.0+VirialNewF+(VirialNew-VirialNewF)/2.0 
                      ELSEIF(dSub.ne.SubBox)THEN 
                         BulkVirial = BulkVirial -VirialOldF/2.0 + VirialNewF/2.0
                      ENDIF             
                   endif
                   SucinMove(dSub) = SucinMove(dSub)+1
                   call particleposition(new,-1)
                   call particleposition(i,1)
               else             
                   MCx(i) = MCx(new)
                   MCy(i) = MCy(new)
                   MCz(i) = MCz(new)
                   do j = 1, NLJ
                     ljx(i,j) = ljx(new,j)
                     ljy(i,j) = ljy(new,j)
                     ljz(i,j) = ljz(new,j)
                   enddo
                   do j = 1, NCL
                     CLx(i,j) = CLx(new,j)
                     CLy(i,j) = CLy(new,j)
                     CLz(i,j) = CLz(new,j)
                   enddo
                endif
      return
      end subroutine SecMove2


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE EBulkChange
!      THIS SUBROUTINE CALCULATES THE ENERGY CHANGE IN THE BULK PHASE PART
!      FOR CMC(NVT) AND MU-CMC ENSEMBLES
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE EBulkChange(switch,iold, inew)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE BUCKINGHAM_M
       USE FLUCTUATION_M 
       USE SUBBOX_M 
       implicit none
       
       integer j,k
       integer switch, iold, inew
       
       real*8  U, UC, CLU
       real*8  PEnergy,PCLEnergy
       real*8  minusE,addE

!      --------------------------------
!      INITIALISE THE ENERGY 
!      --------------------------------
       minusE = 0.0e0
       addE   = 0.0e0

!      ---------------------------------------------------
!      ---------------------------------------------------
       IF(switch .eq. 1)THEN    ! Move out from Bulk bin
           if(steele)then
              call steelepotential(iold,UC)
                  minusE = minusE + UC 
           endif                   
           
           do j=1, Npart 
              if((iold .ne. j).and.(inew.ne.j))then 
                      call PairFFEnergy(iold,j,PEnergy,PCLEnergy)
                      if(indNinSBox(j).EQ.SubBox)then
                         minusE = minusE + PEnergy + PCLEnergy
                      else
                         minusE = minusE + (PEnergy + PCLEnergy)/2.0
                      endif
                      
                      call PairFFEnergy(inew,j,PEnergy,PCLEnergy)  
                      if(indNinSBox(j).eq. SubBox)then 
                         addE = addE + (PEnergy + PCLEnergy)/2.0
                      endif
              endif              
           enddo
       ENDIF
       
       
       IF(switch .eq. 2)THEN  ! Move into Bulk bin
          if(steele)then
             call steelepotential(inew,UC)
                addE   = addE + UC
          endif  
          do j=1, Npart
               if((iold .ne. j).and.(inew.ne.j))then
                      call PairFFEnergy(iold,j,PEnergy,PCLEnergy)
                      if(indNinSBox(j).EQ.SubBox)then
                           minusE = minusE + (PEnergy + PCLEnergy)/2.0
                      endif
                      
                      call PairFFEnergy(inew,j,PEnergy,PCLEnergy)
                      if(indNinSBox(j).eq. SubBox)then  
                         addE = addE + PEnergy + PCLEnergy
                      else
                         addE = addE + (PEnergy + PCLEnergy)/2.0
                      endif
               endif
          enddo
       ENDIF
       
       IF(switch .eq. 3)THEN   !Displace in the Bulk bin
           if(steele)then
             call steelepotential(iold,UC)
                  minusE = minusE + UC 
             call steelepotential(inew,UC)
                  addE   = addE + UC
           endif      
           
           do j=1, Npart
               if((iold .ne. j).and.(inew.ne.j))then
                    call PairFFEnergy(iold,j,PEnergy,PCLEnergy)
                    if(indNinSBox(j).EQ.SubBox)then
                       minusE = minusE + PEnergy + PCLEnergy
                    else
                       minusE = minusE + (PEnergy + PCLEnergy)/2.0
                    endif
                    
                    call PairFFEnergy(inew,j,PEnergy,PCLEnergy) 
                    if(indNinSBox(j).EQ.SubBox)then 
                       addE = addE + PEnergy + PCLEnergy
                    else
                       addE = addE + (PEnergy + PCLEnergy)/2.0
                    endif
               endif
            enddo
       ENDIF
       
       IF(switch .eq. 4)THEN   !Displace in outside of the Bulk bin           
           do k=1, NinSBox(SubBox)
               j = seqNinSBox(SubBox,k)
               if((iold .ne. j).and.(inew.ne.j))then
                    call PairFFEnergy(iold,j,PEnergy,PCLEnergy)
                    minusE = minusE + (PEnergy + PCLEnergy)/2.0
                    
                    call PairFFEnergy(inew,j,PEnergy,PCLEnergy) 
                    addE = addE + (PEnergy + PCLEnergy)/2.0
               endif
            enddo
       ENDIF
       
       EinBulk = EinBulk - minusE + addE  
     return
     END SUBROUTINE EBulkChange 
                    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE PairFFEnergy
!      THIS SUBROUTINE CALCULATES THE ENERGY BETWEEN TWO FLUID PARTICLES
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE PairFFEnergy(n1,n2,PEnergy,PCLEnergy)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE BUCKINGHAM_M
       USE FLUCTUATION_M 
       implicit none
       
       integer n1,n2, k,m
       
       real*8  PEnergy,PCLEnergy
       real*8  xmclj,ymclj,zmclj, xmcn,ymcn,zmcn, mcd2
       real*8  sigma,welldepth
       real*8  xlj, ylj, zlj, d2,U, V
       real*8  xcl, ycl, zcl, cld, CLU, CLV
         
                PEnergy   = 0.0e0
                PCLEnergy = 0.0e0
                xmclj = abs(mcx(n1)-mcx(n2))
                ymclj = abs(mcy(n1)-mcy(n2))
                zmclj = abs(mcz(n1)-mcz(n2))
                xmcn = xmclj
                ymcn = ymclj
                zmcn = zmclj
                if(pbcx)then                
                  if(xmclj > CarbonLengthx1/2.0D0) xmcn = xmclj - CarbonLengthx1
                endif
                if(pbcy)then
                  if(ymclj > CarbonLengthy1/2.0D0) ymcn = ymclj - CarbonLengthy1
                endif
                if(pbcz)then
                   if(zmclj > BoxLengthZ/2.0D0) zmcn = zmclj - BoxLengthZ
                endif
                mcd2 = xmcn*xmcn + ymcn*ymcn + zmcn*zmcn
                if(mcd2.le.r2Cutoff)then !! 1-if
                  do k=1,NLJ
                     do m = 1, NLJ
                        sigma = (sigmaFF(k) + sigmaFF(m))/2.0e0
                        welldepth = sqrt(welldepthFF(k)*welldepthFF(m))
!                       --------------------
!                       DISTANCES IN X, Y, Z
!                       --------------------
                        xlj    = abs(ljx(n1,m)-ljx(n2,k))
                        ylj    = abs(ljy(n1,m)-ljy(n2,k))  
                        zlj    = abs(ljz(n1,m)-ljz(n2,k))   
!                       ----------------------------------------------------
!                        ENSURE THAT WE DEAL WITH THE NEAREST PERIODIC IMAGES
!                        ----------------------------------------------------
                        if(pbcx)then                
                           if(xmclj > CarbonLengthx1/2.0D0) xlj = xlj - CarbonLengthx1
                        endif
                        if(pbcy)then
                          if(ymclj > CarbonLengthy1/2.0D0) ylj = ylj - CarbonLengthy1
                        endif
                        if(pbcz)then
                          if(zmclj > BoxLengthZ/2.0D0) zlj = zlj - BoxLengthZ
                        endif               
!                      -----------------------
!                      INTER-PARTICLE DISTANCE
!                      -----------------------
                       d2 = xlj*xlj + ylj*ylj + zlj*zlj
!                      ----------------
!                      POTENTIAL ENERGY
!                      ----------------
                       call potentialEnergy(d2, sigma, welldepth,U,V)
                          PEnergy = PEnergy + U
                     enddo
                   enddo
                   do k = 1, NCL
                      do m = 1, NCL
                         xcl =  abs(CLx(n1,m) - CLx(n2,k))
                         ycl =  abs(CLy(n1,m) - CLy(n2,k))
                         zcl =  abs(CLz(n1,m) - CLz(n2,k))
                         if(pbcx)then                
                           if(xmclj > CarbonLengthx1/2.0D0) xcl = xcl - CarbonLengthx1
                         endif
                         if(pbcy)then
                           if(ymclj > CarbonLengthy1/2.0D0) ycl = ycl - CarbonLengthy1
                         endif
                         if(pbcz)then
                           if(zmclj > BoxLengthZ/2.0D0) zcl = zcl - BoxLengthZ
                         endif
                         cld = sqrt(xcl*xcl + ycl*ycl + zcl*zcl)
                         call coulombforce(cld, charges(k),charges(m),CLU,CLV)
                         PCLEnergy = PCLEnergy + CLU
                      enddo
                  enddo      
               endif
       RETURN
       END SUBROUTINE PairFFEnergy

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE SHAPECHANGE
!      Adjust the shape of the cross section area (XZ)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     SUBROUTINE ShapeChange(TEnergyOld, EnergyOld, CLEnergyOld,EnergyCOld,TEnergyNew, EnergyNew, CLEnergyNew,EnergyCNew)
     USE Constant_M
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE SUBBOX_M
     implicit none
     
     integer i,j,iSub
     
     real*8 rdn1,rdn2, rdn
     real*8 maxChangeL, NewLinX, NewLinY, NewLinZ
     real*8 ChangeRatio1, ChangeRatio2
     real*8 TEnergyOld, EnergyOld, CLEnergyOld, EnergyCOld
     real*8 TEnergyNew, EnergyNew, CLEnergyNew, EnergyCNew
     real*8 Energy, CLEnergy,EnergyC, SFenergyOld(100000)
     
           
           maxChangeL = 1.0  ! ONE COLLOSION DIAMETER
           call random_number(rdn1)
           IF(DecayESF)rdn1 = 0.0
           IF(rdn1.LE.0.5)THEN  ! CHANGE THE DIMENTION IN X- DIRECTION
              call random_number(rdn2)
              NewLinX = IniCarbonLengthx1  + (rdn2-0.5)*maxChangeL
              ChangeRatio1 = NewLinX/CarbonLengthx1
              NewLinZ = BoxLengthZ*CarbonLengthx1/NewLinX ! Keep the volume constant
              ChangeRatio2 = NewLinZ/BoxLengthZ
              
              DO i=1, Npart
                 MCx(i) = MCx(i)*ChangeRatio1
                 MCz(i) = MCz(i)*ChangeRatio2
                 do j = 1, NLJ
                    LJx(i,j) = LJx(i,j) + MCx(i) - MCx(i)/ChangeRatio1
                    LJz(i,j) = LJz(i,j) + MCz(i) - MCz(i)/ChangeRatio2
                 enddo

                 do j = 1, NCL
                    CLx(i,j) = CLx(i,j) + MCx(i) - MCx(i)/ChangeRatio1
                    CLz(i,j) = CLz(i,j) + MCz(i) - MCz(i)/ChangeRatio2
                 enddo
              ENDDO
              
              CarbonLengthx1 = NewLinX
              BoxLengthZ = NewLinZ
              call ConfigEnergy(EnergyNew,CLEnergyNew,EnergyCNew)
              TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew         
              
              call random_number(rdn)
              if( rdn .lt. exp(-(TEnergyNew-TEnergyOld)/Temperature) )then
                  do iSub = 1, SubBox
                    LayerZL(iSub) = LayerZL(iSub)*ChangeRatio2
                    LayerZH(iSub) = LayerZH(iSub)*ChangeRatio2
                    SStepX(iSub)  = SStepX(iSub)*ChangeRatio1
                    SStepZ(iSub)  = SStepZ(iSub)*ChangeRatio2
                    SubVol(iSub)  = CarbonLengthx1*CarbonLengthy1*(LayerZH(iSub)-LayerZL(iSub))
                    IniStepX(iSub) = CarbonLengthx1/2.0E0
                    IniStepZ(iSub) = (LayerZH(iSub)-LayerZL(iSub))/2.0e0
                  enddo
                  IF(NVT)AccVZH = AccVZH*ChangeRatio2
              else
                  CarbonLengthx1 = CarbonLengthx1/ChangeRatio1
                  BoxLengthZ = BoxLengthZ/ChangeRatio2
                  
                  DO i=1, Npart
                     MCx(i) = MCx(i)/ChangeRatio1
                     MCz(i) = MCz(i)/ChangeRatio2
                    do j = 1, NLJ
                        LJx(i,j) = LJx(i,j) + MCx(i) - MCx(i)*ChangeRatio1
                        LJz(i,j) = LJz(i,j) + MCz(i) - MCz(i)*ChangeRatio2
                    enddo

                    do j = 1, NCL
                        CLx(i,j) = CLx(i,j) + MCx(i) - MCx(i)*ChangeRatio1
                        CLz(i,j) = CLz(i,j) + MCz(i) - MCz(i)*ChangeRatio2
                    enddo
                  ENDDO
                  call ConfigEnergy(EnergyNew,CLEnergyNew,EnergyCNew)
                  TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew  
              endif
                             
           ELSE  ! CHANGE THE DIMENTION IN Y- DIRECTION
         
              call random_number(rdn2)
              NewLinY = IniCarbonLengthy1  + (rdn2-0.5)*maxChangeL
              ChangeRatio1 = NewLinY/CarbonLengthy1
              NewLinZ = BoxLengthZ*CarbonLengthy1/NewLinY ! Keep the volume constant
              ChangeRatio2 = NewLinZ/BoxLengthZ
              
              DO i=1, Npart
                 MCy(i) = MCy(i)*ChangeRatio1
                 MCz(i) = MCz(i)*ChangeRatio2
                 do j = 1, NLJ
                    LJy(i,j) = LJy(i,j) + MCy(i) - MCy(i)/ChangeRatio1
                    LJz(i,j) = LJz(i,j) + MCz(i) - MCz(i)/ChangeRatio2
                 enddo

                 do j = 1, NCL
                    CLy(i,j) = CLy(i,j) + MCy(i) - MCy(i)/ChangeRatio1
                    CLz(i,j) = CLz(i,j) + MCz(i) - MCz(i)/ChangeRatio2    
                 enddo
              ENDDO
              
              BoxLengthZ = NewLinZ
              CarbonLengthy1 = NewLinY
              call ConfigEnergy(EnergyNew,CLEnergyNew,EnergyCNew)
              TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
              
              call random_number(rdn)
              if( rdn .lt. exp(-(TEnergyNew-TEnergyOld)/Temperature) )then
                  do iSub = 1, SubBox
                    LayerZL(iSub) = LayerZL(iSub)*ChangeRatio2
                    LayerZH(iSub) = LayerZH(iSub)*ChangeRatio2
                    SStepZ(iSub)  = SStepZ(iSub)*ChangeRatio2
                    SStepY(iSub)  = SStepY(iSub)*ChangeRatio1
                    SubVol(iSub)  = CarbonLengthx1*CarbonLengthy1*(LayerZH(iSub)-LayerZL(iSub))
                    IniStepY(iSub) = CarbonLengthy1/2.0E0
                    IniStepZ(iSub) = (LayerZH(iSub)-LayerZL(iSub))/2.0e0
                  enddo
                  IF(NVT)AccVZH = AccVZH*ChangeRatio2
                  IF(GCMC)THEN
                     AccVYL        = AccVYL*ChangeRatio1
                     AccVYH        = AccVYH*ChangeRatio1
                  ENDIF
              else
                  CarbonLengthy1 = CarbonLengthy1/ChangeRatio1
                  BoxLengthZ = BoxLengthZ/ChangeRatio2
                  
                  DO i=1, Npart
                     MCy(i) = MCy(i)/ChangeRatio1
                     MCz(i) = MCz(i)/ChangeRatio2
                     SFEnergy(i) = SFEnergyOld(i)
                    do j = 1, NLJ
                        LJy(i,j) = LJy(i,j) + MCy(i) - MCy(i)*ChangeRatio1
                        LJz(i,j) = LJz(i,j) + MCz(i) - MCz(i)*ChangeRatio2
                    enddo

                    do j = 1, NCL
                        CLy(i,j) = CLy(i,j) + MCy(i) - MCy(i)*ChangeRatio1
                        CLz(i,j) = CLz(i,j) + MCz(i) - MCz(i)*ChangeRatio2
                    enddo
                  ENDDO
                  call ConfigEnergy(EnergyNew,CLEnergyNew,EnergyCNew)
                  TEnergyNew = EnergyNew + CLEnergyNew + EnergyCNew
             endif
           ENDIF
      RETURN
      ENDSUBROUTINE ShapeChange
      
      
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE BOXRESET
!      RESET THE SIMULATION BOX WITH THE AVERAGE VALUES
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     SUBROUTINE BOXRESET(SumLX, SumLY, SumLZ,NChange)
     USE Constant_M
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     USE SUBBOX_M
     implicit none
     
     integer NChange, i, j, iSub
     
     real*8 SumLX, SumLY, SumLZ
     real*8 Ave_LX, Ave_LY, Ave_LZ
     real*8 ChangeRatioX, ChangeRatioY, ChangeRatioZ
     
               Ave_LX = SumLX/dble(NChange)
               Ave_LY = SumLY/dble(NChange)
               Ave_LZ = IniCarbonLengthx1*IniCarbonLengthy1*IniBoxlengthZ/Ave_LX/Ave_LY      !!! SumLZ/dble(NChange)
               
               ChangeRatioX = Ave_LX/CarbonLengthx1
               ChangeRatioY = Ave_LY/CarbonLengthy1
               ChangeRatioZ = Ave_LZ/BoxlengthZ
               
               CarbonLengthx1 = Ave_LX
               CarbonLengthy1 = Ave_LY
               BoxlengthZ     = Ave_LZ

               DO i=1, Npart
                 MCx(i) = MCx(i)*ChangeRatioX
                 MCy(i) = MCy(i)*ChangeRatioY
                 MCz(i) = MCz(i)*ChangeRatioZ
                 
                 do j = 1, NLJ
                    LJx(i,j) = LJx(i,j) + MCx(i) - MCx(i)/ChangeRatioX
                    LJy(i,j) = LJy(i,j) + MCy(i) - MCy(i)/ChangeRatioY
                    LJz(i,j) = LJz(i,j) + MCz(i) - MCz(i)/ChangeRatioZ
                 enddo

                 do j = 1, NCL
                    CLx(i,j) = CLx(i,j) + MCx(i) - MCx(i)/ChangeRatioX
                    CLy(i,j) = CLy(i,j) + MCy(i) - MCy(i)/ChangeRatioY
                    CLz(i,j) = CLz(i,j) + MCz(i) - MCz(i)/ChangeRatioZ
                 enddo
              ENDDO

              do iSub = 1, SubBox
                    LayerZL(iSub) = LayerZL(iSub)*ChangeRatioZ
                    LayerZH(iSub) = LayerZH(iSub)*ChangeRatioZ
                    SStepX(iSub)  = SStepX(iSub)*ChangeRatioX
                    SStepY(iSub)  = SStepY(iSub)*ChangeRatioY
                    SStepZ(iSub)  = SStepZ(iSub)*ChangeRatioZ
                    SubVol(iSub)  = CarbonLengthx1*CarbonLengthy1*(LayerZH(iSub)-LayerZL(iSub))
              enddo
              IF(NVT)THEN
                 AccVZH = AccVZH*ChangeRatioZ
                 AccVYL    = AccVYL*ChangeRatioY
                 AccVYH    = AccVYH*ChangeRatioY
              ENDIF
              IF(GCMC)THEN
                 AccVYL    = AccVYL*ChangeRatioY
                 AccVYH    = AccVYH*ChangeRatioY
              ENDIF
              
              TotSurfaceArea  = CarbonLengthx1*CarbonLengthy1              ! for Surface
              BulkVolume      = CarbonLengthx1*CarbonLengthy1*BoxLengthZ   ! the volume of the simulation box 
     
              if(DecayESF)then
                  TotSurfaceArea  = CarbonLengthx1*(CollectHY-CollectLY)  
                  CollectBulkV =  CarbonLengthx1*(CollectHY-CollectLY)*BoxLengthZ 
              endif
              
              if(Nsteele .gt. 1)then
                  TotSurfaceArea  = Nsteele*CarbonLengthx1*CarbonLengthy1
              endif
     
              if(BojanSlit)then
                 TotSurfaceArea  = 2.0e0*CarbonLengthx1*CarbonLengthy1
              endif
              
               call accessiblevolume
               AccVolumeM3    = AccVolume*(SCALELENGTH*1.0E-10)**3
      
        RETURN
        ENDSUBROUTINE BOXRESET
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE decayfcalculation
!      Calculate the parameter f for the decayed SF potential
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
       pure SUBROUTINE decayfcalculation(i,decayf)
       USE FLUID_SOLID_M
       USE MCSETTING_M
       IMPLICIT NONE

       intent(in) :: i
       intent(out) :: decayf
       integer i
       real*8 decayf0, decayf
       
       decayf0 = abs(4.0d0*MCy(i)/CarbonLengthy1-2)
       
!       IF(decayf0.LE.1.0)THEN
!         decayf = 1.0
!       ELSE
!         decayf = 2.0*decayf0**3.0 - 9.0*decayf0**2.0 + 12.0*decayf0 - 4.0
!       ENDIF
       IF(decayf0.LE.1.0d0)THEN
         decayf = 1.0d0
       ELSE
         decayf = 2.0d0*decayf0**3 - 9.0d0*decayf0**2 + 12.0d0*decayf0 - 4.0d0
       ENDIF
       
       RETURN
       ENDSUBROUTINE decayfcalculation
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE ENERGYSINGLEPARTICLESF
!      THIS SUBROUTINE CALCULATES THE ENERGY OF ONE PARTICLE WITH THE SOLID 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine energysingleparticleSF(i,EnergyC)
       USE FLUID_SOLID_M
       USE PHYSICAL_M
       USE MCSETTING_M
       USE BUCKINGHAM_M
       USE FLUCTUATION_M 
       USE SUBBOX_M 
       implicit none
       
       integer i, j,k,m

       real*8  EnergyC
       real*8  dc2, distanceZ
       real*8  UC, xclj, yclj, zclj
       real*8  xcmcij,ycmcij,zcmcij
       real*8  sigmaSF,welldepthSF
       real*8  xmcn,ymcn,zmcn,mcd2,smcy
       real*8  xcmcn,ycmcn,zcmcn,cmcd2
       real*8  BU, gSM
       real*8  deltay,deltay1,cy1
       real*8  LJy1(20), sLJy(20)
       real*8  EnergyCi, EnergyCij   
        
!      --------------------------------
!      INITIALISE THE ENERGY 
!      --------------------------------
       EnergyC  = 0.0E0
!      ---------------------------------------------------
!      INTERACTION ENERGY BETWEEN PARTICLE AND CARBON ATOM
!      ---------------------------------------------------
          if(steele)then
             call steelepotential(i,UC)
             EnergyC = EnergyC + UC 
          endif
          if(CrowellChang)then
             call CrowellChangPotential(i,UC)
             EnergyC = EnergyC + UC 
          endif


          if(Bojan)then
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
          
     RETURN
     ENDSUBROUTINE energysingleparticleSF 
     
     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    SUBROUTINE INITIALNPART
!    GET THE INFORMATION FOR THE INPUTTED NPART IF NOT START WITH AN EMPTY BOX
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
     subroutine InitialNpart
     USE FLUID_SOLID_M
     USE MCSETTING_M
     USE PHYSICAL_M
     USE SUBBOX_M
     implicit none
     
     integer N, N1, N3
     integer i, j, k, iCount, iSelected 
          
     real*8 rdn
     
     N = Npart
     if(Npart .eq. 0)return     
     if(SDisplace)then
        NinSBox = 0
        SeqNinSBox= 0
        indNinSBox = 0
     endif
!    ------------------------------------------------------------
!    NOW RANDOMLY DISTRIBUTE N PARTICLES IN N3 POSSIBLE LOCATIONS
!    ------------------------------------------------------------
     DO i=1, N
        do k = 1, SubBox
           IF(HorizontalBin)THEN
               if(MCz(i).GE.LayerZL(k) .AND. MCz(i).LT.LayerZH(k))then
                   NinSBox(k) = NinSBox(k)+1
                   SeqNinSBox(k,int(NinSBox(k))) = i
                   indNinSBox(i)=k
                   EXIT
               endif
           ELSEIF(VerticalBin)THEN
               if(MCy(i).GE.LayerYL(k) .AND. MCy(i).LT.LayerYH(k))then
                   NinSBox(k) = NinSBox(k)+1
                   SeqNinSBox(k,int(NinSBox(k))) = i
                   indNinSBox(i)=k
                   EXIT
               endif
           ENDIF
        enddo
     ENDDO
     
     RETURN
     ENDSUBROUTINE InitialNpart
     
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
     implicit none
     
     integer i,j,k
     integer switch, NDD
     integer i1, i2, i3, iCount
     integer iBin, inX, inY, inZ, iSec 
     integer MaxBinfAve
     
     real*8 Mindelta
     
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
             BinVolSec(i) = CarbonLengthx1*deltaY2D(i)*deltaZ2D(i)
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
 !         open(unit=26, file='2DLocalD.txt',status='unknown')
 !         rewind(26)
     ENDIF
     
     IF(switch.eq.1)THEN
         sum_Nin2DBin = 0.0e0
         ave_Nin2DBin = 0.0e0
         NDD          = 0
     ENDIF
     
     IF(switch.eq.2)THEN
        NDD = NDD + 1
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
           sum_Nin2DBin(i) = sum_Nin2DBin(i) + Nin2DBin(i)
        ENDDO
     ENDIF

     IF(switch.eq.3)THEN    
          do i = 1, totBin2D
             ave_Nin2DBin(i) = sum_Nin2DBin(i)/dble(NDD)
         enddo 
!         write(26,'(I)')iPC
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
     
     integer i, j, fixUnit, bflag
     
     character*200 filename, filename1, filename2
     
         
     SumD_Bin2D  = 0.0e0
           
     DO i = 1, totBin2D
        DO j = 1, NBinfAve(i)
            SumD_Bin2D(i)  = SumD_Bin2D(i) + ave_Nin2DBin(SeqBinfAve(i,j))
       ENDDO
     ENDDO
     
     unitNo = unitNo + 1
     fixUnit = unitNo
     write(filename1,'(F12.2)') P(iPC)
     write(filename2,'(I3)') iPC
     filename = trim(filename1)//trim(filename2)//'Pa2D.txt'
     open(unit=fixUnit, file=filename,status='unknown')
     rewind(fixUnit)
     write(fixUnit,'(I4,I15)') iPC,totBin2D
     
     Do i = 1, totBin2D    
              ! Remove the bins that not inside of the pore 
              if(BojanSlit)then      
                 bflag = 0
                 do j = 1, NS
                     if( Bin2DY(i).gt.StripLy(j) .AND. Bin2DY(i).le.StripHy(j)  )then
                         if(   Bin2DZ(i) .lt. InLzLyStripZ(j) .and. &
             &                (InLzLyStripZ(j)-Bin2DZ(i))   .GE. (Bin2DY(i)-StripLy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT
                          elseif(   Bin2DZ(i) .gt. InHzLyStripZ(j) .and. &
             &                 (Bin2DZ(i)-InHzLyStripZ(j))   .GE. (Bin2DY(i)-StripLy(j))*tan(StripAngleY(j))   ) then 
                               bflag = 1
                               EXIT  
                         endif
                   endif
                enddo
                if(bflag .eq. 1)CYCLE
             endif                       
             
        Smooth_DDensity2D(i) = ( SumD_Bin2D(i) + ave_Nin2DBin(i) )/(SumV_Bin2D(i) + Bin2DVol(i))
        DDensity2DKmolperM3(i) = Smooth_DDensity2D(i)/6.0230E23/(Scalelength*1.0e-10)**3/1.0e3
!        if(DDensity2DKmolperM3(i).gt.0.0)then
!           IF(GCMC)write(26,'(4E16.6)')P(iPC), Bin2DY(i),Bin2DZ(i),DDensity2DKmolperM3(i)
!            IF(GCMC)write(fixUnit,'(4E16.6)')P(iPC), Bin2DY(i),Bin2DZ(i),DDensity2DKmolperM3(i)
             IF(GCMC)write(fixUnit,'(6E16.6)')P(iPC), Bin2DY(i),Bin2DZ(i),DDensity2DKmolperM3(i), ave_Nin2DBin(i), Bin2DVol(i)
!        endif
     ENDDO

     RETURN
     ENDSUBROUTINE LocalD2DSmooth
     
     
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      SUBROUTINE BojanAdjust_Angle
!      Adjust the paramter for Bojan Strips by considering the angles formed with y- axis    
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
     SUBROUTINE BojanAdjust_Angle
     USE FLUID_SOLID_M
     USE PHYSICAL_M
     USE MCSETTING_M
     implicit none
     
     integer i, j
     
     
 
           do i = 1, NS
              StripWCor(i)     = StripW(i)/cos(StripAngleY(i))
              IF(StripAngleY(i).GE.0.0)THEN
                  InLzLyStripZ(i)  = StripZ0(i)
                  InHzLyStripZ(i)  = BoxLengthZ - StripZ0(i)
                  OutLzLyStripZ(i) = StripZ0(i)- (StriplayerN(i)-1)*Stripgap(i)/cos(StripAngleY(i))
                  OutHzLyStripZ(i) = BoxLengthZ - StripZ0(i) + (StriplayerN(i)-1)*Stripgap(i)/cos(StripAngleY(i))
              ELSEIF(StripAngleY(i).LT.0.0)THEN
                  InLzLyStripZ(i)  = StripZ0(i) + StripW(i)*tan(-StripAngleY(i))
                  InHzLyStripZ(i)  = BoxLengthZ - StripZ0(i) - StripW(i)*tan(-StripAngleY(i))
                  OutLzLyStripZ(i) = StripZ0(i)- (StriplayerN(i)-1)*Stripgap(i)/cos(StripAngleY(i)) + StripW(i)*tan(-StripAngleY(i))
                  OutHzLyStripZ(i) = BoxLengthZ - StripZ0(i) + (StriplayerN(i)-1)*Stripgap(i)/cos(StripAngleY(i)) - StripW(i)*tan(-StripAngleY(i))
              ENDIF                  
              StripLy(i)       = StripC(i)- StripW(i)/2.0e0
              StripHy(i)       = StripC(i)+ StripW(i)/2.0e0
           enddo
         
      return
      ENDSUBROUTINE BojanAdjust_Angle

