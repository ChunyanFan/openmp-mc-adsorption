!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PROGRAM kMCADS1
      USE  Constant_M
      USE  FLUID_SOLID_M
      USE  PHYSICAL_M
      USE  MCSETTING_M
      USE  SUBBOX_M
      USE  POREFIGURE_M
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     DELARATION OF VARIABLES
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
!     --------
!     INTEGERS
!     --------
 integer rank, size, ierror, tag
 
      integer i, cutoffopt,iP !, iden
      integer iExc, Nexc
      integer NC1,NC2
      integer j,k, InterfBin
      integer Part1Bin,Part2Bin,Part3Bin
!     --------------
!     REAL VARIABLES
!     --------------
      real*8 CriticalT
      real*8 AvePoreW
      real*8 LayerH,InterH
      real*8 ExBinZ
      real*8 PerLayerH, InterfL
      double precision:: pretim, entim
      
!     ----------
!     CHARACTERS
!     ----------
      character*200 filename
      character(255) prefix
      character(255) outpath
 include 'mpif.h'
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     CONSTANT
! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

   pretim = MPI_Wtime()
if(rank < 10) then
   write (prefix, "(A5,I1)") "Input", rank
   call chdir(prefix)
else
   write (prefix, "(A5,I2)") "Input", rank
   call chdir(prefix)
endif
      pi             = 4.0E0*ATAN(1.0E0)
      AvogadroNumber = 6.0230E23
      kB             = 1.38066E-23
      Rg             = 8.314E0    ! Gas constant J/(K.mol)
      hv             = 6.626E-34  !Plank number JS/molecular
      permittivity   = 8.854e-12  !F/m or A2S4/Kg.m3
      unite          = 1.602e-19   !C
      GrapheneLayer  = 3.354E0
      eV             = 1.60217646E-19 !J
            
      call random_seed
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     INPUT FILES
! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
      open( unit=62, file="./std_output.txt", status='unknown')
      rewind(1)
      
      write(62,*) '----------------------------------------'
      write(62,*) 'ENTER FILENAME TO STORE OUTPUT DETAILS ='
      write(62,*) '----------------------------------------'
      open( unit=51, file="./output.txt", status='unknown')
      rewind(1)
      
      write(62,*) '------------------------------------------------'
      write(62,*) 'ENTER FILENAME TO STORE POSITION OF PARTICLES = '
      write(62,*) '------------------------------------------------'
      open( unit=52, file='./particleposition.txt', status='unknown')
      rewind(2)
      

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!      INPUTS: OPERATING CONDITIONS
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
!      ----------------------------    
!      INPUTS: ADSORBATE PARAMETERS
!      ----------------------------     
        open(unit=50, file = './MultiFluid.txt',status='unknown')
        rewind(50)
        read(50,*)MW
        read(50,'(f6.2)')CriticalT
        read(50,*)MCx0, MCy0, MCz0
        read(50,*)NLJ
        do j = 1, NLJ
           read(50,*)LJx0(j), LJy0(j),LJz0(j)
           read(50,*)sigmaFF(j),wellDepthFF(j)
        enddo
        read(50,*)NCL
        do k = 1, NCL
           read(50,*)CLx0(k), CLy0(k), CLz0(k)
           read(50,*)Charges(k)
           Charges(k) = Charges(k)*unite
        enddo
        read(50,*)
        read(50,'(L)')OrienD
        read(50,'(L)')SymVector
        read(50,*)fSite
        read(50,*)tSite
        if(OrienD)then
           open(unit=36, file='./OrienD.txt',status='unknown')
           rewind(36)
        endif
        IF(NLJ.GT.1)THEN
           read(50,*)
           read(50,*)thetar     !(K)
           read(50,*)thetav     !(K)
           read(50,*)ChemPD0    !(eV)
           ChemPD0 = ChemPD0*eV  !(J)
           ChemPDe = ChemPD0 + kB*thetav   !(J)
        ENDIF
        
      close(unit=50)   
!      -----------------------------------------------------------
!      SCALING PARAMETERS: USE THOSE OF THE FIRST ATOM ON MOLECULE
!      -----------------------------------------------------------
       SCALELENGTH  = sigmaFF(1)
       SCALEENERGY  = WellDepthFF(1)
      
!      ----------------------------------------
!      DIMENSIONLESS THE  PARAMTER OF ADSORBENT
!      ----------------------------------------   
       GrapheneLayer  = GrapheneLayer/SCALELENGTH
       MCX0           = MCX0/SCALELENGTH
       MCY0           = MCY0/SCALELENGTH
       MCZ0           = MCZ0/SCALELENGTH
       do i = 1, NLJ        
          SigmaFF(i)     = SigmaFF(i)/SCALELENGTH
          WellDepthFF(i) = WellDepthFF(i)/SCALEENERGY    
          LJX0(i)        = LJX0(i)/SCALELENGTH 
          LJY0(i)        = LJY0(i)/SCALELENGTH
          LJZ0(i)        = LJZ0(i)/SCALELENGTH
       enddo
       do j = 1,NCL
          CLX0(j)      =  CLX0(j)/SCALELENGTH
          CLY0(j)      =  CLY0(j)/SCALELENGTH
          CLZ0(j)      =  CLZ0(j)/SCALELENGTH 
       enddo  
      if(OrienD)VectorL = sqrt( (LJX0(fSite)-LJX0(tSite))**2.0+(LJY0(fSite)-LJY0(tSite))**2.0+(LJZ0(fSite)-LJZ0(tSite))**2.0  )

      open(unit=54, file='./operateconditions.txt')
!      --------------------------- 
!      TYPE OF SIMULATION ENSEMBLE
!      ---------------------------
       read(54,*)NVT
       read(54,*)NVTDes
       read(54,*)SURFACE
       read(54,*)PORE
       read(54,*)ExtraBox
       read(54,*)EQUILIBRIUM
       read(54,*)KINETIC
       read(54,*)kVLE
       read(54,*)NBOURBIN
       read(54,*)MassTSFER
       IF(MassTSFER)THEN
          open(unit=11, file='MassTrans.txt',status='unknown')
          rewind(11)
          write(11,'(3A16)')'NCycles',  'ave_rhoBin(i)' !, 'FLUX_bin(i)'
       ENDIF
       read(54,*)ChangeXY
       read(54,*)MeanFPath
       read(54,*)ENTROPY
       read(54,*)DistanceD
       read(54,*)LocalD2D
       read(54,*)LowUi
       read(54,*)HighUi
       read(54,*)Npart
       read(54,*)AddN
       read(54,*)NumOfrho
       Allocate(rho(NumOfrho), SEinput(NumOfrho),P(NumOfrho))
       read(54,*)SEinput(1) 
       SEinput(1) = SEinput(1)*(SCALELENGTH*1.0E-10)**2*AvogadroNumber
       read(54,*)Psaturate
       read(54,'(f6.2)')T
       read(54,*)steele
       read(54,*)Nsteele
       ALLOCATE(Psteele(Nsteele))
       do i = 1, Nsteele
         read(54,*)Psteele(i)
         Psteele(i) = Psteele(i)/SCALELENGTH
       enddo
       read(54,*)Bojan
       read(54,*)Close1End
       read(54,*)Close2Ends
            
       if(Close1End .OR. Close2Ends)then
            open(unit=112,file='CloseEnd.txt',status='unknown')
            rewind(112)
            read(112,*)Clo1Extend
            read(112,*)Clo1Hard
            read(112,*)Clo1LocZ
            read(112,*)Clo1LengthZ
            read(112,*)opClo1LengthZ
            read(112,*)Clo1LayerN
            read(112,*)Clo1gap
            read(112,*)Clo1SigmaSS
            read(112,*)Clo1WelldepthSS
            read(112,*)Clo1rhosperM2
            read(112,*)
            Clo1LocZ = Clo1LocZ/SCALELENGTH
            Clo1LengthZ = Clo1LengthZ/SCALELENGTH
            opClo1LengthZ = opClo1LengthZ/SCALELENGTH
            Clo1gap = Clo1gap/SCALELENGTH
            Clo1SigmaSS = Clo1SigmaSS/SCALELENGTH
            Clo1WelldepthSS = Clo1WelldepthSS/SCALEENERGY
            Clo1rhosperM2 = Clo1rhosperM2*((SCALELENGTH*1.0E-10)**2)    
            Clo1Mid_z = Clo1LocZ + Clo1LengthZ/2.0e0   
          endif
          
          if(Close2Ends)then
            read(112,*)Clo2LocZ
            read(112,*)Clo2LengthZ
            read(112,*)opClo2LengthZ
            read(112,*)Clo2LayerN
            read(112,*)Clo2gap
            read(112,*)Clo2SigmaSS
            read(112,*)Clo2WelldepthSS
            read(112,*)Clo2rhosperM2
            Clo2LocZ = Clo2LocZ/SCALELENGTH
            Clo2LengthZ = Clo2LengthZ/SCALELENGTH
            opClo2LengthZ = opClo2LengthZ/SCALELENGTH
            Clo2gap = Clo2gap/SCALELENGTH
            Clo2SigmaSS = Clo2SigmaSS/SCALELENGTH
            Clo2WelldepthSS = Clo2WelldepthSS/SCALEENERGY
            Clo2rhosperM2 = Clo2rhosperM2*((SCALELENGTH*1.0E-10)**2) 
            Clo2Mid_z = Clo2LocZ + Clo2LengthZ/2.0e0   
          endif
       
       IF(steele)THEN
          open(unit=58, file='./Solid.txt')
          read(58,*)sigmaSS
          read(58,*)WellDepthSS
          SigmaSS = SigmaSS/SCALELENGTH
          WellDepthSS = WellDepthSS/SCALEENERGY
       ENDIF
       
       IF(Bojan)THEN
            open(unit=59, file='./Bojan.txt')
            read(59,*)BojanSlit
            read(59,*)BoPBCy
            read(59,*)NS
            do i = 1, NS
              read(59,*)TopStrip(i)  ! If the strip is on the top of the simulation box or not
              read(59,*)StripZ0(i)
              read(59,*)StripC(i)
              read(59,*)StripW(i)
              read(59,*)StriplayerN(i)
              read(59,*)Stripgap(i)
              read(59,*)BsigmaSS(i)
              read(59,*)BwelldepthSS(i)
              read(59,*)rhosperM2(i)
              read(59,*)
              StripZ0(i)     = StripZ0(i)/SCALELENGTH                 ! The Z coordinate of the top layer of strip i
              StripC(i)      = StripC(i)/SCALELENGTH                  ! The coordinate of the center of the strip i
              StripW(i)      = StripW(i)/SCALELENGTH                  ! The width of strip i
              Stripgap(i)    = Stripgap(i)/SCALELENGTH                ! The gap spacing between two adjacent layers of strip i 
              BSigmaSS(i)    = BSigmaSS(i)/SCALELENGTH
              BWellDepthSS(i)= BWellDepthSS(i)/SCALEENERGY
              rhosperM2(i)   = rhosperM2(i)*((SCALELENGTH*1.0E-10)**2)
            enddo
        ENDIF

!      -----------   
!      TEMPERATURE
!      -----------
       CriticalTemp   = CriticalT/SCALEENERGY
       Temperature    = T/SCALEENERGY  
       deBro   = ((hv*sqrt(AvogadroNumber/(2*pi*MW*T*kB)))/(SCALELENGTH*1.0D-10)) ! thermal de broglie wave length
	   deBro3  = (hv*sqrt(AvogadroNumber/(2*pi*MW*T*kB)))**3/((SCALELENGTH*1D-10)**3) 
!      ----------------------------------
!      PREPARATION FOR DIFFERENT ENSEMBLE
!      ----------------------------------
       IF(SURFACE)THEN
             open(unit=53, file='./SurfBox.txt', status='unknown')
             rewind(53)
             read(53,*)BoxLengthX 
             read(53,*)BoxLengthY
             read(53,*)BoxLengthZ
             read(53,*)AccVZH
             read(53,*)MaxLayer
             read(53,*)LayerH
             read(53,*)InterLayer
             read(53,*)InterH
             read(53,*)InsertZL
             read(53,*)InsertZH
             IF(NSteele.eq.1)THEN
                SubBox  = MaxLayer + InterLayer + 1
                SubFGas = SubBox
             ENDIF
             IF(NSteele.eq.2)THEN
                SubBox  = MaxLayer*2 + InterLayer*2 + 1  !USE BIG SLIT PORE FOR SURFACE
                SubFGas = MaxLayer + InterLayer + 1
             ENDIF
             LayerH = LayerH/SCALELENGTH
             InterH = InterH/SCALELENGTH
             BoxLengthX = BoxLengthX/SCALELENGTH
             BoxLengthY = BoxLengthY/SCALELENGTH
             BoxLengthZ = BoxLengthZ/SCALELENGTH
             AccVZH = AccVZH/SCALELENGTH
             InsertZL = InsertZL/SCALELENGTH
             InsertZH = InsertZH/SCALELENGTH
             ALLOCATE(LayerZL(SubBox+1),LayerZH(SubBox+1), SubVol(SubBox+1))
             ALLOCATE(SeqNinSBox(SubBox+1,50000),totNinSBox(SubBox+1),ave_NinSBox(SubBox+1))
             ALLOCATE(rho_SBox(SubBox+1),rho_SBoxMolPerM3(SubBox+1))
             ALLOCATE(rateSBox(SubBox+1),rateTSBox(SubBox+1), timeSBox(SubBox+1), ChemPSBox(SubBox+1))
             LayerZL(1) = 0.0E0
             LayerZH(1) =  LayerZL(1)+LayerH
             do i = 2, MaxLayer
               LayerZL(i) =  LayerZH(i-1)
               LayerZH(i) =  LayerZL(i)+LayerH             
             enddo
             do i = MaxLayer+1, MaxLayer + InterLayer
                LayerZL(i) =  LayerZH(i-1)
                LayerZH(i) =  LayerZL(i)+InterH    
             enddo
             IF(NSteele.eq.1)THEN
                LayerZL(SubBox) = LayerZH(SubBox-1)
                LayerZH(SubBox) = BoxLengthZ
             ELSEIF(NSteele.eq.2)THEN
                LayerZH(SubBox) = BoxLengthZ
                LayerZL(SubBox) = LayerZH(SubBox)- LayerH
                do i=SubBox-1, SubBox-MaxLayer+1,-1
                   LayerZH(i) = LayerZL(i+1)
                   LayerZL(i) = LayerZH(i)- LayerH
                enddo
                do i=SubBox-MaxLayer, SubBox-MaxLayer-InterLayer+1,-1
                   LayerZH(i) = LayerZL(i+1)
                   LayerZL(i) = LayerZH(i)- InterH
                enddo
                   LayerZH(MaxLayer+InterLayer+1)=LayerZL(MaxLayer+InterLayer+2)
                   LayerZL(MaxLayer+InterLayer+1)=LayerZH(MaxLayer+InterLayer)
            ENDIF
             do i =1, SubBox
                SubVol(i) = BoxLengthX*BoxLengthY*(LayerZH(i)-LayerZL(i))
             enddo
             IF(MassTSFER)THEN
                NFluxSurf = SubBox - 1
                write(11,'(5X,500I5)') (i,i=1,SubBox), (i, i=1,NFluxSurf)
                ALLOCATE(FluxSurf(NFluxSurf),UpFlux(NFluxSurf), DownFlux(NFluxSurf),NetFlux(NFluxSurf))
                do i = 1,SubBox-1
                   FluxSurf(i) = LayerZH(i)
                enddo
             ENDIF            
       ENDIF

       IF(PORE)THEN
             open(unit=56, file='./PoreBox.txt',status='unknown')
             rewind(56)
             read(56,*)BoxLengthX 
             read(56,*)BoxLengthY
             read(56,*)BoxLengthZ
             read(56,*)AccVYL
             read(56,*)AccVYH
             read(56,*)NPartinY
             read(56,*)Part1Ly
             read(56,*)Part2Lz    ! pore width 
             read(56,*)Part1Bin, Part2Bin, Part3Bin
             BoxLengthX = BoxLengthX/SCALELENGTH
             BoxLengthY = BoxLengthY/SCALELENGTH
             BoxLengthZ = BoxLengthZ/SCALELENGTH
             AccVYL = AccVYL/SCALELENGTH
             AccVYH = AccVYH/SCALELENGTH
             Part1Ly = Part1Ly/SCALELENGTH
             Part2Lz = Part2Lz/SCALELENGTH
             SubBox  = Part1Bin+Part2Bin+Part3Bin
             SubFGas = 2  ! FOR CALCULATING THE PROPERTIES OF THE BULK PHASE

             if(ExtraBox)then
                SubBox  = Part1Bin+Part2Bin+Part3Bin + 1
                SubFGas = SubBox
             endif

             IF(NPartinY.GE.2)THEN
                 AccVZL = (BoxLengthZ - Part2Lz)/2.0E0  ! The LOW boundary of the pore section
                 AccVZH = (BoxLengthZ + Part2Lz)/2.0E0  ! The HIGH boundary of the pore section
             ELSE
                AccVZL = 0.0E0
                AccVZH = BoxLengthZ
             ENDIF
             IF(NPartinY.EQ.3)then
               Part3Ly = Part1Ly
               Part2Ly = BoxLengthY - Part1Ly - Part3Ly
             ELSE
               Part2Ly = BoxLengthY - Part1Ly
             ENDIF
             ALLOCATE(LayerXL(SubBox),LayerXH(SubBox),LayerYL(SubBox),LayerYH(SubBox) )
             ALLOCATE(LayerZL(SubBox),LayerZH(SubBox), SubVol(SubBox))
             ALLOCATE(SeqNinSBox(SubBox,10000),totNinSBox(SubBox),ave_NinSBox(SubBox))
             ALLOCATE(rho_SBox(SubBox+1),rho_SBoxMolPerM3(SubBox+1))
             ALLOCATE(rateSBox(SubBox+1),rateTSBox(SubBox+1), timeSBox(SubBox+1), ChemPSBox(SubBox+1))
             LayerXL = 0.0E0
             LayerXH = LayerXL(1)+BoxlengthX
  !           LayerYL(1) = 0.0E0
  !           LayerYH(1) = LayerYL(1) + Part1Ly/dble(Part1Bin)
             DO i = 1, Part1Bin
                LayerZL(i) = 0.0E0
                LayerZH(i) = LayerZL(i)+ BoxlengthZ
  !              if(i.GE.2)then
  !                 LayerYL(i) = LayerYH(i-1)
  !                 LayerYH(i) = LayerYL(i) + Part1Ly/dble(Part1Bin)
  !              endif
                 read(56,*)LayerYL(i), LayerYH(i)
                 LayerYL(i) = LayerYL(i)/SCALELENGTH
                 LayerYH(i) = LayerYH(i)/SCALELENGTH
             ENDDO
             
             DO i = Part1Bin+1, Part1Bin+Part2Bin
                LayerZL(i) = (BoxLengthZ-Part2Lz)/2.0e0
                LayerZH(i) =  LayerZL(i) + Part2Lz
                LayerYL(i) = LayerYH(i-1)
                LayerYH(i) = LayerYL(i) + Part2Ly/dble(Part2Bin)
             ENDDO
             
             DO i = Part1Bin+Part2Bin+1, Part1Bin+Part2Bin+Part3Bin
                LayerZL(i) = 0.0E0
                LayerZH(i) = LayerZL(i)+ BoxlengthZ
                LayerYL(i) = LayerYH(i-1)
                LayerYH(i) = LayerYL(i) + Part3Ly/dble(Part3Bin)
             ENDDO

             if(ExtraBox)then
                read(56,*)ExtraBoxH  ! The height of the ExtraBox
                ExtraBoxH = ExtraBoxH/SCALELENGTH
                LayerZL(SubBox) = BoxLengthZ
                LayerZH(SubBox) = BoxLengthZ + ExtraBoxH
                LayerYL(SubBox) = 0.0E0
                LayerYH(SubBox) = BoxLengthY
             endif

             do i=1, SubBox
                SubVol(i) = (LayerXH(i)-LayerXL(i))*(LayerYH(i)-LayerYL(i))*(LayerZH(i)-LayerZL(i))
             enddo
            IF(MassTSFER)THEN
                NFluxSurf = SubBox - 1
                write(11,'(5X,50I5)') (i,i=1,SubBox), (i, i=1,NFluxSurf)
                ALLOCATE(FluxSurf(NFluxSurf),UpFlux(NFluxSurf), DownFlux(NFluxSurf),NetFlux(NFluxSurf))
                do i = 1,SubBox-1
                   FluxSurf(i) = LayerYH(i)
                enddo
             ENDIF            
       ENDIF

       IF(kVLE)THEN
            open(unit=29,file='./IniVLE.txt',status='unknown')
            rewind(29)
            read(29,*)HomoStart
            read(29,*)ChangeXZ
            read(29,*)VLEradius
            read(29,*)VLEdistance
            read(29,*)WIDOM_VLE
            read(29,*)NGhost
            read(29,*)SolidSlab
            read(29,*)NunitX,NunitZ,NunitY ! The number of unit cells in X-,Y-,Z- directions
            read(29,*)LiqDens   !!! mol/m3
            read(29,*)VLEBoxLengthX
            read(29,*)VLEBoxLengthY
            read(29,*)VLEBoxLengthZ
            read(29,*)LiqLy
            read(29,*)LiqHy
            read(29,*)InterfL
            read(29,*)LiqBin
            read(29,*)GasBin
            read(29,*)InterfBin
            LiqDens = LiqDens*AvogadroNumber*((SCALELENGTH*1.0E-10)**3)
            VLEBoxLengthX = VLEBoxLengthX/SCALELENGTH
            VLEBoxLengthY = VLEBoxLengthY/SCALELENGTH
            VLEBoxLengthZ = VLEBoxLengthZ/SCALELENGTH
            VLEBoxLengthX0 = VLEBoxLengthX
            VLEBoxLengthY0 = VLEBoxLengthY
            VLEBoxLengthZ0 = VLEBoxLengthZ
            LiqLy = LiqLy/SCALELENGTH
            LiqHy = LiqHy/SCALELENGTH 
            PropLy(1) = 0.0e0
            PropHy(1) = 4.0e0
            PropLy(2) = VLEBoxLengthY/2.0e0 - 2.0e0
            PropHy(2) = VLEBoxLengthY/2.0e0 + 2.0e0
            PropLy(3) = VLEBoxLengthY - 4.0e0
            PropHy(3) = VLEBoxLengthY 
            LunitX = VLEBoxLengthX/dble(NunitX)
            LunitZ = VLEBoxLengthZ/dble(NunitZ)
            LunitY = (LiqHy-LiqLy)/dble(NunitY) 
       ENDIF
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!      PARAMETERS FOR NUMBER OF CYCLES, boundary condition
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       read(54,*)pbcx
       read(54,*)pbcy
       read(54,*)pbcz
       read(54,*)Ncycle
       read(54,*)Ncycles
       read(54,*)Dist_Limit
       read(54,*)PairEff_Limit
       read(54,*)PairEsf_Limit
       CALL LimitPotentialEnergy(PairV_Limit,DerivE_Limit)
       Nmove = 1000
                     
       IF(kVLE)THEN
             SubBox = LiqBin + 2*GasBin + 2*InterfBin
             ALLOCATE(LayerXL(SubBox+1),LayerXH(SubBox+1),LayerYL(SubBox+1),LayerYH(SubBox+1) )
             ALLOCATE(LayerZL(SubBox+1),LayerZH(SubBox+1), SubVol(SubBox+1))
             ALLOCATE(IniStepX(SubBox+1),IniStepY(SubBox+1),IniStepZ(SubBox+1))
             ALLOCATE(SStepX(SubBox+1),SStepY(SubBox+1),SStepZ(SubBox+1))
             ALLOCATE(SeqNinSBox(SubBox+1,10000),totNinSBox(SubBox+1),ave_NinSBox(SubBox+1))
             ALLOCATE(AtpinMove(SubBox+1), SucinMove(SubBox+1),SucMovein(SubBox+1),SucMoveout(SubBox+1))
             ALLOCATE(rho_SBox(SubBox+1),rho_SBoxMolPerM3(SubBox+1))
             ALLOCATE(rhoFW(SubBox+1),ChemicalPVLE(SubBox+1), ExcessChemicalPVLE(SubBox+1), IdealChemicalPVLE(SubBox+1))         
             ALLOCATE(ChemPVLETotal(SubBox+1), ExChemPVLETotal(SubBox+1), IdChemPVLETotal(SubBox+1))
             ALLOCATE(ChemPVLEAve(SubBox+1), ExChemPVLEAve(SubBox+1), IdChemPVLEAve(SubBox+1))
             ALLOCATE(ChemPVLEJPerMole(SubBox+1), ExChemPVLEJPerMole(SubBox+1), IdChemPVLEJPerMole(SubBox+1))     

             LayerXL = 0.0E0
             LayerXH = LayerXL(1)+VLEBoxlengthX
             LayerZL = 0.0E0
             LayerZH = LayerZL(1)+VLEBoxlengthZ
             LayerYL(1) = 0.0E0
             perLayerH = (LiqLy-InterfL/2.0)/dble(GasBin)
             if(perLayerH.le.0)then
               write(62,*)'THE INTERFICAL REGION IS TOO LARGE'
               write(51,*)'THE INTERFICAL REGION IS TOO LARGE'
               stop
             endif
             LayerYH(1) = LayerYL(1)+ perLayerH  
             do i=2, GasBin           
                LayerYL(i) = LayerYH(i-1)
                LayerYH(i) = LayerYL(i) + perLayerH  
             enddo
             
             perLayerH = InterfL/dble(InterfBin)
             if(perLayerH.le.0)then
               write(62,*)'THE INTERFICAL REGION IS TOO LARGE'
               write(51,*)'THE INTERFICAL REGION IS TOO LARGE'
               stop
             endif
             do i=GasBin+1, GasBin+InterfBin
                LayerYL(i) = LayerYH(i-1)
                LayerYH(i) = LayerYL(i)+ perLayerH
             enddo
              
             perLayerH = (LiqHy-LiqLy-InterfL)/dble(LiqBin)     
             if(perLayerH.le.0)then
               write(62,*)'THE INTERFICAL REGION IS TOO LARGE'
               write(51,*)'THE INTERFICAL REGION IS TOO LARGE'
               stop
             endif        
             do i =GasBin+InterfBin+1, GasBin+InterfBin+LiqBin
                LayerYL(i) = LayerYH(i-1)
                LayerYH(i) = LayerYL(i)+perLayerH
             enddo
             
             perLayerH = InterfL/dble(InterfBin)
             do i =GasBin+InterfBin+LiqBin+1, GasBin+2*InterfBin+LiqBin
                LayerYL(i) = LayerYH(i-1)
                LayerYH(i) = LayerYL(i)+ perLayerH
             enddo
             
             perLayerH = (VLEBoxlengthY-LiqHy-InterfL/2.0)/dble(GasBin)
             do i =GasBin+2*InterfBin+LiqBin+1, SubBox
                LayerYL(i) = LayerYH(i-1)
                LayerYH(i) = LayerYL(i)+ perLayerH
             enddo
             do i=1, SubBox
                SubVol(i) = (LayerXH(i)-LayerXL(i))*(LayerYH(i)-LayerYL(i))*(LayerZH(i)-LayerZL(i))
                IniStepX(i) = (LayerXH(i)-LayerXL(i))/2.0E0
                IniStepY(i) = (LayerYH(i)-LayerYL(i))/2.0E0
                IniStepZ(i) = (LayerZH(i)-LayerZL(i))/2.0e0
             enddo
       ENDIF

      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!      PRINT THE DIMENSIONAL INPUTS
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
     
      write(62,1119) T,  Ncycle,  Ncycles, Nmove
      write(51,1119) T,  Ncycle,  Ncycles, Nmove
1119        FORMAT(1x, 'INPUT',/, &
    &        1x, '-----',/, &
    &        1x, 'Temperature (K)---------------------------------: ', f20.10,/, &
    &        1x, 'Number of cycles in equilibrium ----------------: ', i20,/, &
    &        1x, 'Number of cycles in sampling -------------------: ', i20,/, & 
    &        1x, 'Number of cycles of move/exchange particles ----: ', i20,/, //)
    
       write(62,1120) NLJ
1120        FORMAT(1x,'The number of LJ site----------------------: ', i10)
       do i = 1, NLJ
       write(62,1121) i,sigmaFF(i),welldepthFF(i)
1121        FORMAT(1X,'NLJ=',I5,5x,'sigma = ',f15.7, 5x,'welldepth = ', f15.7)        
       enddo
       
       write(62,1122) NCL
1122        FORMAT(1x,'The number of CL site---------------------: ', i10)
       do i = 1, NCL
       write(62,1123)i, Charges(i)
1123        FORMAT(1X,'NCL=',I5, 5x,'Charges = ', e15.7)        
       enddo
       

       write(62,1124) sigmaSS, WellDepthSS
       write(51,1124) sigmaSS, WellDepthSS
1124         FORMAT(1x, 'ADSORBENT PARAMETERS ',/, &
    &        1x, '---------------------',/, &  
    &        1x, 'sigmaSS (A)-----------------------------: ', f20.10,/, &
    &        1x, 'Well-DepthSS (K)------------------------: ', f20.10,/,//)
   
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    CALCULATION STARTS HERE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
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
     
         IF(MeanFPath.AND.SURFACE)THEN
           call StopPosition
         ENDIF
         
         if(kVLE)  call kMC_VLE
         if(SURFACE) call MC_SURFACE
         if(PORE.AND.EQUILIBRIUM)call MC_PORE
 !        if(PORE.AND.KINETIC)call kMC_PORE
 
         DEALLOCATE(rho, SEinput,P)
         DEALLOCATE(Psteele)
         IF(SURFACE)THEN
           DEALLOCATE(LayerZL,LayerZH, SubVol)
           DEALLOCATE(SeqNinSBox,totNinSBox,ave_NinSBox)
           DEALLOCATE(rho_SBox,rho_SBoxMolPerM3)
           DEALLOCATE(rateSBox,rateTSBox, timeSBox, ChemPSBox)
         ENDIF
         IF(PORE)THEN
           DEALLOCATE(LayerXL,LayerXH,LayerYL,LayerYH )
           DEALLOCATE(LayerZL,LayerZH, SubVol)
           DEALLOCATE(SeqNinSBox,totNinSBox,ave_NinSBox)
           DEALLOCATE(rho_SBox,rho_SBoxMolPerM3)
           DEALLOCATE(rateSBox,rateTSBox, timeSBox, ChemPSBox)
         ENDIF
         
         
      write(62,*) "Exits normally"
      entim = MPI_Wtime()
      write(62,*) "Time taken by Rank", rank, ":", entim-pretim
      call MPI_FINALIZE(ierror)
         
    END PROGRAM kMCADS1
