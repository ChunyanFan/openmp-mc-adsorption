!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     PROGRAM:   PoreProcess6
!     BASED ON:  PoreProcess1.f90 
!                
!     MODIFICATION: 
!     1. Handle pore with angles (conventional GCMC), positive angle
! 
!     Modification: 27-Feb-2013: remove the 2DLocalD bins that not inside of the pore, Subroutine "Local2DSmooth"
!     08-Oct-2013: 1. Add gap between the pore walls and the closed end 
!                  2. Handle pore with angles, negative angle    
!     DATE FOR THIS PROGRAM: 11-Dec-2012
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Summary of 'Unit' in this program
!          unit=51,   file='output.txt' !renaming 1 as 51
!          unit=62,   file='std_output.txt' !to write standard output written by each process to a seperate file
!          unit=52,   file='particleposition' 
!          unit=53,   file='solidposition.txt'
!          unit=54,   file='operateconditions.txt'
!          unit=50,   file='MultiFluid.txt'
!          unit=60,   file='Solid.txt'
!          unit=57,   file='isothermdata.txt'  
!          unit=58,   file='Bojan.txt'
!          unit=59,   file='PressureGCMC.txt'
!          unit=61,  file='firstsurface.txt'
!          unit=11,  file='firstsurface.txt'
!          unit=12,  file='secondsurface.txt'
!          unit=13,  file='CPfile.txt'
!          unit=14,  file='FCHARGE.txt'
!          unit=15,  file='NVTisotherm.txt'
!          unit=16,  file='StepLength.txt'
!          unit=17,  file='IniConfig.txt'
!          unit=18,  file='ForBuckingham.txt'
!          unit=19,  file='19ForNVT.txt'
!          unit=19,  file='allinsert'
!          unit=20,  file='reasonableinsert'
!          unit=21,  file='accparticle'
!          unit=22,  file='SubBox.txt'
!          unit=23,  file='Move-Cycle.txt'
!          unit=24,  file='IniBin.txt'
!          unit=25,  file='Ini2DLocalD.txt'
!          unit=26,  file='2DLocaD.txt' 
!
!          unit=29,  file='radius.txt'
!
!          unit=31,  file='distanceD'
!          unit=31,  file='firstsurface.txt'
!          unit=32,  file='secondsurface.txt'
!          unit=32,  file='Fluctuation.txt'
!          unit=33,  file='Flu_LocalENex.txt'
!          unit=34,  file='LocalISOH.txt'
!          unit=35,  file='LayerISOH.txt'
!          unit=36,  file='OrienD.txt'
!          unit=41,  file='Sep_distanceD'
!          unit=42,  file='CCompressResult.txt' 
! 
!          unit=111, file='InclineBojan.txt'
!          unit=112, file='CloseEnd.txt'
!          unit=113, file='SeparateLD.txt'
!
!          unit=200+, file='2D.txt'   !!For output the results of 2D Density distribution for every point
!             
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PROGRAM PoreProcess6
!   include 'mpif.h'
      USE  Constant_M
      USE  FLUID_SOLID_M
      USE  PHYSICAL_M
      USE  MCSETTING_M
      USE  BUCKINGHAM_M
      USE  FLUCTUATION_M
      USE  POREFIGURE_M
      USE  FLEXICLINFO_M
      USE  SUBBOX_M
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
      integer modeFinite, modePore,secdsurf, modesteele
      integer j,k, iSub
!     --------------
!     REAL VARIABLES
!     --------------
      real*8 rcutoffset
      real*8 CarbonLengthx2, CarbonLengthy2
 !     real*8 zv
      real*8 SurfaceArea1,SurfaceArea2
      real*8 CriticalT
      real*8 AvePoreW
      real*8 LayerH,InterH
      real*8 ExBinZ
      double precision:: pretim, entim

!     ----------
!     LOGICAL
!     ----------
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
!   call chdir(prefix)
      pi             = 4.0E0*ATAN(1.0E0)
      AvogadroNumber = 6.0230E23
      kB             = 1.38066E-23
      Rg             = 8.314E0    ! Gas constant J/(K.mol)
      hv             = 6.626E-34  !Plank number JS/molecular
      permittivity   = 8.854e-12  !F/m or A2S4/Kg.m3
      unite          = 1.602e-19   !C
      GrapheneLayer  = 3.354E0
      
      ENERGYMATRIX   = 0.0E0
      
      call random_seed
          
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     INPUT FILES
! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
      open( unit=62, file="./std_output.txt", status='unknown')
      rewind(1)
      
      open( unit=51, file="./output.txt", status='unknown')
      rewind(1)

      write(62,*) '----------------------------------------'
      write(62,*) 'ENTER FILENAME TO STORE OUTPUT DETAILS ='
      write(62,*) '----------------------------------------'
!      read(62,*) filename
      
      write(62,*) '------------------------------------------------'
      write(62,*) 'ENTER FILENAME TO STORE POSITION OF PARTICLES = '
      write(62,*) '------------------------------------------------'
!      read(62,*) filename
      open( unit=52, file='./particleposition.txt', status='unknown')
      rewind(2)
      
      open( unit=53, file='./solidposition.txt', status='unknown')
      rewind(3)
      
      
      open( unit=57, file='./isothermdata.txt', status='unknown')
      rewind(7)
      write(57,*)'-----------------------------------------------------------------------------------------------------------------------------------------------------------------'
      write(57,'(A4,18A14)') 'no.','rCutoff','P','P/p0','SurfE','Nb_rho','Ex_Nb_rho','Iso_H','Iso_HFF','Iso_HSF','nIso_H','nIso_HFF','nIso_HSF','nnIso_H','Av_Npart','Av_Ninpore','Ad_Ninpore','CoorN-1','CoorN-2' 
      write(57,'(5x,6A14)')'(-)','(pa)','(-)','(umol/m2)','(mol/m3)','(mol/m3)','(kJ/mol)'
      write(57,*)'-----------------------------------------------------------------------------------------------------------------------------------------------------------------'

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
        read(50,*)MW   ! MOLECULAR WEIGHT
        read(50,'(f6.2)')CriticalT
        read(50,*)MCx0, MCy0, MCz0
        read(50,*)NLJ
        do j = 1, NLJ
           read(50,*)LJx0(j), LJy0(j),LJz0(j)
           read(50,*)sigmaFF(j),wellDepthFF(j)
        enddo
        read(50,'(I)')NCL
        do k = 1, NCL
           read(50,*)CLx0(k), CLy0(k), CLz0(k)
           read(50,*)Charges(k)
           Charges(k) = Charges(k)*unite
        enddo
        read(50,*)
        read(50,'(L)')OrienD
        if(OrienD)then
           open(unit=36, file='./OrienD.txt',status='unknown')
           rewind(36)
           read(50,'(L)')SymVector
           read(50,'(I)')fSite
           read(50,'(I)')tSite
        endif
   close(unit=50)

       open(unit=54, file='./operateconditions.txt')
!      --------------------------- 
!      TYPE OF SIMULATION ENSEMBLE
!      ---------------------------
       read(54,*)GCMC
       read(54,*)INPUTCP
       read(54,*)BULKGCMC
       read(54,*)LOCALF
       IF(LOCALF)BULKGCMC=.TRUE.
       read(54,*)NVT
       read(54,*)BUCKINGHAM
       read(54,*)JOHNSON
       read(54,*)WAGNER
       read(54,'(L)')SDisplace 
       read(54,'(L)')HorizontalBin
       read(54,'(L)')VerticalBin
       read(54,'(L)')TradMC
       read(54,'(L)')ExtraBin
       read(54,'(L)')ChangeXY
       read(54,'(L)')DistanceD
       read(54,'(L)')RadiusD
       read(54,'(L)')LocalD2D
       read(54,'(L)')DecayESF
       if(BUCKINGHAM)then
          open(unit=18, file='./ForBuckingham.txt',status='unknown')
          rewind(18)
          read(18,*)alpha_B
          read(18,*)WelldepthFF(1)
          read(18,*)SigmaFF(1)
          read(18,*)AFF_B
          read(18,*)BFF_B
          read(18,*)alphaFF_B
          read(18,*)WelldepthSS
          read(18,*)ETA
          read(18,*)ASS_B
          read(18,*)BSS_B
          read(18,*)alphaSS_B          
       endif
       
!     -----------------------------------------------------------
!     SCALING PARAMETERS: USE THOSE OF THE FIRST ATOM ON MOLECULE
!     -----------------------------------------------------------
      SCALELENGTH  = sigmaFF(1)
      SCALEENERGY  = WellDepthFF(1)
      
!     ----------------------------------------
!     DIMENSIONLESS THE  PARAMTER OF ADSORBENT
!     ----------------------------------------   
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
      if(BUCKINGHAM)then
          AFF_B = AFF_B/SCALEENERGY/SCALELENGTH**6.0E0
          BFF_B = BFF_B/SCALEENERGY 
          alphaFF_B = alphaFF_B*SCALELENGTH
          ETA = ETA*SCALELENGTH**2.0E0
          ASS_B = ASS_B/SCALEENERGY/SCALELENGTH**6.0E0
          BSS_B = BSS_B/SCALEENERGY 
          alphaSS_B = alphaSS_B*SCALELENGTH
          
          ASF_B = sqrt(AFF_B*ASS_B)
          BSF_B = sqrt(BFF_B*BSS_B)
          alphaSF_B = (alphaFF_B+alphaSS_B)/2.0e0
          
         call BuckinghamExp6
      endif
      
!     ----------------------------------------
!     THE BOX SIZE FOR THE BULKPHASE BOX-GCMC 
!     ----------------------------------------    
       open(unit=61, file='./firstsurface.txt')
       read(61,*)CarbonLengthx1         !! IF GCMC, this is for the bulk phase
       read(61,*)CarbonLengthy1
       read(61,*)BoxLengthZ
       
      
       CarbonLengthx1  =  CarbonLengthx1/SCALELENGTH
       CarbonLengthy1  =  CarbonLengthy1/SCALELENGTH
       BoxLengthZ      =  BoxLengthZ/SCALELENGTH 
       
       IF(NVT)THEN
         read(61,*)AccVYL
         read(61,*)AccVYH
         AccVYL          =  AccVYL/SCALELENGTH 
         AccVYH          =  AccVYH/SCALELENGTH 
         IniCarbonLengthx1 = CarbonLengthx1
         IniCarbonLengthy1 = CarbonLengthy1 
         IniBoxLengthZ     = BoxLengthZ
       ENDIF

       read(54,*)Psaturate
       read(54,'(f7.2)')T
       CriticalTemp   = CriticalT/SCALEENERGY
       Temperature    = T/SCALEENERGY  
       ! ---------------------------------------------
       ! Computation of thermal de Broglie wavelength
       ! ---------------------------------------------
       deBro   = ((hv*sqrt(AvogadroNumber/(2*pi*MW*T*kB)))/(SCALELENGTH*1.0D-10)) ! thermal de broglie wave length
	   deBro3  = (hv*sqrt(AvogadroNumber/(2*pi*MW*T*kB)))**3/((SCALELENGTH*1D-10)**3)

!      ----------------------------------
!      PREPARATION FOR DIFFERENT ENSEMBLE
!      ----------------------------------
       IF(GCMC)then
!       ---------------------------------------------------------------
!       SET PRESSURE POINTS AND CALCULATE THE RHO AND ZACTIVITY USE EOS
!       ---------------------------------------------------------------
        IF(INPUTCP)THEN
             open(unit=13,file='./CPfile.txt',status='unknown')
             rewind(13)
             read(13,*)NumOfPre
         ELSE
              open(unit=59,file='./PressureGCMC.txt',status='unknown')
              rewind(59)
              read(59,*)NumOfPre
         ENDIF
         Allocate(rho(NumOfPre), P(NumOfPre), Pressure(NumOfPre), zactivity(NumOfPre), InrefChemP(NumOfPre))
         Allocate(idealChemicalPotential(NumOfPre), ChemicalPotential(NumOfPre))
         Allocate(ResidualChemicalPotential(NumOfPre), ChemicalPotentialJoulepermol(NumOfPre))
         Allocate(Fun_NN_Bulk(NumOfPre), N_Bulk(NumOfPre), Fun_EN_Bulk(NumOfPre),CFun_NN_Bulk(NumOfPre)) 
         IF(INPUTCP)THEN
              read(13,*)
              do iP = 1 , NumOfPre
                  read(13,*)P(iP), Pressure(iP), rho(iP), InrefChemP(iP)
                  zactivity(iP) = exp(InrefChemP(iP)/Temperature)/deBro3
              enddo
         ELSE
              do iP = 1 , NumOfPre
                 read(59,*)P(iP)
                   Pressure(iP) = P(iP)*((SCALELENGTH*1.0E-10)**3)/(SCALEENERGY*kB)
              enddo
              ! ---------------------------------
              ! CALCULATE DENSITY BY JOHNSON EOS
              ! ---------------------------------
              call EquationOfState 
              call zactivitycalculation
              ! ---------------------------------
              ! CALCULATE THE CHEMICAL POTENTIAL
              ! ---------------------------------
              do iP = 1,NumOfPre
                 idealChemicalPotential(iP) = Temperature*log(deBro3*rho(iP))
                 ChemicalPotential(iP)      = ResidualChemicalPotential(iP) + idealChemicalPotential(iP)
                 ChemicalPotentialJoulepermol(iP) = ChemicalPotential(iP)*SCALEENERGY*Rg
                 write(51,*)'iP,ChemicalPotential,ChemicalPotentialJoulepermol',iP,ChemicalPotential(iP),ChemicalPotentialJoulepermol(iP)
!                ZV = (ChemicalPotential/Temperature) + log(deBro)
              enddo
         ENDIF
         CrowellChang = .FALSE.
         ADSORPTION = .FALSE.
         IF(BULKGCMC)call BULK_GCMC 
         ADSORPTION = .TRUE.
         if(BUCKINGHAM)CrowellChang = .TRUE.
         read(61,*)CarbonLengthx1  !! Simulation box for adsorption in GCMC
         read(61,*)CarbonLengthy1
         read(61,*)BoxLengthZ
         read(61,*)AccVYL
         read(61,*)AccVYH

         CarbonLengthx1  =  CarbonLengthx1/SCALELENGTH
         CarbonLengthy1  =  CarbonLengthy1/SCALELENGTH
         BoxLengthZ      =  BoxLengthZ/SCALELENGTH 
         AccVYL          =  AccVYL/SCALELENGTH 
         AccVYH          =  AccVYH/SCALELENGTH 
         IniCarbonLengthx1 = CarbonLengthx1
         IniCarbonLengthy1 = CarbonLengthy1
         IniBoxLengthZ     = BoxLengthZ
       ENDIF
       
       IF(NVT)then
          open(unit=19,file='./19ForNVT.txt',status='unknown')
          rewind(19)
          read(19,'(L)')NaddNVT
          read(19,'(L)')NVTDes
          read(19,*)addN
          read(19,*)AccVZH
          read(19,*)NumOfrho
          read(19,*)
          AccVZH   = AccVZH/SCALELENGTH
          Allocate(rho(NumOfrho), SEinput(NumOfrho),P(NumOfrho))
          IF(NaddNVT)THEN
            read(19,*)rho(1),SEinput(1)
            rho(1)= rho(1)*(SCALELENGTH*1.0E-10)**3*AvogadroNumber
            SEinput(1) = SEinput(1)*(SCALELENGTH*1.0E-10)**2*AvogadroNumber
          ELSE
            do iden = 1, NumOfrho
                read(19,*)rho(iden),SEinput(iden)   ! mol/m3, mol/m2
                rho(iden)= rho(iden)*(SCALELENGTH*1.0E-10)**3*AvogadroNumber     ! (-)
                SEinput(iden) = SEinput(iden)*(SCALELENGTH*1.0E-10)**2*AvogadroNumber !(-)
            enddo
          ENDIF
       ENDIF
!      ------------------------------ 
!      CONFIGUATION OF SIMULATION BOX
!      ------------------------------
       read(54,*)modePore
      
         InclinePore  =  .FALSE.
         SurfPore     =  .FALSE.
         
       select case (modePore)
       case(1)
         InclinePore  =  .TRUE.
       case(2)
         SurfPore     =  .TRUE.
       endselect
  
       read(54,*)steele
       IF(BUCKINGHAM)steele = .FALSE.
       read(54,*)Nsteele
       ALLOCATE(Psteele(Nsteele))
       do i = 1, Nsteele
         read(54,*)Psteele(i)
       enddo
      
       read(54,'(L)')Bojan
       read(54,'(L)')Close1End
       read(54,'(L)')Close2Ends
 !      IF(Close2Ends)Close1End = .FALSE.
       
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
          endif
       
       
          if(Bojan)then
            open(unit=58, file='./Bojan.txt')
            read(58,*)BojanSlit
            read(58,*)BoPBCy
            read(58,*)NS
            do i = 1, NS
              read(58,*)TopStrip(i)  ! If the strip is on the top of the simulation box or not
              read(58,*)StripZ0(i)
              read(58,*)StripC(i)
              read(58,*)StripW(i)
              read(58,*)StriplayerN(i)
              read(58,*)Stripgap(i)
              read(58,*)BsigmaSS(i)
              read(58,*)BwelldepthSS(i)
              read(58,*)rhosperM2(i)
              read(58,*)StripAngleY(i)
              read(58,*)
            enddo
         endif
         if(Bojan)then
           do i = 1, NS
            StripZ0(i)     = StripZ0(i)/SCALELENGTH                 ! The Z coordinate of the top layer of strip i
            StripC(i)      = StripC(i)/SCALELENGTH                  ! The coordinate of the center of the strip i
            StripW(i)      = StripW(i)/SCALELENGTH                  ! The width of strip i
            Stripgap(i)    = Stripgap(i)/SCALELENGTH                ! The gap spacing between two adjacent layers of strip i 
            BSigmaSS(i)    = BSigmaSS(i)/SCALELENGTH
            BWellDepthSS(i)= BWellDepthSS(i)/SCALEENERGY
            rhosperM2(i)   = rhosperM2(i)*((SCALELENGTH*1.0E-10)**2)
            StripAngleY(i) = StripAngleY(i)/180.0e0*pi                ! Convert the angle formed with Y-axis from 'degree' to 'radian'
            tan_StripAngleY(i) = tan(StripAngleY(i))
            sin_StripAngleY(i) = sin(StripAngleY(i))
            cos_StripAngleY(i) = cos(StripAngleY(i))
           enddo
           CALL BojanAdjust_Angle
         endif
     
       read(54,'(L)')Sep_rho_Dis
       read(54,'(L)')CorCompress
       if(Sep_rho_Dis)then
         open(unit=113,file='SeparateLD.txt',status='unknown')
         rewind(113)
         read(113,*)SepLDN      !Number of sections need to calculate the Density distribution
         do i = 1, SepLDN
            read(113,*)y_rhoBnd(i) !The boundary of y-coordiante for each section
            read(113,*)Sep_Lengthy(i) ! The length of the section in y-direction
            read(113,*)yL_ComBnd(i) !The Low  boundary of the y-coordiante for calculting the Local Compressibility
            read(113,*)yH_ComBnd(i) !The High boundary of the y-coordiante for calculting the Local Compressibility
            read(113,*)zL_ComBnd(i) !The Low  boundary of the z-coordiante for calculting the Local Compressibility
            read(113,*)zH_ComBnd(i) !The High boundary of the z-coordiante for calculting the Local Compressibility
            y_rhoBnd(i)=y_rhoBnd(i)/SCALELENGTH
            Sep_Lengthy(i) = Sep_Lengthy(i)/SCALELENGTH
            yL_ComBnd(i) = yL_ComBnd(i)/SCALELENGTH
            yH_ComBnd(i) = yH_ComBnd(i)/SCALELENGTH
            zL_ComBnd(i) = zL_ComBnd(i)/SCALELENGTH
            zH_ComBnd(i) = zH_ComBnd(i)/SCALELENGTH
         enddo
       endif
       
       read(54,'(L)')SMediation1
       read(54,*)ksi
       read(54,'(L)')SMediation2
       read(54,*)SMLimitZ
       read(54,*)SMLimitZ2
       read(54,*)ksiF
       SMLimitZ = SMLimitZ/SCALELENGTH
       SMLimitZ2 = SMLimitZ2/SCALELENGTH
       read(54,'(L)')FlexiCL     !If the location of charges is flexiable or not
       if(FlexiCL)then
         open(unit=14, file='./FCHARGE.txt', status='unknown')
         rewind(14)
         read(14,*)MinFFLJ   !The minimum Fluid-Fluid interaction between two adsorbate molecules (KJ/MOL)
         read(14,*)PNTfraction !The fraction of PENALTY in the MinFFLJ (-)
         read(14,*)Lequi0      ! the equilibrium bond length
         read(14,*)CLStepRange  ! The total range of the movement of charge (A) = The collosion dimater of the LJ site where the charge stay
         read(14,*)RdmCL        ! The fraction of Changing the location of charges in four different action in ONE SIMULATION CYCLE (-)
         read(14,*)NpickedCL    ! The number of CL sites that change the charges
         do i = 1, NpickedCL
            read(14,*)indexPickedCL(i),indexPickedLJ(i) !THE index number of picked CL sites in the molecular model setting
         enddo
         MinFFLJ = 1.0e3*MinFFLJ/SCALEENERGY/Rg !DIMENTIONLESS THE UNIT
         CLStepRange = CLStepRange/SCALELENGTH !DIMENTIONLESS THE UNIT
         MaxPenalty = -PNTfraction*MinFFLJ ! The maximum penalty potential for changing the charges' location
         MaxCLStep  =  CLStepRange/2.0e0
       endif
       IntraE = 0.0e0
       Lequi = Lequi0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!      PARAMETERS FOR NUMBER OF CYCLES, boundary condition
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
       read(54,'(L)')pbcx
       read(54,'(L)')pbcy
       read(54,'(L)')pbcz
       read(54,*)Ncycle
       read(54,*)Ncycles
       read(54,*)AccRatio
       
       Nmove = 1000
           
       if(SDisplace)then
          IF(HorizontalBin)THEN !Bin divided along the z- direction
             read(54,*)MaxLayer
             read(54,*)LayerH  
             read(54,*)InterH
             SubBox = MaxLayer + 2
             LayerH = LayerH/SCALELENGTH
             InterH = InterH/SCALELENGTH
             ALLOCATE(LayerZL(SubBox+1),LayerZH(SubBox+1), SubVol(SubBox+1))
             ALLOCATE(LayerYL(SubBox+1),LayerYH(SubBox+1))
             ALLOCATE(IniStepX(SubBox+1),IniStepY(SubBox+1),IniStepZ(SubBox+1))
             ALLOCATE(SStepX(SubBox+1),SStepY(SubBox+1),SStepZ(SubBox+1))
             ALLOCATE(NinSBox(SubBox+1),SeqNinSBox(SubBox+1,50000),totNinSBox(SubBox+1),ave_NinSBox(SubBox+1))
             ALLOCATE(SubSE(SubBox+1), SubSEuMolPerM2(SubBox+1))
             ALLOCATE(AtpinMove(SubBox+1), SucinMove(SubBox+1),SucMovein(SubBox+1),SucMoveout(SubBox+1))
             LayerZL(1) = 0.0E0
             LayerZH(1) =  LayerZL(1)+LayerH
             do i = 2, MaxLayer
               LayerZL(i) =  LayerZH(i-1)
               LayerZH(i) =  LayerZL(i)+LayerH             
             enddo
             if(MaxLayer.gt.0)then
               LayerZL(SubBox-1) = LayerZH(MaxLayer)
             else
               LayerZL(SubBox-1) = 0.0e0
             endif
             LayerZH(SubBox-1) = LayerZL(SubBox-1)+InterH
             LayerZL(SubBox) = LayerZH(SubBox-1)
             LayerZH(SubBox) = BoxLengthZ
             
             if(ExtraBin)then
                SubBox=SubBox+1
                read(61,*)ExBinZ
                ExBinZ = ExBinZ/SCALELENGTH
                IF( Bojan .AND. BojanSlit )THEN
                   LayerZL(SubBox)=BoxLengthZ + dble(StripLayerN(1)-1.0)*Stripgap(1)
                ELSE
                   LayerZL(SubBox)=BoxLengthZ
                ENDIF
                LayerZH(SubBox)=BoxLengthZ+ExBinZ
             endif
             do i =1, SubBox
                LayerYL(i)  = 0.0e0
                LayerYH(i)  = CarbonLengthy1
                SubVol(i)   = CarbonLengthx1*CarbonLengthy1*(LayerZH(i)-LayerZL(i))
                IniStepX(i) = CarbonLengthx1/2.0E0
                IniStepY(i) = CarbonLengthy1/2.0E0
                IniStepZ(i) = (LayerZH(i)-LayerZL(i))/2.0e0
             enddo
          ELSEIF(VerticalBin)THEN !Bin divided along the y- direction
             open(unit=24,file='./IniBin.txt',status='unknown')
             rewind(24)
             read(24,*)BinSec
             ALLOCATE(SBoxinBinSec(BinSec),BinSecYL(BinSec),BinSecYH(BinSec))
             ALLOCATE(BinSecZL(BinSec),BinSecZH(BinSec))
             read(24,*)(BinSecYL(i),i=1,BinSec)
             read(24,*)(BinSecYH(i),i=1,BinSec)
             read(24,*)(SBoxinBinSec(i),i=1,BinSec)
             read(24,*)(BinSecZL(i),i=1,BinSec)
             read(24,*)(BinSecZH(i),i=1,BinSec)
             SubBox = 0
             do i = 1, BinSec
                SubBox = SubBox + SBoxinBinSec(i)
                BinSecYL(i) = BinSecYL(i)/SCALELENGTH
                BinSecYH(i) = BinSecYH(i)/SCALELENGTH
                BinSecZL(i) = BinSecZL(i)/SCALELENGTH
                BinSecZH(i) = BinSecZH(i)/SCALELENGTH
             enddo
             ALLOCATE(LayerYL(SubBox+1),LayerYH(SubBox+1), SubVol(SubBox+1))
             ALLOCATE(LayerZL(SubBox+1),LayerZH(SubBox+1))
             ALLOCATE(IniStepX(SubBox+1),IniStepY(SubBox+1),IniStepZ(SubBox+1))
             ALLOCATE(SStepX(SubBox+1),SStepY(SubBox+1),SStepZ(SubBox+1))
             ALLOCATE(NinSBox(SubBox+1),SeqNinSBox(SubBox+1,50000),totNinSBox(SubBox+1),ave_NinSBox(SubBox+1))
             ALLOCATE(AtpinMove(SubBox+1), SucinMove(SubBox+1),SucMovein(SubBox+1),SucMoveout(SubBox+1))
             iSub = 0
             DO i = 1, BinSec
                LayerH = (BinSecYH(i) - BinSecYL(i))/dble(SBoxinBinSec(i))
                do j = 1, SBoxinBinSec(i)
                   iSub = iSub + 1
                   LayerZL(iSub) =  BinSecZL(i)
                   LayerZH(iSub) =  BinSecZH(i)
                   if(j.eq.1)then
                      LayerYL(iSub) = BinSecYL(i)
                      LayerYH(iSub) = LayerYL(iSub) + LayerH
                   else 
                      LayerYL(iSub) = LayerYH(iSub-1)
                      LayerYH(iSub) = LayerYL(iSub) + LayerH
                   endif
                enddo
             ENDDO
             do i =1, SubBox
                SubVol(i) = CarbonLengthx1*(LayerYH(i)-LayerYL(i))*(LayerZH(i)-LayerZL(i))
                IniStepX(i) = CarbonLengthx1/2.0E0
                IniStepY(i) = (LayerYH(i)-LayerYL(i))/2.0E0
                IniStepZ(i) = (LayerZH(i)-LayerZL(i))/2.0e0
             enddo
          ENDIF
          
             open(unit=23,file='./Move-Cycle.txt',status='unknown')
             rewind(23)
             write(23,*)'---------------------------------------------------------------------------'
             write(23,'(4x,4A18)') 'Npart/Pa','NumConfig', 'SucMovein(1-i)', 'SucMoveout(1-i)'
             write(23,*)'---------------------------------------------------------------------------'

             open(unit=16,file='StepLength.txt',status='unknown')
             rewind(16)
             write(16,*)'---------------------------------------------------------------------------'
             write(16,'(4x,4A14)') 'Npart/Pa','StepLX(1-i)A', 'StepLY(1-i)A', 'StepLZ(1-i)A'
             write(16,*)'---------------------------------------------------------------------------'
       endif
       
       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!      INPUTS: ADSORBENT PARAMETERS
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         
!    ------------------------------
!    ASSIGN THE POSITION FOR SOILD
!    ------------------------------      
       if(inclinepore )then
          open(unit=111, file='./InclineBojan.txt',status='unknown')
          read(111,*)AvePoreW
          read(111,*)IncAngle
          read(111,*)IncLayerN
          read(111,*)Incgap
          read(111,*)IncsigmaSS
          read(111,*)IncwelldepthSS
          read(111,*)IncrhosperM2
          AvePoreW = AvePoreW/Scalelength
          IncAngle = IncAngle/180.0e0*pi
          IncCornerZ = AvePoreW-tan(IncAngle)*CarbonLengthy1/2.0e0
          Incgap = Incgap/Scalelength
          IncsigmaSS = IncsigmaSS/Scalelength
          IncwelldepthSS = IncwelldepthSS/ScaleEnergy
          Inc_A = tan(IncAngle)
          Inc_B = -1.0e0
          Inc_C = IncCornerZ
          IncMid_y=CarbonLengthy1/2.0e0
          IncMid_z=IncCornerZ+IncMid_y*tan(IncAngle)
          BoxlengthZ = IncCornerZ + CarbonLengthy1*tan(IncAngle)
          IncLengthY = CarbonLengthy1/cos(IncAngle)
          IncrhosperM2 = IncrhosperM2*((SCALELENGTH*1.0E-10)**2)
       endif
       
      if(Close1End.OR.Close2Ends)then
!         Clo1LengthZ = BoxLengthZ
!         if(inclinepore)Clo1LengthZ = IncCornerZ
         Clo1Mid_z = Clo1LocZ + Clo1LengthZ/2.0e0
         if(inclinepore)then
            Clo1LengthZ = IncCornerZ
            Clo1Mid_z = Clo1LengthZ/2.0e0
         endif
!         if(Clo1Extend)then
!           if(BojanSlit)then
!           Clo1LengthZ = Clo1LengthZ + 2.0e0*abs(StripZ0(1)-(StriplayerN(1)-1)*Stripgap(1))
!           else
!           Clo1LengthZ = Clo1LengthZ + abs(StripZ0(1)-(StriplayerN(1)-1)*Stripgap(1))
!           endif
!         endif
      endif 
      
      if(Close2Ends)then
!         Clo2LengthZ = BoxLengthZ
         Clo2Mid_z = Clo2LocZ + Clo2LengthZ/2.0e0
      endif     

!       if(inclinepore )then
!          open(unit=11, file='firstsurface.txt')
!          read(11,'(f20.10)')CarbonLengthx1
!          read(11,'(f20.10)')CarbonLengthy1
!          read(11,*)NC1
!          do i = 1, NC1
!          read(11,'(3f16.7)')xc(i),yc(i),zc(i)
!          enddo
       
!          open(unit=12, file='secondsurface.txt')
!          read(12,'(f20.10)')CarbonLengthx2
!          read(12,'(f20.10)')CarbonLengthy2
!          read(12,'(f20.10)')PoreWidth
!          read(12,'(f20.10)')CriLengthZ
!          read(12,'(f20.10)')PoreWidth1
!          read(12,'(f20.10)')BoxLengthZ
!          read(12,'(f20.10)')Surftheta
!          read(12,*)NC2
 !            NC3 = NC1 + NC2
!          do i = NC1+1, NC3
!              read(12,'(3f16.7)')xc(i),yc(i),zc(i)
!          enddo
!      endif
    
      if(surfpore)then
     
         open(unit=31, file='firstsurface.txt')
            read(31,'(f20.10)')CarbonLengthx1
            read(31,'(f20.10)')CarbonLengthy1
            read(31,'(f20.10)')PoreWidth2
            read(31,'(f20.10)')BoxLengthZ
            read(31,*)NC1
!            do i = 1, NC1
!               read(31,'(3f16.7)')xc(i),yc(i),zc(i)
!            enddo
       
         open(unit=32, file='secondsurface.txt')
            read(32,'(f20.10)')CarbonLengthy0
            read(32,'(f20.10)')CarbonLengthx2
            read(32,'(f20.10)')CarbonLengthy2
            read(32,'(f20.10)')PoreWidth
            read(32,'(f20.10)')CriLengthZ
            read(32,'(f20.10)')PoreWidth1
            read(32,'(f20.10)')Surftheta
            read(32,*)NC2
            NC3 = NC1 + NC2
!           do i = NC1+1, NC3
!             read(32,'(3f16.7)')xc(i),yc(i),zc(i)
!           enddo
       endif
       read(61,*)NC3
       Allocate(xC(NC3),yC(NC3),zC(NC3))
       
       IF(DecayESF)THEN
         read(61,*) DecayLY
         read(61,*) DecayHY
         read(61,*) CollectLY
         read(61,*) CollectHY
         DecayLY = DecayLY/SCALELENGTH
         DecayHY = DecayHY/SCALELENGTH
         CollectLY = CollectLY/SCALELENGTH
         CollectHY = CollectHY/SCALELENGTH
       ENDIF
       
       if(Steele .or. NVT .or. (NC3.GT.0))then   
           open(unit=60, file='Solid.txt')
           read(60,*)sigmaSS
           read(60,*)WellDepthSS
        endif
      close(unit=60) 
       
!       IF(inclinepore )then
!          do i = 1, NC1
!             read(11,'(3f16.7)')xc(i),yc(i),zc(i)
!          enddo
!          do i = NC1+1, NC3
!              read(12,'(3f16.7)')xc(i),yc(i),zc(i)
!          enddo
!       ELSE 

       IF(surfpore)then
          do i = 1, NC1
             read(31,'(3f16.7)')xc(i),yc(i),zc(i)
          enddo
          do i = NC1+1, NC3
             read(32,'(3f16.7)')xc(i),yc(i),zc(i)
          enddo
       ELSE
          do i = 1, NC3
             read(61,*)xc(i),yc(i),zc(i)
          enddo
       ENDIF

     
!-----------------------------------------------
!    Store the positions of solid atoms
!-----------------------------------------------  
       write(53,*)'NC3 = ', NC3
       do i = 1, NC3
          write(53,'(3f16.7)')xc(i),yc(i),zc(i)
       enddo      
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
   
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     DIMENSIONLESS THE VARIABLES
! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
      SigmaSS        = SigmaSS/SCALELENGTH
      WellDepthSS    = WellDepthSS/SCALEENERGY
      
 !     if(inclinepore)then   
 !       CarbonLengthx2 = CarbonLengthx2/SCALELENGTH
 !       CarbonLengthy2 = CarbonLengthy2/SCALELENGTH
 !       PoreWidth      = PoreWidth/SCALELENGTH
 !       PoreWidth1     = PoreWidth1/SCALELENGTH
 !       CriLengthZ     = CriLengthZ/SCALELENGTH
 !     endif
 
      if(steele)then
        do i = 1, Nsteele
           Psteele(i) = Psteele(i)/SCALELENGTH
        enddo
      endif
      
      if(surfpore)then
        CarbonLengthy0 = CarbonLengthy0/SCALELENGTH
        CarbonLengthx2 = CarbonLengthx2/SCALELENGTH
        CarbonLengthy2 = CarbonLengthy2/SCALELENGTH
        PoreWidth1     = PoreWidth1/SCALELENGTH
        PoreWidth2     = PoreWidth2/SCALELENGTH
        CriLengthZ     = CriLengthZ/SCALELENGTH
      endif

      do i = 1, NC3
         xc(i) = xc(i)/SCALELENGTH
         yc(i) = yc(i)/SCALELENGTH
         zc(i) = zc(i)/SCALELENGTH
      enddo
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!
!    CALCULATION STARTS HERE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!    ------------------------------------------------------  
!    Set rCutOff which can't larger than BoxLength/2
!    ------------------------------------------------------
     BoxLength = CarbonLengthx1
     if(BoxLength .gt. CarbonLengthy1)BoxLength = CarbonLengthy1
     if(BoxLength .gt. BoxLengthZ)BoxLength = BoxLengthZ
     if(pbcx) rcutoff =  CarbonLengthx1/2.0E0
     if(pbcy .AND. rcutoff.GT.CarbonLengthy1/2.0E0)then
        rcutoff = CarbonLengthy1/2.0E0
     endif
     if(pbcz .AND. rcutoff.GT.BoxLengthZ/2.0E0)then
        rcutoff = BoxLengthZ/2.0E0
     endif
     
     rcutoffset = 5.0E0
          
     if(rcutoff.gt.rcutoffset)then
       write(62,*)'---------------------------------------- '
       write(62,*)'Which one would be chosen as the rcutoff:'
       write(62,*)'---------------------------------------- '
       write(62,*)'1-', rcutoff
       write(62,*)'2-', rcutoffset
!       read(62,*) cutoffopt
       if(NLJ.GT.1)then
          cutoffopt=1
       else
          cutoffopt = 2
       endif
       if(cutoffopt.eq.2)rcutoff = rcutoffset
     endif
     
     IF((.not. pbcx) .and. (.not. pbcy) .and. (.not. pbcz))rcutoff = rcutoffset
     IF(Bojan)rcutoff = CarbonLengthx1/2.0E0
     
     write(62,*)'rCutOff=',rCutOff   
     write(51,*)'rCutOff',rCutOff
     r2CutOff = rCutOff*rCutOff
     
     write(62,*)'---------------------------------------------- '
     write(62,*)'ENTER THE MAXIMUM LJ-DISPLACEMENT STEP (0.5) : '
     write(62,*)'---------------------------------------------- '     
     deltaRxInitial = Carbonlengthx1/2.0
     deltaRyInitial = Carbonlengthy1/2.0
     deltaRzInitial = BoxLengthZ/2.0
!    ----------------------------------------------
!    SOLID-FLUID MOLECULAR PARAMETERS
!    ----------------------------------------------
!    Calculated from the Lorentz-Berthelot(LB) rule
!    ----------------------------------------------
!     sigmaSF     = (sigmaFF + sigmaSS)/2.0E0
!     WellDepthSF = sqrt(WellDepthFF*WellDepthSS)
!    --------------------------------
!    SURFACE AREA AND VOLUME VALUES
!    --------------------------------
     IF(DecayESF)THEN
          TotSurfaceArea  = CarbonLengthx1*(CollectHY-CollectLY)             
          BulkVolume      = CarbonLengthx1*CarbonLengthy1*BoxLengthZ 
          CollectBulkV    = CarbonLengthx1*(CollectHY-CollectLY)*BoxLengthZ 
     ELSE
          TotSurfaceArea  = CarbonLengthx1*(AccVYH-AccVYL)             ! for Surface
          BulkVolume      = CarbonLengthx1*CarbonLengthy1*BoxLengthZ   ! the volume of the simulation box 
     ENDIF
     if(inclinePore)then
         SurfaceArea1    = CarbonLengthx1*CarbonLengthy1
         SurfaceArea2    = CarbonLengthx1*IncLengthY
         TotSurfaceArea  = SurfaceArea1 + SurfaceArea2
     endif 
     if(SurfPore)then
         SurfaceArea1    = CarbonLengthx1*CarbonLengthy1
         SurfaceArea2    = CarbonLengthx2*CarbonLengthy2
         TotSurfaceArea  = SurfaceArea1 + SurfaceArea2
     endif
     
     if(Nsteele .gt. 1)then
        TotSurfaceArea  = Nsteele*CarbonLengthx1*CarbonLengthy1
     endif
     
     if(BojanSlit)then
        TotSurfaceArea = 0.0e0
        do i = 1, NS
           TotSurfaceArea = TotSurfaceArea + 2.0e0*CarbonLengthx1*StripW(i)/cos(StripAngleY(i))
        enddo
     endif
     
     if(Close1End.AND.(.NOT. Clo1Hard))then
        TotSurfaceArea = TotSurfaceArea+opClo1LengthZ*CarbonLengthx1
     endif
     
      if(Close2Ends)then
   !      if(Clo1Hard)then
           TotSurfaceArea = TotSurfaceArea+opClo2lengthZ*CarbonLengthx1
   !      else
   !        TotSurfaceArea = TotSurfaceArea+(opClo1LengthZ+opClo2lengthZ)*CarbonLengthx1
   !      endif
     endif
     
     write(62,*)'------------------------------------- '
     write(62,*)'Calculating the Access Volume:        '
     write(62,*)'------------------------------------- '  
     call accessiblevolume
!     AccVolume = (BoxLengthZ-SCALELENGTH)*CarbonLengthx1*CarbonLengthy1
     AccVolumeM3 = AccVolume*(SCALELENGTH*1.0E-10)**3

     write(57,*)'AccVolumeM3',AccVolumeM3
     
     write(62,1130) TotSurfaceArea, BulkVolume, AccVolume
     write(51,1130) TotSurfaceArea, BulkVolume, AccVolume
1130 format( 1x, 'SurfaceArea(-)--------------------------:', f20.10,/, &
  &          1x, 'BulkVolume(-)---------------------------:', f20.10,/, &
  &          1x, 'AccVolume(-)----------------------------:', f20.10,/)  
     
     if(inclinepore.OR. SurfPore)then
     write(62,1141) SurfaceArea1, SurfaceArea2,TotSurfaceArea,BulkVolume, AccVolume
     write(51,1141) SurfaceArea1, SurfaceArea2,TotSurfaceArea,BulkVolume, AccVolume
1141 format( 1x, 'SurfaceArea1(-)-------------------------:', f20.10,/, &
  &          1x, 'SurfaceArea2(-)-------------------------:', f20.10,/, &
  &          1x, 'TotSurfaceArea(-)-----------------------:', f20.10,/, &
  &          1x, 'BulkVolume(-)---------------------------:', f20.10,/, & 
  &          1x, 'AccVolume(-)----------------------------:', f20.10,/)  
     endif
      
!      -----------------------------------------------    
!      INPUTS: THE INITIAL CONFIGURATION OF THE SYSTEM
!      ----------------------------------------------- 
       Npart = 0    
       open( unit=17,  file='IniConfig.txt', status='unknown')
       rewind(17)
       read(17,*)Npart     
       do i = 1, Npart
          read(17,*) mcx(i), mcy(i), mcz(i)       ! Reduced unit
       enddo  
       read(17,*) 
       do i = 1, Npart
          do j = 1, NLJ
             read(17,*)LJX(i,j),LJY(i,j),LJZ(i,j)
          enddo
       enddo
       read(17,*)
       do i = 1, Npart
          do j = 1, NCL
             read(17,*)CLX(i,j),CLY(i,j),CLZ(i,j)
          enddo
       enddo
!    ---------------
!    DO GCMC OR NVT
!    ---------------
     if(GCMC) call ENSEMBLE_GCMC
     if (NVT) call ENSEMBLE_NVT
     
     DEALLOCATE(Psteele)
     DEAllocate(xC,yC,zC)
     write(62,*) "Exits normally"
      entim = MPI_Wtime()
write(62,*) "Time taken by Rank", rank, ":", entim-pretim
    call MPI_FINALIZE(ierror)
    END PROGRAM PoreProcess6
