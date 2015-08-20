 MODULE Constant_M
           real*8 pi, AvogadroNumber, kB, Rg, hv, permittivity, unite
           real*8 deBro, deBro3
           real*8 GrapheneLayer, eV
     END MODULE Constant_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
     MODULE FLUID_SOLID_M
          integer NLJ, NCL, NS
          integer Nsteele
          integer SubFGas, maxAngleBin, fSite, tSite
          integer StriplayerN(100)
         
          real*8 MW, VectorL, deltaAngle   
          real*8 MCx0, MCy0, MCz0
          real*8 MCx(50000), MCy(50000), MCz(50000)   ! MCx(Npart)
          real*8 LJx0(20), LJy0(20), LJz0(20)         ! LJx0(NLJ) 
          real*8 LJx(50000,20),  LJy(50000,20), LJz(50000,20)  !LJx(Npart,NLJ)
          real*8 CLx0(20), CLy0(20), CLz0(20)          ! CLx0(NCL)
          real*8 CLx(50000,20), CLy(50000,20), CLz(50000,20)   ! CLx(Npart,NCL) 
          real*8 sigmaFF(20), wellDepthFF(20),Charges(20)          ! sigmaFF(NLJ)! Charges(NCL)
          real*8 WellDepthSS, sigmaSS  
          real*8 BsigmaSS(100),BwelldepthSS(100)
          real*8 BsigmaSF,BdeltaZ,Py,Ny,subforce
          real*8 StripZ0(100),StripC(100),StripW(100),Stripgap(100),rhosperM2(100)
          real*8 thetav, thetar, ChemPD0, ChemPDe
          
          
          real*8, Dimension(:),POINTER::Psteele 

          logical steele, Bojan, BojanSlit
          logical TopStrip(100), OrienD, SymVector   
     END MODULE FLUID_SOLID_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     MODULE PHYSICAL_M
           integer Npart, AddN, NumUBin
           
           real*8 Temperature, T, CriticalTemp, Psaturate
           real*8 BulkVolume, AccVolume, AccVolumeM3
           real*8 EntropySys,LowUi,HighUi,LowUsys,HighUsys
           real*8, Dimension(:),POINTER::rho, P, Pressure, SEinput
           real*8, Dimension(:),POINTER::idealChemicalPotential, ChemicalPotential
           real*8, Dimension(:),POINTER::ResidualChemicalPotential, ChemicalPotentialJoulepermol
           real*8, Dimension(:),POINTER::UBinL, UBinH, UBinTime, ProbaUBin
           real*8, Dimension(:,:),POINTER::EnergyMatrix,VirialMatrix,DerivEMatrix, DerivEMatrix_Ale
           
     END MODULE PHYSICAL_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     MODULE MCSETTING_M
            integer Ncycle, Ncycles, Nmove, CountN
            integer NGhost, NPartinY, NumOfrho
            integer iden, TotPickN, OutZPickN
            
            real*8 a(8), b(6), G(6), rr
            real*8 Rot(9)
            real*8 rCutOff, r2CutOff
            real*8 TotSurfaceArea, RatioV1, RatioV2
            real*8 SCALELENGTH, SCALEENERGY
            real*8 Virial, VirialC, VC, TVirial, VirialCL, BulkVirial, BulkVirialF
            real*8 DerivE, FirstDerivE
            real*8 Dist_Limit, PairEff_Limit, PairEsf_Limit, PairV_Limit, DerivE_Limit
            real*8 DerivE_Ale, FirstDerivE_Ale
            real*8 BoxLengthX, BoxLengthY, BoxLengthZ, minBoxLength
            real*8 AccVYL, AccVYH, AccVZL, AccVZH
            real*8 Part1Ly, Part2Ly, Part3Ly, Part2Lz, StopZL, StopZH, SingleKineticE
            real*8 InsertZL, InsertZH, ExtraBoxH
            
            logical pbcx,pbcy,pbcz, BoPBCy, GCMC 
            logical NVT, SURFACE, PORE, EQUILIBRIUM, KINETIC, NVTDes, ExtraBox
            logical kVLE, WIDOM_VLE, VLEradius, VLEdistance, ChangeXZ, ChangeXY
            logical SolidSlab, HomoStart, NBOURBIN, MassTSFER, ENTROPY, MeanFPath
            logical DistanceD, LocalD2D
            logical Close1End,Close2Ends, Clo1Hard, Clo1Extend

     END MODULE MCSETTING_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    MODULE POREFIGURE_M
            integer IncLayerN
            integer Clo1LayerN,Clo2LayerN
     
            real*8 IncCornerZ, IncAngle, Incgap
            real*8 IncsigmaSS,IncwelldepthSS
            real*8 Inc_A, Inc_B, Inc_C
            real*8 IncMid_y,IncMid_z,IncLengthY,IncrhosperM2
            real*8 Clo1LocZ,Clo1gap,Clo1rhosperM2
            real*8 Clo1sigmaSS,Clo1welldepthSS
            real*8 Clo1Mid_z, Clo1LengthZ, opClo1LengthZ
            real*8 Clo2LocZ,Clo2gap,Clo2rhosperM2
            real*8 Clo2sigmaSS,Clo2welldepthSS
            real*8 Clo2Mid_z, Clo2LengthZ, opClo2LengthZ
                              
     END MODULE POREFIGURE_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    MODULE SUBBOX_M
            integer SubBox, MaxLayer, InterLayer
            integer LiqBin, GasBin
            integer indNinSBox(50000)
            integer NinPropSub(3), SeqNinPropSub(3,50000), indNinPropSub(50000)
            integer NunitX, NunitY, NunitZ, NFluxSurf
            integer NinSBox(50000)
            
            real*8  rCutoffProp(3),r2CutoffProp(3)
            real*8  VirialProp(3), VirialF(3), EnergyProp(3), EnergyF(3)
            real*8  totEinBulk, EinBulk, EinBulkF 
            real*8  ExtTEnergy, ExtEnergy, ExtCLEnergy, ExtEnergy1, ExtCLEnergy1
            real*8  LiqDens, VLEBoxlengthX, VLEBoxlengthY, VLEBoxlengthZ, LiqLY, LiqHy
            real*8  VLEBoxlengthX0, VLEBoxlengthY0, VLEBoxlengthZ0
            real*8  PropLy(3), PropHy(3),AveProprho(3)
            real*8  LunitX, LunitY, LunitZ
            real*8  TotRate, InterTime, TotTime
            real*8, Dimension(:),POINTER:: LayerXL,LayerXH,LayerYL,LayerYH,LayerZL,LayerZH
            real*8, Dimension(:),POINTER:: SubVol,totNinSBox
            real*8, Dimension(:),POINTER:: IniStepX,IniStepY,IniStepZ
            real*8, Dimension(:),POINTER:: SStepX,SStepY,SStepZ
            real*8, Dimension(:),POINTER:: ave_NinSBox, rho_SBox, rho_SBoxMolPerM3
            real*8, Dimension(:),POINTER:: AtpinMove, SucinMove,SucMovein,SucMoveout
            real*8, Dimension(:),POINTER:: rhoFW,ChemicalPVLE, ExcessChemicalPVLE, IdealChemicalPVLE
            real*8, Dimension(:),POINTER:: ChemPVLETotal, ExChemPVLETotal, IdChemPVLETotal
            real*8, Dimension(:),POINTER:: ChemPVLEAve, ExChemPVLEAve, IdChemPVLEAve
            real*8, Dimension(:),POINTER:: ChemPVLEJPerMole, ExChemPVLEJPerMole, IdChemPVLEJPerMole
            real*8, Dimension(:),POINTER:: rate, AccuR,Energyi, PairEnergyi
            real*8, Dimension(:),POINTER:: UpFlux,DownFlux,NetFlux,FluxSurf
            real*8, Dimension(:,:),POINTER:: SeqNinSBox
            real*8, Dimension(:),POINTER:: rateSBox, rateTSBox, timeSBox, ChemPSBox 
                 
     END MODULE SUBBOX_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       MODULE ROTATE_M
            real*8, Dimension(:,:),POINTER::RLJx,RLJy,RLJz
            real*8, Dimension(:),POINTER::RotEnergyi,Rotrate,AccuRotrate
            real*8, Dimension(:,:),POINTER::RCLx, RCLy, RCLz
            real*8, Dimension(:,:),POINTER::RPairEnergyNew
       END MODULE ROTATE_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    MODULE F2DLOCALD_M
            integer totBin2D, unitNo
            integer Sec2D, Bin2DSec(5)
            integer NinXSec2D, NinYSec2D(5),NinZSec2D(5)
            integer, Dimension(:),POINTER:: Nin2DBin, NBinfAve
            integer, Dimension(:,:),POINTER:: SeqBinfAve
            
            real*8 AveRadius
            real*8 deltaX2D, deltaY2D(5), deltaZ2D(5)  ! the distance for dividing the 2D bin for each section 2DdeltaY(2DSec)
            real*8 SecYL2D(5), SecYH2D(5), SecZL2D(5), SecZH2D(5) !The boundaries to define the section in dividing the 2D bin
            real*8 BinVolSec(5)
            real*8, Dimension(:),POINTER::Sum_Nin2DBin, ave_Nin2DBin
            real*8, Dimension(:),POINTER::Bin2DX,Bin2DY,Bin2DZ
            real*8, Dimension(:),POINTER::DDensity2D,Smooth_DDensity2D,DDensity2DKmolperM3
            real*8, Dimension(:),POINTER::SumD_Bin2D, Bin2DVol, SumV_Bin2D
           
 
    END MODULE F2DLOCALD_M
