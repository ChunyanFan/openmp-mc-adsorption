
     MODULE Constant_M
           real*8 pi, AvogadroNumber, kB, Rg, hv, permittivity, unite
           real*8 deBro, deBro3
           real*8 GrapheneLayer
     END MODULE Constant_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
     MODULE FLUID_SOLID_M
          integer NLJ, NCL, NC3, NS
          integer Nsteele
         
          real*8 MW   
          real*8 MCx0, MCy0, MCz0
          real*8 MCx(50000), MCy(50000), MCz(50000)   ! MCx(Npart)
          real*8 LJx0(20), LJy0(20), LJz0(20)         ! LJx0(NLJ) 
          real*8 LJx(50000,20),  LJy(50000,20), LJz(50000,20)  !LJx(Npart,NLJ)
          real*8 CLx0(20), CLy0(20), CLz0(20)          ! CLx0(NCL)
          real*8 CLx(50000,20), CLy(50000,20), CLz(50000,20)   ! CLx(Npart,NCL) 
          real*8 ECLx(50000,20), ECLy(50000,20), ECLz(50000,20)  ! ECLx(Npart,NCL) Equilibrium position of Charges in Molecular 
          real*8 sigmaFF(20), wellDepthFF(20),Charges(20)          ! sigmaFF(NLJ)! Charges(NCL)
          real*8 WellDepthSS, sigmaSS  
          real*8 BsigmaSS(100),BwelldepthSS(100)
          real*8 StripZ0(100),StripC(100),StripW(100),Stripgap(100),StriplayerN(100),rhosperM2(100), StripAngleY(100)
          real*8 StripWCor(100), OutLzLyStripZ(100), OutHzLyStripZ(100), StripLy(100), StripHy(100)
          real*8 InLzLyStripZ(100), InHzLyStripZ(100)
          real*8 tan_StripAngleY(100), sin_StripAngleY(100), cos_StripAngleY(100)
          real*8 BsigmaSF,BdeltaZ,Py,Ny,subforce

          real*8, Dimension(:),POINTER::xC, yC, zC
          real*8, Dimension(:),POINTER::Psteele     
          
          logical steele, Bojan, BojanSlit   
          logical TopStrip(100)  
     END MODULE FLUID_SOLID_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     MODULE PHYSICAL_M
           integer Npart, addN
           
           real*8 Temperature, T, CriticalTemp, Psaturate
           real*8 AccVolume, BulkVolume, BBulkVolume, AccVolumeM3
           real*8, Dimension(:),POINTER::rho, P, Pressure, zactivity, SEinput, InrefChemP
           real*8, Dimension(:),POINTER::idealChemicalPotential, ChemicalPotential
           real*8, Dimension(:),POINTER::ResidualChemicalPotential, ChemicalPotentialJoulepermol
           real*8, Dimension(:),POINTER::Fun_NN_Bulk, N_Bulk, Fun_EN_Bulk, CFun_NN_Bulk
           
     END MODULE PHYSICAL_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     MODULE MCSETTING_M
            integer iPC, iden, NumOfPre, NumOfrho
            integer Ncycle, Ncycles, Nmove, Ninpore, CollectN
            
            real*8 a(8), b(6), G(6), rr
            real*8 Rot(9)
            real*8 rCutOff, r2CutOff
            real*8 deltaRx,deltaRy,deltaRz
            real*8 deltaRxInitial,deltaRyInitial,deltaRzInitial
            real*8 TotSurfaceArea
            real*8 SCALELENGTH, SCALEENERGY
            real*8 PoreWidth,BoxLengthZ,CarbonLengthx1,CarbonLengthy1,BoxLength
            real*8 IniBoxLengthZ,IniCarbonLengthx1,IniCarbonLengthy1
            real*8 AccVYL, AccVYH, AccVZH
            real*8 Surftheta,PoreWidth1, CriLengthZ, PoreWidth2,CarbonLengthy0
            real*8 Virial, VirialC, VC, BulkVirial, TVirial, VirialCL, VirialF
            real*8 SMLimitZ,SMLimitZ2,ksi,ksiF,SFenergy(100000)
            real*8 CollectBulkV, DecayLY, DecayHY, CollectLY, CollectHY
            
            logical JOHNSON, WAGNER            
            logical inclinePore,Close1End,Close2Ends, Clo1Hard, Clo1Extend
            logical SurfPore
            logical GCMC, NVT, NPT, BULKGCMC, INPUTCP
            logical pbcx,pbcy,pbcz,BoPBCy  
            logical SMediation1,SMediation2
            logical NaddNVT, NVTDes
            logical SDisplace, TradMC, ExtraBin,ChangeXY
            logical DistanceD, RadiusD, DecayESF, LocalD2D
            logical HorizontalBin, VerticalBin

     END MODULE MCSETTING_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     MODULE BUCKINGHAM_M
            real*8 alpha_B, rM, rMin
            real*8 WelldepthFF_B, SigmaFF_B, AFF_B, BFF_B, alphaFF_B
            real*8 WelldepthSS_B, ETA, ASS_B, BSS_B, alphaSS_B
            real*8 ASF_B, BSF_B, alphaSF_B
            
            logical BUCKINGHAM, CrowellChang
     END MODULE BUCKINGHAM_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     MODULE FLUCTUATION_M
            integer NinBin(1000),maxbin 
            integer I_POS(10000),SeqNinBin(1000,1000)
            integer SepLDN, MLBin
            integer maxAngleBin, fSite, tSite
     
            real*8 ENERGYMATRIX(10000,10000), layerB(5)
            real*8 sum_N, sum_NN, CountN
            real*8 deltaZ, VectorL, deltaAngle
            real*8 y_rhoBnd(10), Sep_Lengthy(10)
            real*8 yL_ComBnd(10), yH_ComBnd(10), zL_ComBnd(10), zH_ComBnd(10)
            real*8 Accu_mcBin_ML(1000),Accu_mcBin_SL(1000)
            
            real*8, Dimension(:),POINTER:: UFFinBin,USFinBin
              
            
            logical LOCALF, ADSORPTION, LAYERF, FINALLD  
            logical Sep_rho_Dis, CorCompress, OrienD, SymVector        
     END MODULE FLUCTUATION_M
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
    MODULE FLEXICLINFO_M
            integer NpickedCL
            integer indexPickedCL(20),indexPickedLJ(20)
     
            real*8 MinFFLJ, PNTfraction, CLStepRange, RdmCL 
            real*8 MaxPenalty, MaxCLStep
            real*8 IntraE(10000) !IntraE(NPart)
            real*8 Lequi0,Lequi(10000)
            
            logical FlexiCL                  
     END MODULE FLEXICLINFO_M
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    MODULE SUBBOX_M
            integer SubBox, MaxLayer, BinSec
            integer, Dimension(:),POINTER:: SBoxinBinSec
     
            real*8  totEinBulk, EinBulk, AccRatio
            real*8  indNinSBox(50000)
            real*8  ExtTEnergy, ExtEnergy, ExtCLEnergy, ExtEnergy1, ExtCLEnergy1
            real*8, Dimension(:),POINTER:: LayerZL,LayerZH, SubVol
            real*8, Dimension(:),POINTER:: IniStepX,IniStepY,IniStepZ
            real*8, Dimension(:),POINTER:: SStepX,SStepY,SStepZ
            real*8, Dimension(:),POINTER:: NinSBox,totNinSBox,ave_NinSBox
            real*8, Dimension(:),POINTER:: SubSE, SubSEuMolPerM2
            real*8, Dimension(:),POINTER:: AtpinMove, SucinMove,SucMovein,SucMoveout
            real*8, Dimension(:),POINTER:: BinSecYL, BinSecYH, BinSecZL, BinSecZH
            real*8, Dimension(:),POINTER:: LayerYL,LayerYH
            real*8, Dimension(:,:),POINTER:: SeqNinSBox
                 
     END MODULE SUBBOX_M
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
