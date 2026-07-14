Grid += gr_ptInterface.o gr_ptInit.o gr_ptMove.o \
        gr_ptFinalize.o gr_ptMarkRefineDerefine.o \
        gr_ptMoveOffBlk.o\
        gr_ptFillBlkParticleInfo.o gr_ensureValidNeighborInfo.o \
        gr_solversInit.o gr_solversFinalize.o gr_solversTestInit.o \
	gr_bcInterface.o gr_bcHseInit.o gr_bcApplyToAllBlks.o gr_bcApplyToOneFace.o \
	gr_hgInit.o  gr_mpoleInit.o gr_isoMpoleInit.o gr_pfftInit.o \
	gr_hgFinalize.o  gr_mpoleFinalize.o gr_isoMpoleFinalize.o gr_pfftFinalize.o \
        gr_bcInit.o gr_bcFinalize.o gr_hgInterface.o gr_pfftInterface.o \
	gr_isoInterface.o gr_ptMapInit.o gr_hgPfftInit.o gr_hgPfftFinalize.o \
	gr_interface.o gr_interfaceTypeDecl.o gr_pfftInterfaceTypeDecl.o \
	gr_mgInterface.o gr_mgInit.o gr_mgFinalize.o gr_mgPfftInit.o gr_mgPfftFinalize.o \
        gr_bicgInterface.o gr_bicgInit.o gr_bicgFinalize.o \
        gr_getRegionDataCoordinates.o gr_getBndBox.o

Grid += gr_hg_amr_1blk_bcset_work.o gr_pfftPoissonDirect.o

Grid +=  gr_hypreInterface.o gr_hypreInit.o gr_hypreFinalize.o

Grid +=  gr_amrexMultigridInterface.o gr_amrexMultigridInit.o gr_amrexMultigridFinalize.o

Grid += gr_mpoleCen1Dspherical.o         \
	gr_mpoleCen2Dspherical.o         \
        gr_mpoleCen2Dcylindrical.o       \
	gr_mpoleCen3Dspherical.o         \
        gr_mpoleCen3Dcartesian.o         \
        gr_mpoleCenterOfExpansion.o      \
        gr_mpoleCollectMoments.o         \
        gr_mpoleDeallocateRadialArrays.o \
        gr_mpoleDumpMoments.o            \
        gr_mpoleFinalize.o               \
        gr_mpoleHeapsort.o               \
        gr_mpoleInit.o                   \
        gr_mpoleMom1Dspherical.o         \
	gr_mpoleMom2Dspherical.o         \
        gr_mpoleMom2Dcylindrical.o       \
	gr_mpoleMom3Dspherical.o         \
        gr_mpoleMom3Dcartesian.o         \
        gr_mpoleMomBins1Dspherical.o     \
	gr_mpoleMomBins2Dspherical.o     \
        gr_mpoleMomBins2Dcylindrical.o   \
	gr_mpoleMomBins3Dspherical.o     \
        gr_mpoleMomBins3Dcartesian.o     \
        gr_mpoleMoments.o                \
        gr_mpolePot1Dspherical.o         \
	gr_mpolePot2Dspherical.o         \
        gr_mpolePot2Dcylindrical.o       \
	gr_mpolePot3Dspherical.o         \
        gr_mpolePot3Dcartesian.o         \
        gr_mpolePotentials.o             \
        gr_mpolePrintRadialInfo.o        \
        gr_mpoleRad1Dspherical.o         \
	gr_mpoleRad2Dspherical.o         \
        gr_mpoleRad2Dcylindrical.o       \
	gr_mpoleRad3Dspherical.o         \
        gr_mpoleRad3Dcartesian.o         \
        gr_mpoleRadialSampling.o         \
        gr_mpoleSetInnerZoneGrid.o       \
        gr_mpoleSetOuterZoneGrid.o       \
        gr_mpoleSetRadialBinData.o

Grid += gr_sbInterface.o gr_sbInit.o gr_sbSendBoundBox.o\
	gr_sbSendParticles.o gr_sbCreateParticles.o gr_sbGetProcBlock.o gr_sbSendPosn.o gr_sbStoreParticlesPerProc.o \
	gr_sbSendParticleCount.o gr_sbSendParticles.o gr_sbUpdateForces.o gr_sbSendForces.o gr_sbFinalize.o

Grid += gr_bhInterface.o gr_bhInit.o gr_bhFinalize.o

Grid += gr_xyzToBlockLevel.o gr_xyzToBlock.o

Grid += gr_forceDerefInRadius.o
Grid += gr_markRefineDerefine.o gr_estimateBlkError.o

Grid += gr_blockMatch.o
Grid += gr_ptWritePCs.o

gr_interface.o     : gr_interfaceTypeDecl.o
gr_pfftInterface.o : gr_pfftInterfaceTypeDecl.o
