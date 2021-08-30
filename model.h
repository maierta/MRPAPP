#ifndef MODEL_H
#define MODEL_H

#ifdef USE_SRRUO
#include "SrRuO_SO.h"
#elif USE_SRRUO3D
#include "SrRuO_SO_3D.h"
#elif USE_1BANDWSPIN
#include "1band_wSpin.h"
#elif USE_BILAYER_FESC
#include "bilayerFESC.h"
#elif USE_BILAYER_1BAND
#include "bilayer.h"
#elif USE_ORTHOIIBILAYER
#include "orthoIIBilayer.h"
#elif USE_BSCCOBILAYER
#include "BSCCObilayer.h"
#elif USE_BILAYER_FESC
#include "bilayer.h"
#elif USE_BAFEAS
#include "BaFeAs_5orb.h"
#elif USE_KFE2SE2
#include "KFe2Se2.h"
#elif USE_FOURORBITAL
#include "FourOrbital.h"
#elif USE_TBFILE
#include "tbFromFile.h"
#elif USE_COUPLEDLADDERS
#include "coupledLadders.h"
#elif USE_NDNIO2
#include "NdNiO2.h"
#endif
#elif USE_MODELFROMFILESO
#include "modelFromFileSO.h"
#endif

