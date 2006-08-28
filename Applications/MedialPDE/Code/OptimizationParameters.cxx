#include "OptimizationParameters.h"

OptimizationParameters
::OptimizationParameters()
{
  // Initialize the enum mappings if necessary
  xOptimizerRegMap.AddPair(CONJGRAD, "ConjGrad");
  xOptimizerRegMap.AddPair(GRADIENT, "Gradient");
  xOptimizerRegMap.AddPair(EVOLUTION, "EvolStrat");

  xMappingRegMap.AddPair(AFFINE, "Affine");
  xMappingRegMap.AddPair(COARSE_TO_FINE, "CoarseFine");
  xMappingRegMap.AddPair(IDENTITY, "Identity");
  xMappingRegMap.AddPair(PCA, "PCA");

  xImageMatchRegMap.AddPair(VOLUME, "VolumeOverlap");
  xImageMatchRegMap.AddPair(BOUNDARY, "BoundaryIntegral");

  xPenaltyTermRegMap.AddPair(BOUNDARY_JACOBIAN, "BoundaryJacobianEnergyTerm");
  xPenaltyTermRegMap.AddPair(MEDIAL_REGULARITY, "MedialRegularityTerm");
  xPenaltyTermRegMap.AddPair(ATOM_BADNESS, "AtomBadnessTerm");
  xPenaltyTermRegMap.AddPair(RADIUS, "RadiusPenaltyTerm");

  xCTFSettingsRegMap.AddPair(COSINE_BASIS_PDE, "CosineBasisPDE");
  xCTFSettingsRegMap.AddPair(LOOP_SUBDIVISION_PDE, "LoopSubdivisionPDE");
  xCTFSettingsRegMap.AddPair(NONE, "NONE");

  // Clear the settings pointer
  xCTFSettings = NULL;

  // Set default values
  xOptimizer = CONJGRAD;
  xImageMatch = VOLUME;
  xMapping = IDENTITY;
}

OptimizationParameters
::~OptimizationParameters()
{
  if(xCTFSettings != NULL) delete xCTFSettings;
}

void
OptimizationParameters
::ReadFromRegistry(Registry &R)
{
  // Read the enums
  xOptimizer = R["Optimizer"].GetEnum(xOptimizerRegMap, CONJGRAD);
  xImageMatch = R["ImageMatch"].GetEnum(xImageMatchRegMap, VOLUME);
  xMapping = R["Mapping"].GetEnum(xMappingRegMap, IDENTITY);

  // Read the weights of the different terms in the optimization
  xTermWeights[BOUNDARY_JACOBIAN] = R[xPenaltyTermRegMap(BOUNDARY_JACOBIAN)][0.5];
  xTermWeights[MEDIAL_REGULARITY] = R[xPenaltyTermRegMap(MEDIAL_REGULARITY)][1.0];
  xTermWeights[ATOM_BADNESS] = R[xPenaltyTermRegMap(ATOM_BADNESS)][0.01];
  xTermWeights[RADIUS] = R[xPenaltyTermRegMap(RADIUS)][0.1];

  // Read the coarse-to-fine specification
  if(xCTFSettings != NULL) delete xCTFSettings;
  CTFSettings c2f = R["CoarseToFine.Mode"].GetEnum(xCTFSettingsRegMap, NONE);
  if(c2f == COSINE_BASIS_PDE)
    xCTFSettings = new FourierCoarseToFineSettings();
  else xCTFSettings = NULL;

  // Read the coarse to fine settings if they are there
  if(xCTFSettings != NULL)
    xCTFSettings->ReadFromRegistry(R.Folder("CoarseToFine"));

  // Read the PCA settings if they are there
  xPCAFileName = R["PCA.FileName"][""];
  nPCAModes = R["PCA.NumberOfModes"][10];
}