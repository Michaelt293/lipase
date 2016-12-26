{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TemplateHaskell #-}
module FattyAcid where

import qualified Isotope as I
import Isotope hiding (monoisotopicMass)
import Spectra
import Data.Maybe
import Control.Lens

newtype NumCarbons = NumCarbons { _getNumCarbons :: Int }
  deriving (Show, Eq, Ord, Num)

makeClassy ''NumCarbons

instance ShowVal NumCarbons where
  showVal (NumCarbons n) = show n

newtype NumDoubleBonds = NumDoubleBonds { _getNumDoubleBonds :: Int }
  deriving (Show, Eq, Ord, Num)

makeClassy ''NumDoubleBonds

instance ShowVal NumDoubleBonds where
  showVal (NumDoubleBonds n) = show n

data FattyAcyl = FattyAcyl {
    _fattyAcylCarbons     :: NumCarbons
  , _fattyAcylDoubleBonds :: NumDoubleBonds
} deriving (Eq, Ord)

makeClassy ''FattyAcyl

instance HasNumCarbons FattyAcyl where
  numCarbons = fattyAcylCarbons

instance HasNumDoubleBonds FattyAcyl where
  numDoubleBonds = fattyAcylDoubleBonds

instance Show FattyAcyl where
  show (FattyAcyl cs dbs) = showVal cs ++ ":" ++ showVal dbs

instance ToMolecularFormula FattyAcyl where
  toMolecularFormula (FattyAcyl cs dbs) = mkMolecularFormula
    [ (C, _getNumCarbons cs)
    , (H, _getNumCarbons cs * 2 + 1 - _getNumDoubleBonds dbs * 2)
    , (O, 2)
    ]

instance ToElementalComposition FattyAcyl where
  toElementalComposition = toElementalComposition . toMolecularFormula

c10_0 = FattyAcyl 10 0

c12_0 = FattyAcyl 12 0

c14s = FattyAcyl 14 <$> [0, 1]

c15s = FattyAcyl 15 <$> [0, 1]

c16s = FattyAcyl 16 <$> [0, 1, 2]

c17s = FattyAcyl 17 <$> [0, 1]

c18s = FattyAcyl 18 <$> [0, 1, 2, 3, 4]

c19s = FattyAcyl 19 <$> [0, 1]

c20s = FattyAcyl 20 <$> [0, 1, 2, 3, 4, 5]

c21s = FattyAcyl 21 <$> [0, 1]

c22s = FattyAcyl 22 <$> [0, 1, 2, 3, 4, 5, 6]

evenChainFAs :: [FattyAcyl]
evenChainFAs = [c10_0, c12_0] ++ c14s ++ c16s ++ c18s ++ c20s ++ c22s

oddChainFAs :: [FattyAcyl]
oddChainFAs = c15s ++ c17s ++ c19s ++ c21s

fattyAcyls :: [FattyAcyl]
fattyAcyls = evenChainFAs ++ oddChainFAs

faMonoisotopicMass :: FattyAcyl -> MonoisotopicMass
faMonoisotopicMass fa = I.monoisotopicMass fa |+| I.monoisotopicMass H

fattyAcidMonoisotopicMasses :: [MonoisotopicMass]
fattyAcidMonoisotopicMasses = faMonoisotopicMass <$> fattyAcyls

data AssignedFA = AssignedFA {
    _getAssignedFA :: Maybe FattyAcyl
  , _assignedFAIonInfo :: IonInfo
} deriving (Show, Eq, Ord) -- Maybe write my own Show instance.

makeClassy ''AssignedFA

instance HasIonInfo AssignedFA where
  ionInfo = assignedFAIonInfo

data AssignedFAs = AssignedFAs {
    _assignedFAsPrecIon :: Mz
  , _getAssignedFAs :: [AssignedFA]
} deriving (Eq, Ord)

makeClassy ''AssignedFAs

instance HasMonoisotopicMass AssignedFAs where
  monoisotopicMass = assignedFAsPrecIon

instance Show AssignedFAs where
  show (AssignedFAs p fas) =
    "Precursor ion: " ++ show p ++ "\n" ++
    showList' fas

assignFA :: Double -> MonoisotopicMass -> [FattyAcyl] -> Maybe FattyAcyl
assignFA n m fas =
  case fas of
    [] -> Nothing
    f:fs -> if withinTolerance n m (faMonoisotopicMass f)
              then Just f
              else assignFA n m fs

neutralLossRowToFA :: NeutralLossRow -> AssignedFA
neutralLossRowToFA (NeutralLossRow nl i) =
  AssignedFA (assignFA 0.3 nl fattyAcyls) i

toAssignedFAs :: NeutralLossSpectrum -> AssignedFAs
toAssignedFAs (NeutralLossSpectrum p nls) =
  AssignedFAs p (neutralLossRowToFA <$> nls)

collectFAs :: AssignedFAs -> [FattyAcyl]
collectFAs fas =
  catMaybes $ _getAssignedFA <$> fas ^. getAssignedFAs
