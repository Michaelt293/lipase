{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TemplateHaskell #-}
module FattyAcid where

import Isotope
import Spectra
import Data.Maybe
import Control.Lens

newtype NumCarbons = NumCarbons { getNumCarbons :: Int }
  deriving (Eq, Ord, Num)

instance Show NumCarbons where
  show (NumCarbons n) = show n

newtype NumDoubleBonds = NumDoubleBonds { getNumDoubleBonds :: Int }
  deriving (Eq, Ord, Num)

instance Show NumDoubleBonds where
  show (NumDoubleBonds n) = show n

data FattyAcyl = FattyAcyl {
    _numCarbons     :: NumCarbons
  , _numDoubleBonds :: NumDoubleBonds
} deriving (Eq, Ord)

makeClassy ''FattyAcyl

instance Show FattyAcyl where
  show (FattyAcyl cs dbs) = show cs ++ ":" ++ show dbs

instance ToMolecularFormula FattyAcyl where
  toMolecularFormula (FattyAcyl cs dbs) = mkMolecularFormula
    [ (C, getNumCarbons cs)
    , (H, getNumCarbons cs * 2 + 1 - getNumDoubleBonds dbs * 2)
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
faMonoisotopicMass fa = monoisotopicMass fa |+| monoisotopicMass H

fattyAcidMonoisotopicMasses :: [MonoisotopicMass]
fattyAcidMonoisotopicMasses = faMonoisotopicMass <$> fattyAcyls

data TentativelyAssignedFA = TentativelyAssignedFA {
    _getTentativelyAssignedFA :: Maybe FattyAcyl
  , _tentativelyAssignedFAIonInfo :: IonInfo
} deriving (Show, Eq, Ord) -- Maybe write my own Show instance.

makeClassy ''TentativelyAssignedFA

instance HasIonInfo TentativelyAssignedFA where
  ionInfo = tentativelyAssignedFAIonInfo

data TentativelyAssignedFAs = TentativelyAssignedFAs {
    _tentativelyAssignedFAsPrecIon :: Mz
  , _getTentativelyAssignedFAs :: [TentativelyAssignedFA]
} deriving (Eq, Ord)

makeClassy ''TentativelyAssignedFAs

instance Show TentativelyAssignedFAs where
  show (TentativelyAssignedFAs p fas) =
    "Precursor ion: " ++ show p ++ "\n" ++
    showList' fas

assignFA :: Double -> MonoisotopicMass -> [FattyAcyl] -> Maybe FattyAcyl
assignFA n m fas =
  case fas of
    [] -> Nothing
    f:fs -> if withinTolerance n m (faMonoisotopicMass f)
              then Just f
              else assignFA n m fs

neutralLossRowToFA :: NeutralLossRow -> TentativelyAssignedFA
neutralLossRowToFA (NeutralLossRow nl i) =
  TentativelyAssignedFA (assignFA 0.3 nl fattyAcyls) i

toTentativelyAssignedFAs :: NeutralLossSpectrum -> TentativelyAssignedFAs
toTentativelyAssignedFAs (NeutralLossSpectrum p nls) =
  TentativelyAssignedFAs p (neutralLossRowToFA <$> nls)

collectFAs :: TentativelyAssignedFAs -> [FattyAcyl]
collectFAs fas =
  catMaybes $ _getTentativelyAssignedFA <$> fas ^. getTentativelyAssignedFAs
