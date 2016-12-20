{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module FattyAcid where

import Isotope
import Spectra

newtype NumCarbons = NumCarbons { getNumCarbons :: Int }
  deriving (Eq, Ord, Num)

instance Show NumCarbons where
  show (NumCarbons n) = show n

newtype NumDoubleBonds = NumDoubleBonds { getNumDoubleBonds :: Int }
  deriving (Eq, Ord, Num)

instance Show NumDoubleBonds where
  show (NumDoubleBonds n) = show n

data FattyAcyl = FattyAcyl {
    numCarbons     :: NumCarbons
  , numDoubleBonds :: NumDoubleBonds
} deriving (Eq, Ord)

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

fattyAcidMonoisotopicMasses :: [MonoisotopicMass]
fattyAcidMonoisotopicMasses =
  (\fa -> monoisotopicMass fa |+| monoisotopicMass H) <$> fattyAcyls

data TentativelyAssignedFA = TentativelyAssignedFA {
    tentativelyAssignedFA :: Maybe FattyAcyl
  , tentativelyAssignedFAIonInfo :: IonInfo
} deriving (Show, Eq, Ord) -- Maybe write my own Show instance.

type TentativelyAssignedFAs = [TentativelyAssignedFA]
