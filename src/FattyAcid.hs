{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TemplateHaskell #-}
module FattyAcid where

import qualified Isotope as I
import Isotope hiding (monoisotopicMass)
import Spectra
import Data.Maybe
import Data.Monoid
import Control.Lens

newtype NumCarbons = NumCarbons { _getNumCarbons :: Int }
  deriving (Show, Eq, Ord, Num)

makeClassy ''NumCarbons

instance Monoid NumCarbons where
  mempty = 0
  mappend = (+)

instance ShowVal NumCarbons where
  showVal (NumCarbons n) = show n

newtype NumDoubleBonds = NumDoubleBonds { _getNumDoubleBonds :: Int }
  deriving (Show, Eq, Ord, Num)

makeClassy ''NumDoubleBonds

instance Monoid NumDoubleBonds where
  mempty = 0
  mappend = (+)

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
  show (FattyAcyl cs dbs) = showVal cs <> ":" <> showVal dbs

instance ToMolecularFormula FattyAcyl where
  toMolecularFormula (FattyAcyl cs dbs) = mkMolecularFormula
    [ (C, cs^.getNumCarbons)
    , (H, cs^.getNumCarbons.to (\x -> x * 2 - 1) - dbs^.getNumDoubleBonds.to (*2))
    , (O, 2)
    ]

instance ToElementalComposition FattyAcyl where
  toElementalComposition = toElementalComposition . toMolecularFormula

newtype FattyAcid = FattyAcid { _getFattyAcyl :: FattyAcyl }
  deriving (Eq, Ord)

makeLenses ''FattyAcid

instance Show FattyAcid where
  show (FattyAcid a) = show a

instance ToMolecularFormula FattyAcid where
  toMolecularFormula (FattyAcid fa) =
    toMolecularFormula fa |+| mkMolecularFormula [(H, 1)]

instance ToElementalComposition FattyAcid where
  toElementalComposition = toElementalComposition . toMolecularFormula

c10_0 = FattyAcid (FattyAcyl 10 0)

c12_0 = FattyAcid (FattyAcyl 12 0)

c14s = FattyAcid . FattyAcyl 14 <$> [0, 1]

c15s = FattyAcid . FattyAcyl 15 <$> [0, 1]

c16s = FattyAcid . FattyAcyl 16 <$> [0, 1, 2]

c17s = FattyAcid . FattyAcyl 17 <$> [0, 1]

c18s = FattyAcid . FattyAcyl 18 <$> [0, 1, 2, 3, 4]

c19s = FattyAcid . FattyAcyl 19 <$> [0, 1]

c20s = FattyAcid . FattyAcyl 20 <$> [0, 1, 2, 3, 4, 5]

c21s = FattyAcid . FattyAcyl 21 <$> [0, 1]

c22s = FattyAcid . FattyAcyl 22 <$> [0, 1, 2, 3, 4, 5, 6]

evenChainFAs :: [FattyAcid]
evenChainFAs = [c10_0, c12_0] <> c14s <> c16s <> c18s <> c20s <> c22s

oddChainFAs :: [FattyAcid]
oddChainFAs = c15s <> c17s <> c19s <> c21s

fattyAcids :: [FattyAcid]
fattyAcids = evenChainFAs <> oddChainFAs

-- faMonoisotopicMass :: FattyAcyl -> MonoisotopicMass
-- faMonoisotopicMass fa = I.monoisotopicMass fa |+| I.monoisotopicMass H

fattyAcidMonoisotopicMasses :: [MonoisotopicMass]
fattyAcidMonoisotopicMasses = I.monoisotopicMass <$> fattyAcids

-- data AssignedFA = AssignedFA {
--     _getAssignedFA         :: Maybe FattyAcyl
--   , _afIntensity           :: Intensity
--   , _afRelativeAbundance   :: RelativeAbundance
--   , _afNormalisedAbundance :: NormalisedAbundance
-- } deriving (Show, Eq, Ord) -- Maybe write my own Show instance.
--
-- makeClassy ''AssignedFA
--
-- instance HasIntensity AssignedFA where
--   intensity = afIntensity
--
-- instance HasRelativeAbundance AssignedFA where
--   relativeAbundance = afRelativeAbundance
--
-- instance HasNormalisedAbundance AssignedFA where
--   normalisedAbundance = afNormalisedAbundance
--
-- data AssignedFAs = AssignedFAs {
--     _assignedFAsPrecIon :: Mz
--   , _getAssignedFAs :: [AssignedFA]
-- } deriving (Eq, Ord)
--
-- makeClassy ''AssignedFAs
--
-- instance HasMonoisotopicMass AssignedFAs where
--   monoisotopicMass = assignedFAsPrecIon
--
-- instance Show AssignedFAs where
--   show (AssignedFAs p fas) =
--     "Precursor ion: " <> show p <> "\n" <>
--     showList' fas

assignFAsFromNeutralLoss :: MS2Spectrum a MonoisotopicMass -> MS2Spectrum a (Maybe FattyAcid)
assignFAsFromNeutralLoss =
  fmap (\nl -> assignFA 0.4 nl fattyAcids)

assignFA :: MonoisotopicMass -> MonoisotopicMass -> [FattyAcid] -> Maybe FattyAcid
assignFA n m fas =
  case fas of
    [] -> Nothing
    fa:fas -> if withinTolerance n m (I.monoisotopicMass fa)
              then Just fa
              else assignFA n m fas

-- toAssignedFAs :: NeutralLossSpectrum -> AssignedFAs
-- toAssignedFAs (NeutralLossSpectrum p nls) =
--   AssignedFAs p (neutralLossRowToFA <$> nls)

-- collectFAs :: AssignedFAs -> [FattyAcyl]
-- collectFAs fas =
--   catMaybes $ fas ^.. getAssignedFAs.traverse.getAssignedFA
