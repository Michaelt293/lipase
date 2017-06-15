{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TemplateHaskell #-}
module FattyAcid where

import Control.Lens

import Spectra
import qualified Isotope as I
import Isotope hiding (monoisotopicMass)
import Data.Monoid ((<>))

newtype NumCarbons = NumCarbons { _getNumCarbons :: Int }
  deriving (Eq, Ord, Num)

makeClassy ''NumCarbons

instance Monoid NumCarbons where
  mempty = 0
  mappend = (+)

instance Show NumCarbons where
  show (NumCarbons n) = show n

newtype NumDoubleBonds = NumDoubleBonds { _getNumDoubleBonds :: Int }
  deriving (Eq, Ord, Num)

makeClassy ''NumDoubleBonds

instance Monoid NumDoubleBonds where
  mempty = 0
  mappend = (+)

instance Show NumDoubleBonds where
  show (NumDoubleBonds n) = show n

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
  show (FattyAcyl cs dbs) = show cs <> ":" <> show dbs

instance ToElementalComposition FattyAcyl where
  toElementalComposition (FattyAcyl cs dbs) = mkElementalComposition
    [ (C, cs^.getNumCarbons)
    , (H, cs^.getNumCarbons.to (\x -> x * 2 - 1) -
          dbs^.getNumDoubleBonds.to (*2))
    , (O, 2)
    ]
  charge _ = Just 0

newtype FattyAcid = FattyAcid { _getFattyAcyl :: FattyAcyl }
  deriving (Eq, Ord)

makeLenses ''FattyAcid

instance Show FattyAcid where
  show (FattyAcid a) = "FA " <> show a

instance ToElementalComposition FattyAcid where
  toElementalComposition (FattyAcid fa) =
    toElementalComposition fa <> mkElementalComposition [(H, 1)]
  charge _ = Just 0


c10_0 :: FattyAcid
c10_0 = FattyAcid (FattyAcyl 10 0)

c12_0 :: FattyAcid
c12_0 = FattyAcid (FattyAcyl 12 0)

c14s :: [FattyAcid]
c14s = FattyAcid . FattyAcyl 14 <$> [0, 1]

c15s :: [FattyAcid]
c15s = FattyAcid . FattyAcyl 15 <$> [0, 1]

c16s :: [FattyAcid]
c16s = FattyAcid . FattyAcyl 16 <$> [0, 1, 2]

c17s :: [FattyAcid]
c17s = FattyAcid . FattyAcyl 17 <$> [0, 1]

c18s :: [FattyAcid]
c18s = FattyAcid . FattyAcyl 18 <$> [0, 1, 2, 3, 4]

c19s :: [FattyAcid]
c19s = FattyAcid . FattyAcyl 19 <$> [0, 1]

c20s :: [FattyAcid]
c20s = FattyAcid . FattyAcyl 20 <$> [0, 1, 2, 3, 4, 5]

c21s :: [FattyAcid]
c21s = FattyAcid . FattyAcyl 21 <$> [0, 1]

c22s :: [FattyAcid]
c22s = FattyAcid . FattyAcyl 22 <$> [0, 1, 2, 3, 4, 5, 6]

evenChainFAs :: [FattyAcid]
evenChainFAs = [c10_0, c12_0] <> c14s <> c16s <> c18s <> c20s <> c22s

oddChainFAs :: [FattyAcid]
oddChainFAs = c15s <> c17s <> c19s <> c21s

fattyAcids :: [FattyAcid]
fattyAcids = evenChainFAs <> oddChainFAs

fattyAcidMonoisotopicMasses :: [MonoisotopicMass]
fattyAcidMonoisotopicMasses = I.monoisotopicMass <$> fattyAcids

assignFAsFromNeutralLoss
  :: MS2Spectrum a MonoisotopicMass -> MS2Spectrum a (Maybe FattyAcid)
assignFAsFromNeutralLoss =
  fmap (\nl -> assignFA 0.4 nl fattyAcids)

assignFA
  :: MonoisotopicMass -> MonoisotopicMass -> [FattyAcid] -> Maybe FattyAcid
assignFA n m fas =
  case fas of
    [] -> Nothing
    fa:fas' -> if withinTolerance n m (I.monoisotopicMass fa)
                 then Just fa
                 else assignFA n m fas'
