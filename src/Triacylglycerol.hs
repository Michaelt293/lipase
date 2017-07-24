{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
module Triacylglycerol where

import Control.Lens
import FattyAcid
import Spectra
import Isotope hiding (monoisotopicMass)
import Isotope.Ion
import Data.List (sort, intercalate, nub)
import Data.Monoid ((<>))
import Data.Csv (ToField(..))
import Data.ByteString.Char8 (pack)
import Control.Monad (guard)


data Triacylglycerol = Triacylglycerol {
    _fa1 :: !FattyAcyl
  , _fa2 :: !FattyAcyl
  , _fa3 :: !FattyAcyl
}

makeClassy ''Triacylglycerol

allFAs :: Applicative f =>
  (FattyAcyl -> f FattyAcyl) -> Triacylglycerol -> f Triacylglycerol
allFAs f (Triacylglycerol fa1' fa2' fa3') =
  Triacylglycerol <$> f fa1' <*> f fa2' <*> f fa3'

instance Eq Triacylglycerol where
  tg1 == tg2 = sort (tg1^..allFAs) == sort (tg2^..allFAs)

instance Show Triacylglycerol where
  show tg =
    "TG " <> intercalate "_" (tg^..allFAs.to show)

instance ToElementalComposition Triacylglycerol where
  toElementalComposition tg =
    mkElementalComposition [(C, 3), (H, 5)] <>
    foldMap toElementalComposition (tg^..allFAs)
  charge _ = Just 0

instance ToField [CondensedTriacylglycerol] where
  toField tags =
    if null tags
      then "no triacylglycerols identified"
      else pack . intercalate ", " $ show <$> tags

instance ToField [Triacylglycerol] where
  toField tags =
    if null tags
      then "no triacylglycerols identified"
      else pack . intercalate ", " $ show <$> tags

instance IsSaturated Triacylglycerol where
  isSaturated tag = foldOf (allFAs.numDoubleBonds) tag == 0

instance IsMonounsaturated Triacylglycerol where
  isMonounsaturated tag = foldOf (allFAs.numDoubleBonds) tag == 1

instance IsPolyunsaturated Triacylglycerol where
  isPolyunsaturated tag = foldOf (allFAs.numDoubleBonds) tag > 1

tagFAs :: Triacylglycerol -> [FattyAcyl]
tagFAs tg = nub $ tg^..allFAs

-- Brute force function to find possible TAGs. Write a more efficient
-- implementation if there are performance issues.
possibleTAGs :: Mz -> Mz -> [FattyAcid] -> [Triacylglycerol]
possibleTAGs n prec fas = nub $ do
  FattyAcid _fa1 <- fas
  FattyAcid _fa2 <- fas
  FattyAcid _fa3 <- fas
  guard . withinTolerance n prec .
    mz $ Protonated Triacylglycerol {..}
  return Triacylglycerol {..}

findPossibleTAGs
  :: MSSpectrum Mz
  -> MS2Spectrum Mz (Maybe FattyAcid)
  -> MS2Spectrum (Mz, Intensity, [Triacylglycerol]) (Maybe FattyAcid)
findPossibleTAGs msSpec ms2Spec =
  precursorIon .~ (precursorIon', intensity', tags) $ ms2Spec
    where
      precursorIon' = ms2Spec^.precursorIon
      tags = possibleTAGs 0.4 precursorIon'
               (ms2Spec^..ms2Spectrum.traverse.ion._Just)
      intensity' = findPrecursorIon 0.4 precursorIon' msSpec^._Just.intensity

data CondensedTriacylglycerol = CondensedTriacylglycerol {
    _totalFACarbons        :: !NumCarbons
  , _totalFANumDoubleBonds :: !NumDoubleBonds
} deriving (Eq, Ord)

makeClassy ''CondensedTriacylglycerol

instance Show CondensedTriacylglycerol where
  show (CondensedTriacylglycerol cs dbs) = "TG " <> show cs <> ":" <> show dbs

tagToCondensedTriacylglycerol :: Triacylglycerol -> CondensedTriacylglycerol
tagToCondensedTriacylglycerol tg =
  CondensedTriacylglycerol
    (foldOf (folded.numCarbons) (tg^..allFAs))
    (foldOf (folded.numDoubleBonds) (tg^..allFAs))

instance Ord Triacylglycerol where
  tg1 `compare` tg2 =
    case tagToCondensedTriacylglycerol tg1 `compare` tagToCondensedTriacylglycerol tg2 of
      LT -> LT
      EQ -> sort (tg1^..allFAs) `compare` sort (tg2^..allFAs)
      GT -> GT

tagsToCondensedTags :: [Triacylglycerol] -> [CondensedTriacylglycerol]
tagsToCondensedTags tags = nub $ tagToCondensedTriacylglycerol <$> tags

totalTagIntensity :: [MS2Spectrum (Mz, Intensity, [Triacylglycerol]) a] -> Intensity
totalTagIntensity = foldOf (folded.precursorIon._2)

correctionRatio :: Intensity -> Intensity -> Double
correctionRatio i total = (i / total)^.getIntensity
