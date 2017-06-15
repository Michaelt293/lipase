{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleInstances #-}
module Triacylglycerol where

import Control.Lens

import FattyAcid
import Spectra
import Isotope hiding (monoisotopicMass)
import Isotope.Ion
import Data.List (sort, intercalate, nub)
import Data.Maybe (fromJust, isJust)
import Data.Monoid ((<>))
import Data.Map (Map)
import qualified Data.Map as M
import Control.Monad (guard)


data Triacylglycerol = Triacylglycerol {
    _fa1 :: FattyAcyl
  , _fa2 :: FattyAcyl
  , _fa3 :: FattyAcyl
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
    "TG(" <> intercalate "_" (tg^..allFAs.to show) <> ")"

instance ToElementalComposition Triacylglycerol where
  toElementalComposition tg =
    mkElementalComposition [(C, 3), (H, 5)] <>
    foldMap toElementalComposition (tg^..allFAs)
  charge _ = Just 0

makeLenses ''Protonated

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
  :: MSSpectrum Mz -> MS2Spectrum Mz (Maybe FattyAcid) -> MS2Spectrum (Mz, Intensity, [Triacylglycerol]) (Maybe FattyAcid)
findPossibleTAGs msSpec ms2Spec =
  precursorIon .~ (precursorIon', intensity', tags) $ ms2Spec
    where
      precursorIon' = ms2Spec^.precursorIon
      tags = possibleTAGs 0.4 precursorIon'
               (ms2Spec^..ms2Spectrum.traverse.ion._Just)
      intensity' = findPrecursorIon 0.4 precursorIon' msSpec^._Just.intensity

data CondensedTriacylglycerol = CondensedTriacylglycerol {
    _totalFACarbons        :: NumCarbons
  , _totalFANumDoubleBonds :: NumDoubleBonds
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

data FinalResult = FinalResult {
    _finalResultMz                    :: Mz
  , _finalResultMzNormalisedAbundance :: NormalisedAbundance
  , _finalResultFAs                   :: Map FattyAcid NormalisedAbundance
  , _finalResultTags                  :: [Triacylglycerol]
} deriving (Show, Eq, Ord)

makeLenses ''FinalResult

instance HasNormalisedAbundance FinalResult where
  normalisedAbundance = finalResultMzNormalisedAbundance
--
toFinalResult :: Intensity -> MS2Spectrum (Mz, Intensity, [Triacylglycerol]) (Maybe FattyAcid) -> FinalResult
toFinalResult totalIntensity spec =
  FinalResult (spec^.precursorIon._1)
              (NormalisedAbundance correctionRatio')
              (M.fromList (fmap (\x -> (x^.ion.to fromJust, x^.normalisedAbundance)) spectrumRows))
              (spec^.precursorIon._3)
  where
    spectrumRows :: [SpectrumRow (Maybe FattyAcid)]
    spectrumRows = spec^..ms2Spectrum.folded.filtered (\x -> x^.ion.to isJust)
    correctionRatio' :: Double
    correctionRatio' = correctionRatio (spec^.precursorIon._2) totalIntensity

--
toFinalResults :: [MS2Spectrum (Mz, Intensity, [Triacylglycerol]) (Maybe FattyAcid)] -> [FinalResult]
toFinalResults specs = toFinalResult totalIntensity <$> specs
  where
    totalIntensity = foldOf (folded.precursorIon._2) specs

totalNormalisedAbundance :: Map FattyAcid NormalisedAbundance -> NormalisedAbundance
totalNormalisedAbundance fas = foldMap (^._2) $ M.toList fas

reCalNormalisedAbundance ::
  Double -> Map FattyAcid NormalisedAbundance -> Map FattyAcid NormalisedAbundance
reCalNormalisedAbundance n fas = M.fromList $
  mapped._2.getNormalisedAbundance %~ (*n) $ M.toList fas

identifiedFAs :: [FinalResult] -> [(FattyAcid, NormalisedAbundance)]
identifiedFAs frs =
  sort . M.toList $ foldr
    (M.unionWith mappend)
    M.empty
    reCalNormalisedAbundances
    where
      correctionRatio' fr = fr^.finalResultMzNormalisedAbundance.getNormalisedAbundance
      reCalNormalisedAbundance' fr = reCalNormalisedAbundance (correctionRatio' fr) (fr^.finalResultFAs)
      reCalNormalisedAbundances = fmap reCalNormalisedAbundance' frs

renderPairNormalisedAbundance :: (Show a) => (a, NormalisedAbundance) -> String
renderPairNormalisedAbundance (a, na) = show a <> ", " <> show na

identifiedTags :: [FinalResult] -> [Triacylglycerol]
identifiedTags frs = sort . nub $ frs^..traverse.finalResultTags.traverse

identifiedTAGSummary :: [FinalResult] -> [([CondensedTriacylglycerol], [Triacylglycerol])]
identifiedTAGSummary frs = filter (/= ([],[])) $
  zip (fmap tagsToCondensedTags finalTags) finalTags
      where
        finalTags = frs^..traverse.finalResultTags

allTagFAs :: [FinalResult] -> [FattyAcid]
allTagFAs frs = sort . nub . concat $ fmap M.keys (frs^..traverse.finalResultFAs)

identifiedCondensedTags :: [FinalResult] -> [([CondensedTriacylglycerol], NormalisedAbundance)]
identifiedCondensedTags frs =
  zip (fmap tagsToCondensedTags (frs^..traverse.finalResultTags))
      (frs^..traverse.normalisedAbundance)
