{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleInstances #-}
module Triacylglycerol where

import FattyAcid
import Spectra
import Isotope hiding (monoisotopicMass)
import qualified Isotope as I
import Data.List
import Data.List.HIUtils
import Data.Maybe
import Data.Monoid
import Data.Ord
import Data.Map (Map)
import qualified Data.Map as M
import Control.Monad
import Control.Lens
import Control.Arrow

data Triacylglycerol = Triacylglycerol {
    _fa1 :: FattyAcyl
  , _fa2 :: FattyAcyl
  , _fa3 :: FattyAcyl
}

makeClassy ''Triacylglycerol

allFAs :: Applicative f =>
  (FattyAcyl -> f FattyAcyl) -> Triacylglycerol -> f Triacylglycerol
allFAs f (Triacylglycerol fa1 fa2 fa3) =
  Triacylglycerol <$> f fa1 <*> f fa2 <*> f fa3

instance Eq Triacylglycerol where
  tg1 == tg2 = sort (tg1^..allFAs) == sort (tg2^..allFAs)

instance Ord Triacylglycerol where
  tg1 <= tg2 = sort (tg1^..allFAs) <= sort (tg2^..allFAs)

instance Show Triacylglycerol where
  show tg =
    "TG(" <> intercalate "_" (tg^..allFAs.to show) <> ")"

instance ToMolecularFormula Triacylglycerol where
  toMolecularFormula tg =
    mkMolecularFormula [(C, 3), (H, 5)] <>
    foldMap toMolecularFormula (tg^..allFAs)

instance ToElementalComposition Triacylglycerol where
  toElementalComposition = toElementalComposition . toMolecularFormula

newtype Protonated a = Protonated { _getProtonatedIon :: a }
  deriving (Eq, Ord)

instance ToMolecularFormula a => ToMolecularFormula (Protonated a) where
  toMolecularFormula (Protonated p) =
    toMolecularFormula p |+| mkMolecularFormula [(H, 1)]

instance ToElementalComposition a => ToElementalComposition (Protonated a) where
  toElementalComposition (Protonated p) =
    toElementalComposition p |+| mkElementalComposition [(H, 1)]

instance ToElementalComposition a => Ion (Protonated a) where
  charge _ = 1

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
  show (CondensedTriacylglycerol cs dbs) = "TAG " <> showVal cs <> ":" <> showVal dbs

tagToCondensedTriacylglycerol tg =
  CondensedTriacylglycerol
    (foldOf (folded.numCarbons) (tg^..allFAs))
    (foldOf (folded.numDoubleBonds) (tg^..allFAs))

tagsToCondensedTags :: [Triacylglycerol] -> [CondensedTriacylglycerol]
tagsToCondensedTags tags = nub $ tagToCondensedTriacylglycerol <$> tags

--data AssignedTAGs a = AssignedTAGs {
--    _tagSpectrumRow :: Maybe (SpectrumRow a)
--  , _tags           :: [Triacylglycerol]
--  , _tagAssignedFAs :: AssignedFAs
--} deriving (Show, Eq, Ord)
--
--makeLenses ''AssignedTAGs
--
--instance HasAssignedFAs (AssignedTAGs a) where
--  assignedFAs = tagAssignedFAs
--
--instance HasMonoisotopicMass (AssignedTAGs a) where
--  monoisotopicMass = assignedFAs.monoisotopicMass
--
--assignTAGs :: MSSpectrum Mz -> AssignedFAs -> AssignedTAGs Mz
--assignTAGs spec fas =
--  AssignedTAGs (findPrecursorIon 0.3 (fas ^. monoisotopicMass) spec)
--               (findPossibleTAGs fas)
--               fas
--
---- allTentativelyAssignedFAs :: [AssignedTAGs] -> [FattyAcyl]
---- allTentativelyAssignedFAs tgs =
----   sort . nub $
----     tgs^..traverse.tagAssignedFAs.getAssignedFAs.traverse.getAssignedFA._Just
--

--
--allTagFAs :: [AssignedTAGs a] -> [FattyAcyl]
--allTagFAs tgs = sort . nub $ allAssignedTAGs tgs >>= tagFAs
--
-- totalTagIntensity :: [MS2Spectrum (Mz, Intensity, [Triacylglycerol]) a] -> Intensity
-- totalTagIntensity =  foldOf (folded.precursorIon._2)
--    intensity' (AssignedTAGs specRow _ _) =
--      maybe (Intensity 0) (^.intensity) specRow
--
--normalisedAbundanceFAsIndentified :: AssignedTAGs a -> NormalisedAbundance
--normalisedAbundanceFAsIndentified tags =
--  sum normalisedAbundList
--  where
--    normalisedAbundList =
--      _afNormalisedAbundance <$>
--      filter (isJust . _getAssignedFA) (tags ^. getAssignedFAs)
--
--formatNormalisedAbundanceFAsIndentified tags =
--  tags^.monoisotopicMass.to showVal <>
--  ", " <>
--  (showVal . normalisedAbundanceFAsIndentified $ tags)
--

totalTagIntensity :: [MS2Spectrum (Mz, Intensity, [Triacylglycerol]) a] -> Intensity
totalTagIntensity = foldOf (folded.precursorIon._2)

correctionRatio :: Intensity -> Intensity -> Double
correctionRatio i total = (i / total)^.getIntensity

---- maybe 0 (\x -> (x ^. getIntensity) / i) r
--
--relativeAbundanceOfTags :: [AssignedTAGs a] -> [Double]
--relativeAbundanceOfTags tgs =
--  correctionRatio (totalTagIntensity tgs) <$> tgs
--
--collectAssignedTagsFas :: AssignedTAGs a -> [(FattyAcyl, NormalisedAbundance)]
--collectAssignedTagsFas tgs = [ (fa, na) | (Just fa, na) <- fattyAcylNormAbun]
--  where
--    assignedFAList = tgs^.tagAssignedFAs.getAssignedFAs
--    fattyAcylNormAbun =
--       zip (assignedFAList^..traverse.getAssignedFA) -- _Just
--           (assignedFAList^..traverse.normalisedAbundance)
--
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
--
--tagMzNormalisedAbundances frs =
--  (\x -> show (tagToCondensedTriacylglycerol <$> (x^.finalResultTags)) <> ": " <> take 4 (show (100 * x^.finalResultMzrelativeAbundance)))
--    <$> frs
--
--sumShouldEqual1 frs = sum $ frs^..traverse.finalResultMzrelativeAbundance
--
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
    --(\x -> reCalNormalisedAbundance (x^.finalResultMzNormalisedAbundance.getNormalisedAbundance))
    --(frs^..traverse.finalResultFAs)

-- mapped._2 %~ sum $ aggResult
--  where
--    aggResult = aggregateAL . concat $
--      zipWith reCalNormalisedAbundance
--              (frs^..traverse.finalResultMzNormalisedAbundance)
--              (frs^..traverse.finalResultFAs)
--
renderPairNormalisedAbundance :: (Show a) => (a, NormalisedAbundance) -> String
renderPairNormalisedAbundance (a, na) = show a <> ", " <> showVal na

identifiedTags :: [FinalResult] -> [Triacylglycerol]
identifiedTags frs = sort . nub $ frs^..traverse.finalResultTags.traverse

allTagFAs :: [FinalResult] -> [FattyAcid]
allTagFAs frs = sort . nub . concat $ fmap M.keys (frs^..traverse.finalResultFAs)

identifiedCondensedTags :: [FinalResult] -> [([CondensedTriacylglycerol], NormalisedAbundance)]
identifiedCondensedTags frs =
  zip (fmap tagsToCondensedTags (frs^..traverse.finalResultTags))
      (frs^..traverse.normalisedAbundance)
