{-# LANGUAGE TemplateHaskell #-}
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
    mkMolecularFormula [(C, 3), (H, 6)] |+|
    foldMap toMolecularFormula (tg^..allFAs)

instance ToElementalComposition Triacylglycerol where
  toElementalComposition = toElementalComposition . toMolecularFormula


tagFAs :: Triacylglycerol -> [FattyAcyl]
tagFAs tg = nub $ tg^..allFAs

-- Brute force function to find possible TAGs. Write a more efficient
-- implementation if there are performance issues.
possibleTAGs :: Double -> MonoisotopicMass -> [FattyAcyl] -> [Triacylglycerol]
possibleTAGs n prec fas = nub $ do
  fa1 <- fas
  fa2 <- fas
  fa3 <- fas
  guard (withinTolerance n prec (I.monoisotopicMass (Triacylglycerol fa1 fa2 fa3)))
  return $ Triacylglycerol fa1 fa2 fa3

findPossibleTAGs :: AssignedFAs -> [Triacylglycerol]
findPossibleTAGs fas =
  possibleTAGs 0.3
               (fas ^. monoisotopicMass)
               (collectFAs fas)

data AssignedTAGs = AssignedTAGs {
    _tagSpectrumRow :: Maybe SpectrumRow
  , _tags           :: [Triacylglycerol]
  , _tagAssignedFAs :: AssignedFAs
} deriving (Show, Eq, Ord)

makeLenses ''AssignedTAGs

instance HasAssignedFAs AssignedTAGs where
  assignedFAs = tagAssignedFAs

instance HasMonoisotopicMass AssignedTAGs where
  monoisotopicMass = assignedFAs.monoisotopicMass

assignTAGs :: MSSpectrum -> AssignedFAs -> AssignedTAGs
assignTAGs spec fas =
  AssignedTAGs (findPrecursorIon 0.3 (fas ^. monoisotopicMass) spec)
               (findPossibleTAGs fas)
               fas

-- allTentativelyAssignedFAs :: [AssignedTAGs] -> [FattyAcyl]
-- allTentativelyAssignedFAs tgs =
--   sort . nub $
--     tgs^..traverse.tagAssignedFAs.getAssignedFAs.traverse.getAssignedFA._Just

allAssignedTAGs :: [AssignedTAGs] -> [Triacylglycerol]
allAssignedTAGs tgs = sort . nub . concat $ tgs^..traverse.tags

allTagFAs :: [AssignedTAGs] -> [FattyAcyl]
allTagFAs tgs = sort . nub $ allAssignedTAGs tgs >>= tagFAs

totalTagIntensity :: [AssignedTAGs] -> Intensity
totalTagIntensity = foldMap intensity'
  where
    intensity' (AssignedTAGs specRow _ _) =
      maybe (Intensity 0) (^.intensity) specRow

normalisedAbundanceFAsIndentified :: AssignedTAGs -> NormalisedAbundance
normalisedAbundanceFAsIndentified tags =
  sum normalisedAbundList
  where
    normalisedAbundList =
      _afNormalisedAbundance <$>
      filter (isJust . _getAssignedFA) (tags ^. getAssignedFAs)

formatNormalisedAbundanceFAsIndentified tags =
  tags^.monoisotopicMass.to showVal <>
  ", " <>
  (showVal . normalisedAbundanceFAsIndentified $ tags)

correctionRatio :: Intensity -> AssignedTAGs -> Double
correctionRatio (Intensity i) (AssignedTAGs r _ _) =
  maybe 0 (\x -> (x ^. getIntensity) / i) r

relativeAbundanceOfTags :: [AssignedTAGs] -> [Double]
relativeAbundanceOfTags tgs =
  correctionRatio (totalTagIntensity tgs) <$> tgs

collectAssignedTagsFas :: AssignedTAGs -> [(FattyAcyl, NormalisedAbundance)]
collectAssignedTagsFas tgs = [ (fa, na) | (Just fa, na) <- fattyAcylNormAbun]
  where
    assignedFAList = tgs^.tagAssignedFAs.getAssignedFAs
    fattyAcylNormAbun =
       zip (assignedFAList^..traverse.getAssignedFA)
           (assignedFAList^..traverse.normalisedAbundance)

data FinalResult = FinalResult {
    _finalResultMz                  :: Mz
  , _finalResultMzrelativeAbundance :: Double
  , _finalResultFAs                 :: [(FattyAcyl, NormalisedAbundance)]
  , _finalResultTags                :: [Triacylglycerol]
} deriving (Show, Eq, Ord)

makeLenses ''FinalResult

toFinalResult :: Intensity -> AssignedTAGs -> FinalResult
toFinalResult ti tgs =
  FinalResult (tgs^.monoisotopicMass)
              (correctionRatio ti tgs)
              (collectAssignedTagsFas tgs)
              (tgs^.tags)

finalResults :: [AssignedTAGs] -> [FinalResult]
finalResults tgs = toFinalResult (totalTagIntensity tgs) <$> tgs

tagMzNormalisedAbundances frs =
  (\x -> showVal (x^.finalResultMz) <> ": " <> take 4 (show (100 * x^.finalResultMzrelativeAbundance)))
    <$> frs

sumShouldEqual1 frs = sum $ frs^..traverse.finalResultMzrelativeAbundance

reCalNormalisedAbundance ::
  Double -> [(FattyAcyl, NormalisedAbundance)] -> [(FattyAcyl, NormalisedAbundance)]
reCalNormalisedAbundance n = mapped._2.getNormalisedAbundance %~ (*n)

accumulateNormalisedAbundance :: [FinalResult] -> [(FattyAcyl, NormalisedAbundance)]
accumulateNormalisedAbundance frs = mapped._2 %~ sum $ aggResult
  where
    aggResult = aggregateAL . concat $
      zipWith reCalNormalisedAbundance
              (frs^..traverse.finalResultMzrelativeAbundance)
              (frs^..traverse.finalResultFAs)

renderFattyAcylNormalisedAbundance :: (FattyAcyl, NormalisedAbundance) -> String
renderFattyAcylNormalisedAbundance (fa, na) = show fa <> ", " <> showVal na
