{-# LANGUAGE TemplateHaskell #-}
module Triacylglycerol where

import FattyAcid
import Spectra
import Isotope hiding (monoisotopicMass)
import qualified Isotope as I
import Data.List
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
  Triacylglycerol a1 b1 c1 == Triacylglycerol a2 b2 c2 =
    sort [a1, b1, c1] == sort [a2, b2, c2]

instance Ord Triacylglycerol where
  Triacylglycerol a1 b1 c1 <= Triacylglycerol a2 b2 c2 =
    sort [a1, b1, c1] <= sort [a2, b2, c2]

instance Show Triacylglycerol where
  show (Triacylglycerol a b c) =
    "TG(" <> intercalate "_" (show <$> sort [a, b, c]) <> ")"

instance ToMolecularFormula Triacylglycerol where
  toMolecularFormula (Triacylglycerol a b c) =
    mkMolecularFormula [(C, 3), (H, 6)] |+| foldMap toMolecularFormula [a, b, c]

instance ToElementalComposition Triacylglycerol where
  toElementalComposition = toElementalComposition . toMolecularFormula

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

tagFAs :: Triacylglycerol -> [FattyAcyl]
tagFAs (Triacylglycerol fa1 fa2 fa3) = nub [fa1, fa2, fa3]

data AssignedTAGs = AssignedTAGs {
    _tagSpectrumRow :: Maybe SpectrumRow
  , _tags :: [Triacylglycerol]
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

allTentativelyAssignedFAs :: [AssignedTAGs] -> [FattyAcyl]
allTentativelyAssignedFAs tgs =
  sort . nub . catMaybes $
    tgs^..traverse.tagAssignedFAs.getAssignedFAs.traverse.getAssignedFA

allAssignedTAGs :: [AssignedTAGs] -> [Triacylglycerol]
allAssignedTAGs tgs = sort . nub . concat $ _tags <$> tgs

allTagFAs :: [AssignedTAGs] -> [FattyAcyl]
allTagFAs tgs = sort . nub . concat $ tagFAs <$> allAssignedTAGs tgs

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
      _ionInfoNormalisedAbundance . _assignedFAIonInfo <$>
      filter (isJust . _getAssignedFA) (tags ^. getAssignedFAs)

formatNormalisedAbundanceFAsIndentified tags =
  tags^.monoisotopicMass.to showVal <> ", " <> (showVal . normalisedAbundanceFAsIndentified $ tags)

correctionRatio :: Intensity -> AssignedTAGs -> Double
correctionRatio (Intensity i) (AssignedTAGs r _ _) =
  maybe 0 (\x -> (x ^. getIntensity) / i) r

relativeAbundanceOfTags :: [AssignedTAGs] -> [Double]
relativeAbundanceOfTags tgs =
  correctionRatio (totalTagIntensity tgs) <$> tgs

--collectAssignedTagsFas :: AssignedTAGs -> [(FattyAcyl, NormalisedAbundance)]
collectAssignedTagsFas tgs = first fromJust <$> filter (\(fa, _) -> isJust fa)
  (zip (assignedFAList^..traverse.getAssignedFA)
       (assignedFAList^..traverse.assignedFAIonInfo.ionInfoNormalisedAbundance))
  where assignedFAList = tgs^.tagAssignedFAs.getAssignedFAs

data FinalResult = FinalResult {
    _finalResultMz :: Mz
  , _finalResultMzrelativeAbundance :: Double
  , _finalResultFAs :: [(FattyAcyl, NormalisedAbundance)]
  , _finalResultTags :: [Triacylglycerol]
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

-- | Sort, then group
-- from Data.List.HIUtils
aggregateBy :: (a -> a -> Ordering) -> [a] -> [[a]]
aggregateBy x = groupBy (\a b -> x a b == EQ) . sortBy x

-- | Aggregate an association list, such that keys become unique.
-- from Data.List.HIUtils
aggregateAL :: (Ord a) => [(a,b)] -> [(a,[b])]
aggregateAL = fmap (fst . head &&& fmap snd) . aggregateBy (comparing fst)

reCalNormalisedAbundance :: Double -> [(FattyAcyl, NormalisedAbundance)] -> [(FattyAcyl, NormalisedAbundance)]
reCalNormalisedAbundance n = mapped._2.getNormalisedAbundance %~ (*n)

accumulateNormalisedAbundance :: [FinalResult] -> [(FattyAcyl, NormalisedAbundance)]
accumulateNormalisedAbundance frs = mapped._2 %~ sum $ aggResult
  where
    aggResult = aggregateAL . concat $
      zipWith reCalNormalisedAbundance (frs^..traverse.finalResultMzrelativeAbundance) (frs^..traverse.finalResultFAs)

renderFattyAcylNormalisedAbundance :: (FattyAcyl, NormalisedAbundance) -> String
renderFattyAcylNormalisedAbundance (fa, na) = show fa <> ", " <> showVal na 
