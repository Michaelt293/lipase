
module Triacylglycerol where

import FattyAcid
import Spectra
import Isotope
import Data.List
import Data.Maybe
import Control.Monad

data Triacylglycerol = Triacylglycerol {
    fa1 :: FattyAcyl
  , fa2 :: FattyAcyl
  , fa3 :: FattyAcyl
}

instance Eq Triacylglycerol where
  Triacylglycerol a1 b1 c1 == Triacylglycerol a2 b2 c2 =
    sort [a1, b1, c1] == sort [a2, b2, c2]

instance Ord Triacylglycerol where
  Triacylglycerol a1 b1 c1 <= Triacylglycerol a2 b2 c2 =
    sort [a1, b1, c1] <= sort [a2, b2, c2]

instance Show Triacylglycerol where
  show (Triacylglycerol a b c) =
    "TG(" ++ intercalate "_" (show <$> sort [a, b, c]) ++ ")"

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
  guard (withinTolerance n prec (monoisotopicMass (Triacylglycerol fa1 fa2 fa3)))
  return $ Triacylglycerol fa1 fa2 fa3

findPossibleTAGs :: TentativelyAssignedFAs -> [Triacylglycerol]
findPossibleTAGs fas =
  possibleTAGs 0.3
               (monoisotopicMass'
               (tentativelyAssignedFAsPrecIon fas))
               (collectFAs fas)

tagFAs :: Triacylglycerol -> [FattyAcyl]
tagFAs (Triacylglycerol fa1 fa2 fa3) = nub [fa1, fa2, fa3]

data AssignedTAGs = AssignedTAGs {
    spectrumRow :: Maybe SpectrumRow
  , tags :: [Triacylglycerol]
  , assignedFAs :: TentativelyAssignedFAs
} deriving (Show, Eq, Ord)

assignTAGs :: MSSpectrum -> TentativelyAssignedFAs -> AssignedTAGs
assignTAGs spec fas =
  AssignedTAGs (findPrecursorIon 0.3 (tentativelyAssignedFAsPrecIon fas) spec)
               (findPossibleTAGs fas)
               fas

allTentativelyAssignedFAs :: [AssignedTAGs] -> [FattyAcyl]
allTentativelyAssignedFAs tgs =
  sort . nub . catMaybes $
    tentativelyAssignedFA <$> concat (tentativelyAssignedFAs . assignedFAs <$> tgs)

allAssignedTAGs :: [AssignedTAGs] -> [Triacylglycerol]
allAssignedTAGs tgs = sort . nub . concat $ tags <$> tgs

allTagFAs :: [AssignedTAGs] -> [FattyAcyl]
allTagFAs tgs = sort . nub . concat $ tagFAs <$> allAssignedTAGs tgs

totalTagIntensity :: [AssignedTAGs] -> Intensity
totalTagIntensity = foldMap intensity'
  where
    intensity' (AssignedTAGs (Just (SpectrumRow _ (IonInfo i _ _))) _ _) = i
    intensity' (AssignedTAGs Nothing _ _) = Intensity 0

normalisedAbundanceFAsIndentified :: AssignedTAGs -> NormalisedAbundance
normalisedAbundanceFAsIndentified (AssignedTAGs _ _ fas) =
  if null normalisedAbundList
    then 0
    else sum normalisedAbundList
  where
    normalisedAbundList =
      normalisedAbundance . tentativelyAssignedFAIonInfo <$>
      filter (isJust . tentativelyAssignedFA) (tentativelyAssignedFAs fas)

correctionRatio :: AssignedTAGs -> Intensity -> Double
correctionRatio (AssignedTAGs r _ _) (Intensity i) =
  case r of
    Nothing -> 0
    Just r' -> (getIntensity . intensity . spectrumRowIonInfo) r' / i
