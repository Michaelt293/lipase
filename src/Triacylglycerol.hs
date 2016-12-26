{-# LANGUAGE TemplateHaskell #-}
module Triacylglycerol where

import FattyAcid
import Spectra
import Isotope hiding (monoisotopicMass)
import qualified Isotope as I
import Data.List
import Data.Maybe
import Control.Monad
import Control.Lens

data Triacylglycerol = Triacylglycerol {
    _fa1 :: FattyAcyl
  , _fa2 :: FattyAcyl
  , _fa3 :: FattyAcyl
}

makeClassy ''Triacylglycerol

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
  guard (withinTolerance n prec (I.monoisotopicMass (Triacylglycerol fa1 fa2 fa3)))
  return $ Triacylglycerol fa1 fa2 fa3

findPossibleTAGs :: AssignedFAs -> [Triacylglycerol]
findPossibleTAGs fas =
  possibleTAGs 0.3
               (fas ^. assignedFAsPrecIon)
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

assignTAGs :: MSSpectrum -> AssignedFAs -> AssignedTAGs
assignTAGs spec fas =
  AssignedTAGs (findPrecursorIon 0.3 (fas ^. assignedFAsPrecIon) spec)
               (findPossibleTAGs fas)
               fas

allTentativelyAssignedFAs :: [AssignedTAGs] -> [FattyAcyl]
allTentativelyAssignedFAs tgs =
  sort . nub . catMaybes $
    _getAssignedFA <$> concat (_getAssignedFAs . _tagAssignedFAs <$> tgs)

allAssignedTAGs :: [AssignedTAGs] -> [Triacylglycerol]
allAssignedTAGs tgs = sort . nub . concat $ _tags <$> tgs

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
      _ionInfoNormalisedAbundance . _assignedFAIonInfo <$>
      filter (isJust . _getAssignedFA) (fas ^. getAssignedFAs)

correctionRatio :: AssignedTAGs -> Intensity -> Double
correctionRatio (AssignedTAGs r _ _) (Intensity i) =
  case r of
    Nothing -> 0
    Just r' ->  (r' ^. getIntensity) / i
