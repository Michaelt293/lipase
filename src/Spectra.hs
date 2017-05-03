{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TemplateHaskell #-}

module Spectra where

import Control.Lens

import Isotope ( MonoisotopicMass
               , ToElementalComposition(..)
               , MonoisotopicMass(..)
               )
import Isotope.Ion (Mz(..))
import Data.List (intercalate)
import Data.Foldable (maximumBy)
import Data.Ord (comparing)
import Data.Monoid ((<>))
import Data.Csv (FromNamedRecord(..), decodeByName, (.:))
import Data.Vector (toList)
import qualified Data.ByteString.Lazy

makeClassy ''MonoisotopicMass

deriving instance Num MonoisotopicMass

deriving instance Fractional MonoisotopicMass

deriving instance Num Mz

deriving instance Fractional Mz

showMonoisotopicMass :: MonoisotopicMass -> String
showMonoisotopicMass (MonoisotopicMass v) = show v

newtype Intensity = Intensity { _getIntensity :: Double }
  deriving (Eq, Ord, Num, Fractional)

makeClassy ''Intensity

instance Monoid Intensity where
  mempty = 0
  mappend = (+)

instance Show Intensity where
  show (Intensity v) = show v

newtype RelativeAbundance = RelativeAbundance {
  _getRelativeAbundance :: Double
} deriving (Eq, Ord, Num, Fractional) -- write show instance, 2 d.p, add "%"

makeClassy ''RelativeAbundance

instance Monoid RelativeAbundance where
  mempty = 0
  mappend = (+)

instance Show RelativeAbundance where
  show (RelativeAbundance v) = show v

newtype NormalisedAbundance =  NormalisedAbundance {
  _getNormalisedAbundance :: Double -- write show instance, 2 d.p, add "%"
} deriving (Eq, Ord, Num, Fractional)

makeClassy ''NormalisedAbundance

instance Monoid NormalisedAbundance where
  mempty = 0
  mappend = (+)

instance Show NormalisedAbundance where
  show (NormalisedAbundance v) = show v

data SpectrumRow a = SpectrumRow {
    _ion                   :: a
  , _srIntensity           :: Intensity
  , _srRelativeAbundance   :: RelativeAbundance
  , _srNormalisedAbundance :: NormalisedAbundance
} deriving (Eq, Ord, Functor, Foldable)

makeClassy ''SpectrumRow

instance HasIntensity (SpectrumRow a) where
  intensity = srIntensity

instance HasNormalisedAbundance (SpectrumRow a) where
  normalisedAbundance = srNormalisedAbundance

instance HasRelativeAbundance (SpectrumRow a) where
  relativeAbundance = srRelativeAbundance

instance Show a => Show (SpectrumRow a) where
  show (SpectrumRow m i r n) =
    intercalate ", " [show m, show i, show r, show n]

newtype MSSpectrum a = MSSpectrum
  { _msSpectrum :: [SpectrumRow a] }
  deriving (Eq, Ord, Functor, Foldable)

makeLenses ''MSSpectrum

data SpectrumCsv = SpectrumCsv {
    _csvMz        :: Mz
  , _csvIntensity :: Intensity
}

makeLenses ''SpectrumCsv

instance FromNamedRecord SpectrumCsv where
  parseNamedRecord m =
    (\mz' i -> SpectrumCsv (Mz mz') (Intensity i))
      <$> m .: "m/z"
      <*> m .: "Intensity"

-- Adapted from https://github.com/Gabriel439/slides/blob/master/lambdaconf/data/exercises/02/Main.hs
process :: FilePath -> IO [SpectrumCsv]
process file = do
    bytes <- Data.ByteString.Lazy.readFile file
    case decodeByName bytes of
        Left   err        -> fail err
        Right (_, vector) -> return $ toList vector

insertAbundances :: [SpectrumCsv] -> [SpectrumRow Mz]
insertAbundances spec =
  fmap (\csv ->
    SpectrumRow
      (csv^.csvMz)
      (csv^.csvIntensity)
      (RelativeAbundance (csv^.csvIntensity.getIntensity * 100 / maxIntensity))
      (NormalisedAbundance ((csv^.csvIntensity.getIntensity) * 100 / totalIntensity)))
  spec
  where
    totalIntensity = foldOf (folded.csvIntensity) spec ^.getIntensity
    maxIntensity =
      maximumBy (comparing _csvIntensity) spec ^.csvIntensity.getIntensity

mkMSSpectrum :: [SpectrumCsv] -> MSSpectrum Mz
mkMSSpectrum = MSSpectrum . insertAbundances

showList' :: (Show a) => [a] -> String -- maybe delete this function
showList' l = intercalate "\n" $ show <$> l

instance Show a => Show (MSSpectrum a) where
  show (MSSpectrum s) = "MS spectrum \n" <> showList' s

data MS2Spectrum a b = MS2Spectrum {
    _precursorIon :: a
  , _ms2Spectrum :: [SpectrumRow b]
} deriving (Show, Eq, Ord, Functor, Foldable)

makeLenses ''MS2Spectrum

lookupIon :: (Num b, Ord b) => b -> b -> MS2Spectrum a b -> Maybe (SpectrumRow b)
lookupIon n i (MS2Spectrum _ rs) = loop n i rs
  where
    loop n' i' rs' =
      case rs' of
        [] -> Nothing
        r':rs'' -> if withinTolerance n' i' (r'^.ion)
                     then Just r'
                     else loop n' i' rs''

removePrecursorIon :: Mz -> [SpectrumCsv] -> [SpectrumCsv]
removePrecursorIon prec =
  filter (^.csvMz.to (< prec - 2))

mkMS2SpectrumRemovePrecursor :: Mz -> [SpectrumCsv] -> MS2Spectrum Mz Mz
mkMS2SpectrumRemovePrecursor mz' rs =
  MS2Spectrum mz' . insertAbundances $ removePrecursorIon mz' rs

filterByRelativeAbundance ::
  (RelativeAbundance -> Bool) -> MS2Spectrum a b -> MS2Spectrum a b
filterByRelativeAbundance p =
  ms2Spectrum %~ filter (^.relativeAbundance.to p)

filterByNormalisedAbundance ::
  (NormalisedAbundance -> Bool) -> MS2Spectrum a b -> MS2Spectrum a b
filterByNormalisedAbundance p =
  ms2Spectrum %~ filter (^.normalisedAbundance.to p)

calNeutralLoss :: Mz -> Mz -> MonoisotopicMass
calNeutralLoss prec frag = MonoisotopicMass . getMz $ prec - frag

calNeutralLosses :: MS2Spectrum Mz Mz -> MS2Spectrum Mz MonoisotopicMass
calNeutralLosses spec = calNeutralLoss (spec^.precursorIon) <$> spec

withinTolerance :: (Num a, Ord a) => a -> a -> a -> Bool
withinTolerance n v1 v2 = n > abs (v1 - v2)

findPrecursorIon :: Mz -> Mz -> MSSpectrum Mz -> Maybe (SpectrumRow Mz)
findPrecursorIon n m (MSSpectrum spec) =
  maximumByOf traverse (comparing (^.intensity)) filteredSpectrumRows
  where
    filteredSpectrumRows :: [SpectrumRow Mz]
    filteredSpectrumRows = filter
                             (\x -> withinTolerance n m (x ^. ion))
                             spec
