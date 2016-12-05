{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Spectra where

import Data.List
import Data.Ord

newtype Mz = Mz { getMz :: Double }
  deriving (Eq, Ord)

instance Show Mz where
  show (Mz n) = "m/z " ++ show n

newtype Dalton = Dalton { getDalton :: Double }
  deriving (Eq, Ord, Num, Fractional)

instance Show Dalton where
  show (Dalton n) = show n ++ " Da"

newtype Intensity = Intensity { getIntensity :: Double }
  deriving (Eq, Ord, Num, Fractional)

instance Show Intensity where
  show (Intensity n) = "Intensity " ++ show n

data Abundance = RelativeAbundance Double
               | NormalisedAbundance Double
               deriving (Eq, Ord)

instance Show Abundance where
  show (RelativeAbundance a) = show a ++ "%"
  show (NormalisedAbundance a) = show a ++ "%"

class Spectrum a where
  getSpectrum :: a -> [(Mz, Intensity, Abundance)]
  sumIntensities :: a -> Intensity
  mostAbundantIon :: a -> (Mz, Intensity, Abundance)
  toNormalisedAbundances :: a -> [(Mz, Intensity, Abundance)]
  toRelativeAbundances :: a -> [(Mz, Intensity, Abundance)]
  filterByAbundance :: Abundance -> a -> [(Mz, Intensity, Abundance)]
  filterSpectrumByAbundance :: Abundance -> a -> a
  mostAbundantIon spec = maximumBy (comparing (\(_, i, _) -> i)) (getSpectrum spec)
  sumIntensities spec = foldr (\(_, i, _) acc -> acc + i) 0 (getSpectrum spec)
  toNormalisedAbundances spec =
    (\(mz, i, _) -> (mz, i,
      NormalisedAbundance
        (getIntensity i * 100 / getIntensity (sumIntensities spec))))
        <$> getSpectrum spec
  toRelativeAbundances spec =
    (\(mz, i, _) -> (mz, i,
      RelativeAbundance
        (getIntensity i * 100 / getIntensity ((\(_, a, _) -> a) (mostAbundantIon spec)))))
        <$> getSpectrum spec
  filterByAbundance abun spec =
    case abun of
      (NormalisedAbundance _) -> filter (\(_, _, a) -> a > abun) (toNormalisedAbundances spec)
      (RelativeAbundance _) -> filter (\(_, _, a) -> a > abun) (toRelativeAbundances spec)
  {-# MINIMAL (getSpectrum, filterSpectrumByAbundance) #-}

data MSSpectrum = MSSpectrum { getMSSpectrum :: [(Mz, Intensity, Abundance)] }
  deriving (Eq, Ord)

instance Spectrum MSSpectrum where
  getSpectrum = getMSSpectrum
  filterSpectrumByAbundance a spec = MSSpectrum (filterByAbundance a spec)

instance Show MSSpectrum where
  show (MSSpectrum s) = "MS spectrum \n" ++ renderTupleList s

renderTupleList :: (Show a, Show b, Show c) => [(a, b, c)] -> String
renderTupleList ps = intercalate "\n" $ renderPair <$> ps
  where
    renderPair (a, b, c) = show a ++ ", " ++ show b ++ ", " ++ show c

data MS2Spectrum = MS2Spectrum {
    mS2precursorMz :: Mz
  , getMS2Spectrum :: [(Mz, Intensity, Abundance)]
} deriving (Eq, Ord)

instance Show MS2Spectrum where
  show (MS2Spectrum p s) =
     "MS2 spectrum \n" ++ "precursor ion: " ++ show p ++ renderTupleList s

instance Spectrum MS2Spectrum where
  getSpectrum = getMS2Spectrum
  filterSpectrumByAbundance a spec =
    MS2Spectrum (mS2precursorMz spec)(filterByAbundance a spec)

calNeutralLosses :: MS2Spectrum -> [Dalton]
calNeutralLosses (MS2Spectrum p s) =
  (\(mz, _, _) -> Dalton (getMz p - getMz mz)) <$> s

data NeutralLossSpectrum = NeutralLossSpectrum {
    nLPrecursorMz          :: Mz
  , getNeutralLossSpectrum :: [(Dalton, Intensity, Abundance)]
  } deriving (Eq, Ord)

instance Show NeutralLossSpectrum where
  show (NeutralLossSpectrum p s) =
     "Neutral loss spectrum \n" ++
     "precursor ion: " ++ show p ++ "\n" ++
     renderTupleList s ++ "\n"

toNeutralLossSpectrum :: MS2Spectrum -> NeutralLossSpectrum
toNeutralLossSpectrum (MS2Spectrum prec spec) = NeutralLossSpectrum prec
 ((\(mz, i, a) -> (Dalton (getMz prec - getMz mz), i, a)) <$> spec)

-- Reads spectrum from CSV file
readSpectrum :: [[String]] -> [(Mz, Intensity)]
readSpectrum spec = ionAndAbundance <$> spec
  where
    ionAndAbundance :: [String] -> (Mz, Intensity)
    ionAndAbundance line = case line of
      ion : abund : _ -> (Mz (read ion), Intensity (read abund) )
      _ -> error "not a valid mass spectrum"

sumIntensities' :: [(Mz, Intensity)] -> Intensity
sumIntensities' = foldr (\(_, i) acc -> acc + i) 0

normalizedAbundances :: [(Mz, Intensity)] -> [Abundance]
normalizedAbundances spec =
  (\(_, i) ->
      NormalisedAbundance
        (getIntensity i * 100 / getIntensity (sumIntensities' spec)))
        <$> spec

insertAbundances :: [(Mz, Intensity)] -> [(Mz, Intensity, Abundance)]
insertAbundances spec =
  zipWith (\(a, b) c -> (a, b, c)) spec (normalizedAbundances spec)
