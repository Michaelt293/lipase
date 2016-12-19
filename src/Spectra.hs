{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Spectra where

import Data.List
import Data.Ord
import Numeric

newtype Mz = Mz { getMz :: Double }
  deriving (Eq, Ord, Num)

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

-- formatFloatN floatNum numOfDecimal = showFFloat (Just numOfDecimals) floatNum ""

newtype RelativeAbundance = RelativeAbundance Double
  deriving (Show, Eq, Ord, Num, Fractional) -- write show instance, 2 d.p, add "%"

newtype NormalisedAbundance =  NormalisedAbundance Double -- write show instance, 2 d.p, add "%"
  deriving (Show, Eq, Ord, Num, Fractional)

class Spectrum a where
  smap :: ([SpectrumRow] -> [SpectrumRow]) -> a -> a
  filterByRelativeAbundance :: RelativeAbundance -> a -> a
  filterByNormalisedAbundance :: NormalisedAbundance -> a -> a
  selectRange :: Mz -> Mz -> a -> a
  recalculateAbundances :: a -> a
  filterByRelativeAbundance a = smap (filterRelativeAbundance a)
  filterByNormalisedAbundance a = smap (filterNormalisedAbundance a)
  selectRange min' max' spec =
    recalculateAbundances $ smap (filter (\(SpectrumRow m _ _ _) -> m < max' && m > min')) spec
  recalculateAbundances = smap (insertAbundances . fmap (\(SpectrumRow m i _ _) -> (m, i)))
  {-# MINIMAL (smap) #-}
    --insertAbundances (\(SpectrumRow m i _ _) -> (m, i) $

filterNormalisedAbundance :: NormalisedAbundance -> [SpectrumRow] -> [SpectrumRow]
filterNormalisedAbundance abun =
  filter (\(SpectrumRow _ _ _ n) -> n > abun)

filterRelativeAbundance :: RelativeAbundance -> [SpectrumRow] -> [SpectrumRow]
filterRelativeAbundance abun =
  filter (\(SpectrumRow _ _ a _) -> a > abun)

data SpectrumRow = SpectrumRow {
    mz :: Mz
  , intensity :: Intensity
  , relativeAbundance :: RelativeAbundance
  , normalisedAbundance :: NormalisedAbundance
} deriving (Eq, Ord)

instance Show SpectrumRow where
  show (SpectrumRow m i r n) = intercalate ", " [show m, show i, show r, show n]

data MSSpectrum = MSSpectrum
  { getMSSpectrum :: [SpectrumRow] }
  deriving (Eq, Ord)

instance Spectrum MSSpectrum where
  smap f spec = MSSpectrum $ f (getMSSpectrum spec)

instance Show MSSpectrum where
  show (MSSpectrum s) = "MS spectrum \n" ++ show s

-- Delete this function if it remains unused
renderTupleList :: (Show a, Show b, Show c, Show d) => [(a, b, c, d)] -> String
renderTupleList ps = intercalate "\n" (render <$> ps)
  where
    render (a, b, c, d) = intercalate ", " [show a, show b, show c, show d]

data MS2Spectrum = MS2Spectrum {
    mS2precursorMz :: Mz
  , getMS2Spectrum :: [SpectrumRow]
} deriving (Eq, Ord)

instance Show MS2Spectrum where
  show (MS2Spectrum p s) =
     "MS2 spectrum \n" ++ "precursor ion: " ++ show p ++ show s

instance Spectrum MS2Spectrum where
  smap f spec = MS2Spectrum (mS2precursorMz spec) $ f (getMS2Spectrum spec)

removePrecursorIon :: MS2Spectrum -> MS2Spectrum
removePrecursorIon spec = selectRange 0 (mS2precursorMz spec - Mz 2) spec

calNeutralLosses :: MS2Spectrum -> [Dalton]
calNeutralLosses (MS2Spectrum p s) =
  (\(SpectrumRow mz _ _ _) -> Dalton (getMz p - getMz mz)) <$> s

data NeutralLossSpectrum = NeutralLossSpectrum {
    nLPrecursorMz          :: Mz
  , getNeutralLossSpectrum :: [(Dalton, Intensity, RelativeAbundance, NormalisedAbundance)]
  } deriving (Eq, Ord)

instance Show NeutralLossSpectrum where
  show (NeutralLossSpectrum p s) =
     "Neutral loss spectrum \n" ++
     "precursor ion: " ++ show p ++ "\n" ++
     show s

toNeutralLossSpectrum :: MS2Spectrum -> NeutralLossSpectrum
toNeutralLossSpectrum (MS2Spectrum prec spec) = NeutralLossSpectrum prec
 ((\(SpectrumRow mz i r n) -> (Dalton (getMz prec - getMz mz), i, r, n)) <$> spec)

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

normalizedAbundances :: [(Mz, Intensity)] -> [NormalisedAbundance]
normalizedAbundances spec =
  (\(_, i) ->
      NormalisedAbundance
        (getIntensity i * 100 / getIntensity (sumIntensities' spec)))
        <$> spec

relativeAbundances :: [(Mz, Intensity)] -> [RelativeAbundance]
relativeAbundances spec =
  (\(_, i) ->
    RelativeAbundance
      (getIntensity i * 100 / getIntensity (snd (maximumBy (comparing snd) spec))))
      <$> spec

insertAbundances :: [(Mz, Intensity)] -> [SpectrumRow]
insertAbundances spec =
  zipWith3 (\(a, b) c d -> SpectrumRow a b c d)
    spec
    (relativeAbundances spec)
    (normalizedAbundances spec)
