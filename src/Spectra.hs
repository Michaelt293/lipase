{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Spectra where

import Isotope
import Data.List
import Data.Ord
import Numeric

newtype Mz = Mz { monoisotopicMass' :: MonoisotopicMass }
           deriving (Eq, Ord, Operators)

instance Show Mz where
  show (Mz n) = "m/z " ++ show (getMonoisotopicMass n)

mkMz :: Double -> Mz
mkMz = Mz . MonoisotopicMass

newtype Intensity = Intensity { getIntensity :: Double }
  deriving (Eq, Ord, Num, Fractional)

instance Monoid Intensity where
  mempty = 0
  mappend = (+)

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
    recalculateAbundances $ smap (filter (\(SpectrumRow m _) -> m < max' && m > min')) spec
  recalculateAbundances = smap (insertAbundances . fmap (\(SpectrumRow m i) -> (m, intensity i)))
  {-# MINIMAL (smap) #-}
    --insertAbundances (\(SpectrumRow m i _ _) -> (m, i) $

filterNormalisedAbundance :: NormalisedAbundance -> [SpectrumRow] -> [SpectrumRow]
filterNormalisedAbundance abun =
  filter (\(SpectrumRow _ i) -> normalisedAbundance i > abun)

filterRelativeAbundance :: RelativeAbundance -> [SpectrumRow] -> [SpectrumRow]
filterRelativeAbundance abun =
  filter (\(SpectrumRow _ i) -> relativeAbundance i > abun)

data IonInfo = IonInfo {
    intensity :: Intensity
  , relativeAbundance :: RelativeAbundance
  , normalisedAbundance :: NormalisedAbundance
} deriving (Eq, Ord)

instance Show IonInfo where
  show (IonInfo i r n) = intercalate ", " [show i, show r, show n]

data SpectrumRow = SpectrumRow {
    mz :: Mz
  , spectrumRowIonInfo :: IonInfo
} deriving (Eq, Ord)

instance Show SpectrumRow where
  show (SpectrumRow m i) = intercalate ", " [show m, show i]

data MSSpectrum = MSSpectrum
  { getMSSpectrum :: [SpectrumRow] }
  deriving (Eq, Ord)

instance Spectrum MSSpectrum where
  smap f spec = MSSpectrum $ f (getMSSpectrum spec)

instance Show MSSpectrum where
  show (MSSpectrum s) = "MS spectrum \n" ++ showList' s

showList' :: (Show a) => [a] -> String
showList' l = intercalate "\n" $ show <$> l

data MS2Spectrum = MS2Spectrum {
    mS2precursorMz :: Mz
  , getMS2Spectrum :: [SpectrumRow]
} deriving (Eq, Ord)

instance Show MS2Spectrum where
  show (MS2Spectrum p s) =
     "MS2 spectrum \n" ++ "precursor ion: " ++ show p ++ showList' s

instance Spectrum MS2Spectrum where
  smap f spec = MS2Spectrum (mS2precursorMz spec) $ f (getMS2Spectrum spec)

removePrecursorIon :: MS2Spectrum -> MS2Spectrum
removePrecursorIon spec =
  selectRange (mkMz 0)
              (mS2precursorMz spec |-| mkMz 2) spec

calNeutralLosses :: MS2Spectrum -> [MonoisotopicMass]
calNeutralLosses (MS2Spectrum p s) =
  (\(SpectrumRow mz _) -> monoisotopicMass' (p |-| mz)) <$> s

data NeutralLossRow = NeutralLossRow {
    neutralLoss :: MonoisotopicMass
  , neutralLossIonInfo :: IonInfo
} deriving (Eq, Ord)

instance Show NeutralLossRow where
  show (NeutralLossRow n i) = intercalate ", " [show n, show i]

data NeutralLossSpectrum = NeutralLossSpectrum {
    nLPrecursorMz       :: Mz
  , neutralLossSpectrum :: [NeutralLossRow]
  } deriving (Eq, Ord)

instance Show NeutralLossSpectrum where
  show (NeutralLossSpectrum n s) =
     "Neutral loss spectrum \n" ++
     "precursor ion: " ++ show n ++ "\n" ++
     showList' s

toNeutralLossSpectrum :: MS2Spectrum -> NeutralLossSpectrum
toNeutralLossSpectrum (MS2Spectrum prec spec) = NeutralLossSpectrum prec
 ((\(SpectrumRow m i) ->
   NeutralLossRow (monoisotopicMass' (prec |-| m)) i) <$> spec)

-- Reads spectrum from CSV file
readSpectrum :: [[String]] -> [(Mz, Intensity)]
readSpectrum spec = ionAndAbundance <$> spec
  where
    ionAndAbundance :: [String] -> (Mz, Intensity)
    ionAndAbundance line = case line of
      ion : abund : _ -> (mkMz (read ion), Intensity (read abund))
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
  zipWith3 (\(a, b) c d -> SpectrumRow a (IonInfo b c d))
    spec
    (relativeAbundances spec)
    (normalizedAbundances spec)

withinTolerance :: Double -> MonoisotopicMass -> MonoisotopicMass -> Bool
withinTolerance n v1 v2 = n > abs (getMonoisotopicMass (v1 |-| v2))

findPrecursorIon :: Double -> Mz -> MSSpectrum -> Maybe SpectrumRow
findPrecursorIon n m (MSSpectrum spec) =
  if null filteredSpectrumRows
    then Nothing
    else Just $
      maximumBy (comparing (intensity . spectrumRowIonInfo)) filteredSpectrumRows
  where
    filteredSpectrumRows = filter (\x -> withinTolerance n (monoisotopicMass' m)
                                  ((monoisotopicMass' . mz) x))
                                  spec
