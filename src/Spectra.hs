{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TemplateHaskell #-}

module Spectra where

import Isotope hiding (monoisotopicMass)
import Data.List
import Data.Ord
import Numeric
import Control.Lens

makeClassy ''MonoisotopicMass

class ShowVal a where
  showVal :: a -> String

type Mz = MonoisotopicMass

instance ShowVal MonoisotopicMass where
  showVal (MonoisotopicMass v) = show v

newtype Intensity = Intensity { getIntensity :: Double }
  deriving (Show, Eq, Ord, Num, Fractional)

makeClassy ''Intensity

instance Monoid Intensity where
  mempty = 0
  mappend = (+)

instance ShowVal Intensity where
  showVal (Intensity v) = show v

-- formatFloatN floatNum numOfDecimal = showFFloat (Just numOfDecimals) floatNum ""

newtype RelativeAbundance = RelativeAbundance {
  getRelativeAbundance :: Double
} deriving (Show, Eq, Ord, Num, Fractional) -- write show instance, 2 d.p, add "%"

makeClassy ''RelativeAbundance

instance Monoid RelativeAbundance where
  mempty = 0
  mappend = (+)

instance ShowVal RelativeAbundance where
  showVal (RelativeAbundance v) = show v

newtype NormalisedAbundance =  NormalisedAbundance {
  getNormalisedAbundance :: Double -- write show instance, 2 d.p, add "%"
} deriving (Show, Eq, Ord, Num, Fractional)

makeClassy ''NormalisedAbundance

instance Monoid NormalisedAbundance where
  mempty = 0
  mappend = (+)

instance ShowVal NormalisedAbundance where
  showVal (NormalisedAbundance v) = show v

data IonInfo = IonInfo {
    _ionInfoIntensity :: Intensity
  , _ionInfoRelativeAbundance :: RelativeAbundance
  , _ionInfoNormalisedAbundance :: NormalisedAbundance
} deriving (Eq, Ord)

makeClassy ''IonInfo

instance HasIntensity IonInfo where
  intensity = ionInfoIntensity

instance HasRelativeAbundance IonInfo where
  relativeAbundance = ionInfoRelativeAbundance

instance HasNormalisedAbundance IonInfo where
  normalisedAbundance = ionInfoNormalisedAbundance

instance Show IonInfo where
  show (IonInfo i r n) = intercalate ", " [show i, show r, show n]

data SpectrumRow = SpectrumRow {
    _mz :: Mz
  , _spectrumRowIonInfo :: IonInfo
} deriving (Eq, Ord)

makeClassy ''SpectrumRow

instance HasMonoisotopicMass SpectrumRow where
  monoisotopicMass = mz

instance HasIonInfo SpectrumRow where
  ionInfo = spectrumRowIonInfo

instance HasNormalisedAbundance SpectrumRow where
  normalisedAbundance = ionInfo.normalisedAbundance

instance HasRelativeAbundance SpectrumRow where
  relativeAbundance = ionInfo.relativeAbundance

instance HasIntensity SpectrumRow where
  intensity = ionInfo.intensity

instance Show SpectrumRow where
  show (SpectrumRow m i) = intercalate ", " [show m, show i]

newtype SpectrumRows = SpectrumRows { _getSpectrumRows :: [SpectrumRow] }
  deriving (Show, Eq, Ord)

makeClassy ''SpectrumRows

data MSSpectrum = MSSpectrum
  { _msSpectrumSpectrumRows :: SpectrumRows }
  deriving (Eq, Ord)

makeClassy ''MSSpectrum

instance HasSpectrumRows MSSpectrum where
  spectrumRows = msSpectrumSpectrumRows

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

insertAbundances :: [(Mz, Intensity)] -> SpectrumRows
insertAbundances spec = SpectrumRows $
  zipWith3 (\(a, b) c d -> SpectrumRow a (IonInfo b c d))
    spec
    (relativeAbundances spec)
    (normalizedAbundances spec)

filterNormalisedAbundance :: NormalisedAbundance -> SpectrumRows -> SpectrumRows
filterNormalisedAbundance abun spec = SpectrumRows $
  filter (\r -> r ^. normalisedAbundance > abun) (spec ^. getSpectrumRows)

filterRelativeAbundance :: RelativeAbundance -> SpectrumRows -> SpectrumRows
filterRelativeAbundance abun spec = SpectrumRows $
  filter (\r -> r ^. relativeAbundance > abun) (spec ^. getSpectrumRows)

showList' :: (Show a) => [a] -> String
showList' l = intercalate "\n" $ show <$> l

class Spectrum a where
  smap :: (SpectrumRows -> SpectrumRows) -> a -> a
  filterByRelativeAbundance :: RelativeAbundance -> a -> a
  filterByNormalisedAbundance :: NormalisedAbundance -> a -> a
  selectRange :: Mz -> Mz -> a -> a
  recalculateAbundances :: a -> a
  filterByRelativeAbundance a = smap (filterRelativeAbundance a)
  filterByNormalisedAbundance a = smap (filterNormalisedAbundance a)
  selectRange min' max' spec =
    recalculateAbundances $
      smap (\x -> SpectrumRows
                  (filter (\(SpectrumRow m _) -> m < max' && m > min')
                  (x ^. getSpectrumRows))) spec
  recalculateAbundances =
    smap (\x -> insertAbundances
                (fmap (\(SpectrumRow m i) -> (m, i ^. intensity))
                                             (x ^. getSpectrumRows)))
  {-# MINIMAL (smap) #-}

instance Spectrum MSSpectrum where
  smap f spec = MSSpectrum $ f (spec ^. spectrumRows)

instance Show MSSpectrum where
  show (MSSpectrum s) = "MS spectrum \n" ++ showList' (s ^. getSpectrumRows)

data MS2Spectrum = MS2Spectrum {
    _precursorMz :: Mz
  , _ms2SpectrumRows :: SpectrumRows
} deriving (Eq, Ord)

makeClassy ''MS2Spectrum

instance HasMonoisotopicMass MS2Spectrum where
  monoisotopicMass = precursorMz

instance HasSpectrumRows MS2Spectrum where
  spectrumRows = ms2SpectrumRows

instance Show MS2Spectrum where
  show (MS2Spectrum p s) =
     "MS2 spectrum \n" ++ "precursor ion: " ++ show p ++ "\n" ++
     showList' (s ^. getSpectrumRows)

instance Spectrum MS2Spectrum where
  smap f spec = MS2Spectrum (spec ^. precursorMz) $ f (spec ^. spectrumRows)

removePrecursorIon :: MS2Spectrum -> MS2Spectrum
removePrecursorIon spec =
  selectRange (MonoisotopicMass 0)
              (spec ^. monoisotopicMass |-| MonoisotopicMass 2) spec

calNeutralLosses :: MS2Spectrum -> [MonoisotopicMass]
calNeutralLosses (MS2Spectrum p s) =
  fmap (\r -> p |-| r ^. mz) (s ^. getSpectrumRows)

data NeutralLossRow = NeutralLossRow {
    _neutralLoss :: MonoisotopicMass
  , _neutralLossIonInfo :: IonInfo
} deriving (Eq, Ord)

makeClassy ''NeutralLossRow

instance Show NeutralLossRow where
  show (NeutralLossRow n i) = intercalate ", " [show n, show i]

data NeutralLossSpectrum = NeutralLossSpectrum {
    _nLPrecursorMz :: Mz
  , _getNeutralLossSpectrum :: [NeutralLossRow]
  } deriving (Eq, Ord)

makeClassy ''NeutralLossSpectrum

instance HasMonoisotopicMass NeutralLossSpectrum where
  monoisotopicMass = nLPrecursorMz

instance Show NeutralLossSpectrum where
  show (NeutralLossSpectrum n s) =
     "Neutral loss spectrum \n" ++
     "precursor ion: " ++ show n ++ "\n" ++
     showList' s

toNeutralLossSpectrum :: MS2Spectrum -> NeutralLossSpectrum
toNeutralLossSpectrum (MS2Spectrum prec spec) = NeutralLossSpectrum prec
 ((\(SpectrumRow m i) ->
   NeutralLossRow (prec |-| m) i) <$> (spec ^. getSpectrumRows))

-- Reads spectrum from CSV file
readSpectrum :: [[String]] -> [(Mz, Intensity)]
readSpectrum spec = ionAndAbundance <$> spec
  where
    ionAndAbundance :: [String] -> (Mz, Intensity)
    ionAndAbundance line = case line of
      ion : abund : _ -> (MonoisotopicMass
       (read ion), Intensity (read abund))
      _ -> error "not a valid mass spectrum"

withinTolerance :: Double -> MonoisotopicMass -> MonoisotopicMass -> Bool
withinTolerance n v1 v2 = n > abs (getMonoisotopicMass (v1 |-| v2))

findPrecursorIon :: Double -> Mz -> MSSpectrum -> Maybe SpectrumRow
findPrecursorIon n m (MSSpectrum spec) =
  if null filteredSpectrumRows
    then Nothing
    else Just $
      maximumBy (comparing (^. intensity)) filteredSpectrumRows
  where
    filteredSpectrumRows :: [SpectrumRow]
    filteredSpectrumRows = filter (\x -> withinTolerance n m (x ^. mz)) (spec ^. getSpectrumRows)
