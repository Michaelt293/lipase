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

newtype Intensity = Intensity { _getIntensity :: Double }
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

newtype MSSpectrum = MSSpectrum
  { _msSpectrum :: [SpectrumRow] }
  deriving (Eq, Ord)

makeLenses ''MSSpectrum

sumIntensities' :: [(Mz, Intensity)] -> Intensity
sumIntensities' = foldr (\(_, i) acc -> acc + i) 0

normalizedAbundances :: [(Mz, Intensity)] -> [NormalisedAbundance]
normalizedAbundances spec =
  (\(_, i) ->
      NormalisedAbundance
        (_getIntensity i * 100 / _getIntensity (sumIntensities' spec)))
        <$> spec

relativeAbundances :: [(Mz, Intensity)] -> [RelativeAbundance]
relativeAbundances spec =
  (\(_, i) ->
    RelativeAbundance
      (_getIntensity i * 100 / maximumBy (comparing snd) spec ^._2 ^. getIntensity))
      <$> spec

insertAbundances :: [(Mz, Intensity)] -> [SpectrumRow]
insertAbundances spec =
  zipWith3 (\(a, b) c d -> SpectrumRow a (IonInfo b c d))
    spec
    (relativeAbundances spec)
    (normalizedAbundances spec)

filterNormalisedAbundance :: NormalisedAbundance -> [SpectrumRow] -> [SpectrumRow]
filterNormalisedAbundance abun = filter (\r -> r ^. normalisedAbundance > abun)

filterRelativeAbundance :: RelativeAbundance -> [SpectrumRow] -> [SpectrumRow]
filterRelativeAbundance abun = filter (\r -> r ^. relativeAbundance > abun)

showList' :: (Show a) => [a] -> String
showList' l = intercalate "\n" $ show <$> l

class Spectrum a where
  smap :: ([SpectrumRow] -> [SpectrumRow]) -> a -> a
  filterByRelativeAbundance :: RelativeAbundance -> a -> a
  filterByNormalisedAbundance :: NormalisedAbundance -> a -> a
  selectRange :: Mz -> Mz -> a -> a
  recalculateAbundances :: a -> a
  filterByRelativeAbundance a = smap (filterRelativeAbundance a)
  filterByNormalisedAbundance a = smap (filterNormalisedAbundance a)
  selectRange min' max' spec =
    recalculateAbundances $
      smap (filter (\(SpectrumRow m _) -> m < max' && m > min')) spec
  recalculateAbundances =
    smap (insertAbundances . fmap (\(SpectrumRow m i) -> (m, i ^. intensity)))

instance Spectrum MSSpectrum where
  smap = over msSpectrum

instance Show MSSpectrum where
  show (MSSpectrum s) = "MS spectrum \n" ++ showList' s

data MS2Spectrum = MS2Spectrum {
    _precursorMz :: Mz
  , _ms2Spectrum :: [SpectrumRow]
} deriving (Eq, Ord)

makeLenses ''MS2Spectrum

instance HasMonoisotopicMass MS2Spectrum where
  monoisotopicMass = precursorMz

instance Show MS2Spectrum where
  show (MS2Spectrum p s) =
     "MS2 spectrum \n" ++ "precursor ion: " ++ show p ++ "\n" ++
     showList' s

instance Spectrum MS2Spectrum where
  smap = over ms2Spectrum

removePrecursorIon :: MS2Spectrum -> MS2Spectrum
removePrecursorIon spec =
  selectRange (MonoisotopicMass 0)
              (spec ^. monoisotopicMass |-| MonoisotopicMass 2) spec

calNeutralLosses :: MS2Spectrum -> [MonoisotopicMass]
calNeutralLosses (MS2Spectrum p s) =
  fmap (\r -> p |-| r ^. mz) s

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

makeLenses ''NeutralLossSpectrum

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
   NeutralLossRow (prec |-| m) i) <$> spec)

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
    filteredSpectrumRows = filter
                             (\x -> withinTolerance n m (x ^. mz))
                             spec
