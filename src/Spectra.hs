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

import qualified Isotope as I
import Isotope hiding (monoisotopicMass)
import Data.List
import Data.Ord
import Data.Monoid
import Data.Bifunctor
import Data.Csv
import Data.Text (Text)
import Data.Vector (Vector, toList)
import Numeric
import Control.Lens

import qualified Data.ByteString.Lazy
import qualified Data.Csv

makeClassy ''MonoisotopicMass

class ShowVal a where
  showVal :: a -> String

newtype Charge = Charge { getCharge :: Int }
  deriving (Show, Eq, Ord, Num)

newtype Mz = Mz { _getMz :: Double }
  deriving (Show, Eq, Ord, Num, Fractional)

makeLenses ''Mz

class ToElementalComposition a => Ion a where
  charge                  :: a -> Charge
  mz                      :: a -> Mz
  mz a = Mz $ monoisotopicMass' / charge'
    where
      monoisotopicMass' =
        getMonoisotopicMass . I.monoisotopicMass . toElementalComposition $ a
      charge' = fromIntegral . getCharge . charge $ a
  {-# MINIMAL charge #-}

deriving instance Num MonoisotopicMass

deriving instance Fractional MonoisotopicMass

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
  _getRelativeAbundance :: Double
} deriving (Show, Eq, Ord, Num, Fractional) -- write show instance, 2 d.p, add "%"

makeClassy ''RelativeAbundance

instance Monoid RelativeAbundance where
  mempty = 0
  mappend = (+)

instance ShowVal RelativeAbundance where
  showVal (RelativeAbundance v) = show v

newtype NormalisedAbundance =  NormalisedAbundance {
  _getNormalisedAbundance :: Double -- write show instance, 2 d.p, add "%"
} deriving (Show, Eq, Ord, Num, Fractional)

makeClassy ''NormalisedAbundance

instance Monoid NormalisedAbundance where
  mempty = 0
  mappend = (+)

instance ShowVal NormalisedAbundance where
  showVal (NormalisedAbundance v) = show v

data SpectrumRow a = SpectrumRow {
    _ion                   :: a
  , _srIntensity           :: Intensity
  , _srRelativeAbundance   :: RelativeAbundance
  , _srNormalisedAbundance :: NormalisedAbundance
} deriving (Eq, Ord, Functor, Foldable)

makeClassy ''SpectrumRow

-- instance HasMonoisotopicMass (SpectrumRow a) where
--   monoisotopicMass = mz

instance HasIntensity (SpectrumRow a) where
  intensity = srIntensity

instance HasNormalisedAbundance (SpectrumRow a) where
  normalisedAbundance = srNormalisedAbundance

instance HasRelativeAbundance (SpectrumRow a) where
  relativeAbundance = srRelativeAbundance

instance Show a => Show (SpectrumRow a) where
  show (SpectrumRow m i r n) =
    intercalate ", " [show m, show i, show m, show n]

newtype MSSpectrum a = MSSpectrum
  { _msSpectrum :: [SpectrumRow a] }
  deriving (Eq, Ord, Functor, Foldable)

makeLenses ''MSSpectrum

-- sumIntensities :: [(Mz, Intensity)] -> Intensity
-- sumIntensities = foldOf (folded._2) -- foldr (\(_, i) acc -> acc + i) 0

data SpectrumCsv = SpectrumCsv {
    _csvMz        :: Mz
  , _csvIntensity :: Intensity
}

makeLenses ''SpectrumCsv

instance FromNamedRecord SpectrumCsv where
  parseNamedRecord m =
    (\mz i -> SpectrumCsv (Mz mz) (Intensity i))
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


-- relativeAbundances :: [(Mz, Intensity)] -> [RelativeAbundance]
-- relativeAbundances spec =
--   (\(_, i) ->
--     <$> spec
--   where
--
--
-- insertAbundances spec =
--   zipWith3 (\(a, b) c d -> SpectrumRow a b c d)
--     spec
--     (relativeAbundances spec)
--     (normalizedAbundances spec)


showList' :: (Show a) => [a] -> String -- maybe delete this function
showList' l = intercalate "\n" $ show <$> l

-- class Spectrum a where
--   smap :: ([SpectrumRow] -> [SpectrumRow]) -> a -> a
--   filterByRelativeAbundance :: RelativeAbundance -> a -> a
--   filterByNormalisedAbundance :: NormalisedAbundance -> a -> a
--   recalculateAbundances :: a -> a
--   filterByRelativeAbundance a = smap (filterRelativeAbundance a)
--   filterByNormalisedAbundance a = smap (filterNormalisedAbundance a)
--   recalculateAbundances =
--     smap (insertAbundances . fmap (\(SpectrumRow m i _ _) -> (m, i ^. intensity)))

-- instance Spectrum MSSpectrum where
--   smap = over msSpectrum

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
    loop n i rs =
      case rs of
        [] -> Nothing
        r:rs' -> if withinTolerance n i (r^.ion)
                   then Just r
                   else loop n i rs'

removePrecursorIon :: Mz -> [SpectrumCsv] -> [SpectrumCsv]
removePrecursorIon prec =
  filter (^.csvMz.to (< prec - 2))

mkMS2SpectrumRemovePrecursor :: Mz -> [SpectrumCsv] -> MS2Spectrum Mz Mz
mkMS2SpectrumRemovePrecursor mz rs =
  MS2Spectrum mz . insertAbundances $ removePrecursorIon mz rs

filterByRelativeAbundance ::
  (RelativeAbundance -> Bool) -> MS2Spectrum a b -> MS2Spectrum a b
filterByRelativeAbundance p =
  ms2Spectrum %~ filter (^.relativeAbundance.to p)

filterByNormalisedAbundance ::
  (NormalisedAbundance -> Bool) -> MS2Spectrum a b -> MS2Spectrum a b
filterByNormalisedAbundance p =
  ms2Spectrum %~ filter (^.normalisedAbundance.to p)



-- instance HasMonoisotopicMass (MS2Spectrum a b) where
--   monoisotopicMass = ion

--instance Show (MS2Spectrum a) where
--  show (MS2Spectrum p s) =
--     "MS2 spectrum \n" <>
--     "precursor ion: " <> show p <> "\n" <>
--     showList' s

-- instance Spectrum MS2Spectrum where
--   smap = over ms2Spectrum

--removePrecursorIon :: MS2Spectrum Mz Mz -> MS2Spectrum Mz Mz
--removePrecursorIon spec =
--  ms2Spectrum %~ filter
--    (^.ion.to (< spec^.precursorIon |-| MonoisotopicMass 2)) $ spec
--  --(MonoisotopicMass 0)
  --            (spec ^. monoisotopicMass |-| MonoisotopicMass 2) spec


--calNeutralLosses :: MS2Spectrum Mz Mz -> [MonoisotopicMass]
--calNeutralLosses (MS2Spectrum p s) =
--  fmap (\r -> p |-| r ^. ion) s
--
-- data NeutralLossRow = NeutralLossRow {
--     _neutralLoss           :: MonoisotopicMass
--   , _nlIntensity           :: Intensity
--   , _nlRelativeAbundance   :: RelativeAbundance
--   , _nlNormalisedAbundance :: NormalisedAbundance
-- } deriving (Eq, Ord)
--
-- makeClassy ''NeutralLossRow
--
-- instance Show NeutralLossRow where
--   show (NeutralLossRow m i r n) =
--     intercalate ", " [show n, show i, show r, show n]
--
-- data NeutralLossSpectrum = NeutralLossSpectrum {
--     _nlPrecursorMz :: Mz
--   , _getNeutralLossSpectrum :: [NeutralLossRow]
--   } deriving (Eq, Ord)
--
-- makeLenses ''NeutralLossSpectrum
--
-- instance HasMonoisotopicMass NeutralLossSpectrum where
--   monoisotopicMass = nlPrecursorMz
--
-- instance Show NeutralLossSpectrum where
--   show (NeutralLossSpectrum n s) =
--      "Neutral loss spectrum \n" <>
--      "precursor ion: " <> show n <> "\n" <>
--      showList' s
--
-- toNeutralLossSpectrum :: MS2Spectrum Mz Mz -> NeutralLossSpectrum
-- toNeutralLossSpectrum (MS2Spectrum prec spec) = NeutralLossSpectrum prec
--  ((\(SpectrumRow m i r n) ->
--    NeutralLossRow (prec |-| m) i r n) <$> spec)

calNeutralLoss :: Mz -> Mz -> MonoisotopicMass
calNeutralLoss prec frag = (prec - frag)^.getMz.to MonoisotopicMass

calNeutralLosses :: MS2Spectrum Mz Mz -> MS2Spectrum Mz MonoisotopicMass
calNeutralLosses spec = calNeutralLoss (spec^.precursorIon) <$> spec

-- Reads spectrum from CSV file
-- readSpectrum :: [[String]] -> [(Mz, Intensity)]
-- readSpectrum spec = ionAndAbundance <$> spec
--   where
--     ionAndAbundance :: [String] -> (Mz, Intensity)
--     ionAndAbundance line = case line of
--       ion : abund : _ -> (MonoisotopicMass
--         (read ion), Intensity (read abund))
--       _ -> error "not a valid mass spectrum"

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
