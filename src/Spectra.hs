{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Spectra where

import Data.List

newtype Mz = Mz { getMz :: Double }
  deriving (Eq, Ord)

instance Show Mz where
  show (Mz n) = "m/z " ++ show n

newtype Dalton = Dalton { getDalton :: Double }
 deriving (Eq, Ord, Num, Fractional)

instance Show Dalton where
  show (Dalton n) = show n ++ " Da"

data IonAbundance = RelativeAbundance Double
                  | Intensity Double
                  deriving (Eq, Ord)

instance Show IonAbundance where
  show (RelativeAbundance n) = "Relative Abundance " ++ show n ++ "%"
  show (Intensity n) = "Intensity " ++ show n

data MSSpectrum = MSSpectrum { getMSSpectrum :: [(Mz, IonAbundance)] }
  deriving (Eq, Ord)

instance Show MSSpectrum where
  show (MSSpectrum s) = "MS spectrum \n" ++ renderPairList s

renderPairList :: (Show a, Show b) => [(a, b)] -> String
renderPairList ps = intercalate "\n" $ renderPair <$> ps
  where
    renderPair (a, b) = show a ++ ", " ++ show b

data MS2Spectrum = MS2Spectrum {
    mS2precursorMz :: Mz
  , getMS2Spectrum :: [(Mz, IonAbundance)]
} deriving (Eq, Ord)

instance Show MS2Spectrum where
  show (MS2Spectrum p s) =
     "MS2 spectrum \n" ++ "precursor ion: " ++ show p ++ renderPairList s

data NeutralLossSpectrum = NeutralLossSpectrum {
   nLPrecursorMz          :: Mz
 , getNeutralLossSpectrum :: [(Mz, IonAbundance)]
} deriving (Eq, Ord)

instance Show NeutralLossSpectrum where
 show (NeutralLossSpectrum p s) =
    "Neutral loss spectrum \n" ++
    "precursor ion: " ++ show p ++ "\n" ++
    renderPairList s

ionAndAbundance :: [String] -> (Mz, IonAbundance)
ionAndAbundance line = case line of
  ion : abund : _ -> (Mz (read ion), Intensity (read abund) )
  _ -> error "not a valid mass spectrum"

readSpectrum :: [[String]] -> [(Mz, IonAbundance)]
readSpectrum spec = ionAndAbundance <$> spec
