{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
module Results where

import Control.Lens hiding ((.=))
import FattyAcid
import Triacylglycerol
import Spectra
import Isotope.Ion
import Data.Csv ( ToNamedRecord(..)
                , ToField(..)
                , DefaultOrdered(..)
                , namedRecord
                , (.=)
                , header
                )
import Data.Map (Map)
import qualified Data.Map as Map
import Data.Maybe (fromJust, isJust)
import Data.List (sort, nub)


data FinalResult = FinalResult {
    _finalResultMz                    :: !Mz
  , _finalResultMzNormalisedAbundance :: !NormalisedAbundance
  , _finalResultFAs                   :: Map FattyAcid NormalisedAbundance
  , _finalResultTags                  :: ![Triacylglycerol]
} deriving (Show, Eq, Ord)

makeLenses ''FinalResult

instance HasNormalisedAbundance FinalResult where
  normalisedAbundance = finalResultMzNormalisedAbundance
--
toFinalResult
  :: Intensity
  -> MS2Spectrum (Mz, Intensity, [Triacylglycerol]) (Maybe FattyAcid)
  -> FinalResult
toFinalResult totalIntensity spec =
  FinalResult (spec^.precursorIon._1)
              (NormalisedAbundance correctionRatio')
              (Map.fromList
                (fmap
                  (\x -> (x^.ion.to fromJust, x^.normalisedAbundance))
                spectrumRows))
              (spec^.precursorIon._3)
  where
    spectrumRows :: [SpectrumRow (Maybe FattyAcid)]
    spectrumRows = spec^..ms2Spectrum.folded.filtered (\x -> x^.ion.to isJust)

    correctionRatio' :: Double
    correctionRatio' = correctionRatio (spec^.precursorIon._2) totalIntensity

--
toFinalResults
  :: [MS2Spectrum (Mz, Intensity, [Triacylglycerol]) (Maybe FattyAcid)]
  -> [FinalResult]
toFinalResults specs = toFinalResult totalIntensity <$> specs
  where
    totalIntensity = foldOf (folded.precursorIon._2) specs

data IdentifiedFA = IdentifiedFA {
    identifiedFA            :: !FattyAcid
  , identifiedFAAbund       :: !NormalisedAbundance
  , presentInIdentifiedTAGs :: !Bool
} deriving (Show, Eq, Ord)

instance ToNamedRecord IdentifiedFA where
  toNamedRecord IdentifiedFA{..} =
    namedRecord
      [ "Fatty Acid"                 .= identifiedFA
      , "Normalised Abundance"       .= identifiedFAAbund
      , "Present in identified TAGs" .= presentInIdentifiedTAGs
      ]

instance ToField Bool where
  toField b =
    if b then "yes" else "no"

instance DefaultOrdered IdentifiedFA where
  headerOrder _ =
    header
      [ "Fatty Acid"
      , "Normalised Abundance"
      , "Present in identified TAGs"
      ]

instance IsSaturated IdentifiedFA where
  isSaturated = isSaturated . identifiedFA

instance IsMonounsaturated IdentifiedFA where
  isMonounsaturated = isMonounsaturated . identifiedFA

instance IsPolyunsaturated IdentifiedFA where
  isPolyunsaturated = isPolyunsaturated . identifiedFA

identifiedFAs :: [FinalResult] -> (NormalisedAbundance, [IdentifiedFA])
identifiedFAs frs =
  (percentageIonsIdentified,
  (\(fa, na) -> IdentifiedFA fa na (fa `elem` allTagFAs)) <$>
  Map.toList normalisedFattyAcidAbundances)
    where
      getRatio :: FinalResult -> Double
      getRatio fr = fr^.getNormalisedAbundance

      reCalNormalisedAbundance :: FinalResult -> Map FattyAcid NormalisedAbundance
      reCalNormalisedAbundance fr =
        mapped.getNormalisedAbundance %~ (* getRatio fr) $ fr^.finalResultFAs

      fattyAcidAbundances :: Map FattyAcid NormalisedAbundance
      fattyAcidAbundances =
        foldr
          (Map.unionWith mappend)
          Map.empty
          (reCalNormalisedAbundance <$> frs)

      percentageIonsIdentified :: NormalisedAbundance
      percentageIonsIdentified = sum fattyAcidAbundances

      normalisedFattyAcidAbundances :: Map FattyAcid NormalisedAbundance
      normalisedFattyAcidAbundances =
        mapped.getNormalisedAbundance %~
          (\n ->  n * 100 / percentageIonsIdentified^.getNormalisedAbundance) $
          fattyAcidAbundances

      allTagFAs :: [FattyAcid]
      allTagFAs = FattyAcid <$>
                  (sort . nub $ frs^..traverse.finalResultTags.traverse.allFAs)

data TAGSummary = TAGSummary {
    msIon                        :: !Mz
  , tagSummayNormalisedAbundance :: !NormalisedAbundance
  , condensedTriacylglycerols    :: ![CondensedTriacylglycerol]
  , triacylglycerols             :: ![Triacylglycerol]
} deriving (Show, Eq, Ord)

instance ToNamedRecord TAGSummary where
  toNamedRecord TAGSummary{..} =
    namedRecord
      [ "M/z"                            .= msIon
      , "Normalised Abundance"           .= tagSummayNormalisedAbundance
      , "Species-level Triacylglycerols" .= condensedTriacylglycerols
      , "Alkyl-level Triacylglycerols"   .= triacylglycerols
      ]

instance ToField Mz where
  toField (Mz n) = toField n

instance DefaultOrdered TAGSummary where
  headerOrder _ =
    header
      [ "M/z"
      , "Normalised Abundance"
      , "Species-level Triacylglycerols"
      , "Alkyl-level Triacylglycerols"
      ]

tagSummary :: [FinalResult] -> [TAGSummary]
tagSummary =
  fmap (\fr ->
    TAGSummary
    (fr^.finalResultMz)
    (fr^.normalisedAbundance.to (*100))
    (nub . tagsToCondensedTags $ fr^.finalResultTags)
    (fr^.finalResultTags))
