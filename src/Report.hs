{-# LANGUAGE OverloadedStrings #-}
module Report where

import FattyAcid
import Spectra
import Results
import Numeric (showFFloat)
import Data.List (sortBy, maximumBy, intersperse)
import Data.Ord (comparing, Down(..))
import Data.Monoid ((<>))
import Data.ByteString (ByteString)
import Data.ByteString.Char8 (pack)

nMostAbundantFAs :: Int -> [IdentifiedFA] -> [IdentifiedFA]
nMostAbundantFAs n =
  take n . sortBy (comparing (Down . identifiedFAAbund))

mostAbundantFA' :: (IdentifiedFA -> Bool) -> [IdentifiedFA] -> Maybe IdentifiedFA
mostAbundantFA' f fas =
  if null fas
    then Nothing
    else Just . maximumBy (comparing identifiedFAAbund) $ filter f fas

mostAbundantSaturatedFA :: [IdentifiedFA] -> Maybe IdentifiedFA
mostAbundantSaturatedFA =
  mostAbundantFA' isSaturated

mostAbundantMonounsaturatedFA :: [IdentifiedFA] -> Maybe IdentifiedFA
mostAbundantMonounsaturatedFA =
  mostAbundantFA' isMonounsaturated

mostAbundantFAPolyunsaturated :: [IdentifiedFA] -> Maybe IdentifiedFA
mostAbundantFAPolyunsaturated =
  mostAbundantFA' isPolyunsaturated

spectraWithIdentifiedTAGs :: [TAGSummary] -> [TAGSummary]
spectraWithIdentifiedTAGs =
  filter (not . null . condensedTriacylglycerols)

mostAbundantTAG :: [TAGSummary] -> Maybe TAGSummary
mostAbundantTAG tags =
  if null tags
    then Nothing
    else Just $
      maximumBy
      (comparing tagSummayNormalisedAbundance)
      (spectraWithIdentifiedTAGs tags)

normalisedAbundanceTotalTAGs :: [TAGSummary] -> NormalisedAbundance
normalisedAbundanceTotalTAGs =
  sum . fmap tagSummayNormalisedAbundance . spectraWithIdentifiedTAGs

-- Adapted from https://stackoverflow.com/questions/1559590/haskell-force-floats-to-have-two-decimals
formatFloatN :: RealFloat a => Int -> a-> String
formatFloatN numOfDecimals floatNum = showFFloat (Just numOfDecimals) floatNum ""

formatListToString :: Show a => [a] -> String
formatListToString l =
  case reverse (intersperse ", " (show <$> l)) of
    [] -> ""
    [x] -> x
    (x:_:xs) -> concat . reverse $ x : " and " : xs

formatListToString' :: [String] -> String
formatListToString' l =
  case reverse (intersperse ", " l) of
    [] -> ""
    [x] -> x
    (x:_:xs) -> concat . reverse $ x : " and " : xs

abundanceAssignedFAs :: RealFloat a => a -> ByteString
abundanceAssignedFAs a =
  "Of the neutral losses observed, " <>
  n <>
  "% (by normalised abundance) were assigned to fatty neutral losses. "
  where
    n = pack $ formatFloatN 2 a

faDescription :: [IdentifiedFA] -> ByteString
faDescription fas =
  case topThree of
    [] -> "No fatty acids were identified."
    [_] ->
      "The most abundant fatty acid identified was " <>
      topThreeFA <>
      " with a normalised abundance of " <>
      topThreeAbun <>
      "%, respectively. "
    _ ->
      "The " <>
      n <>
      " most abundant fatty acids identified were " <>
      topThreeFA <>
      " with normalised abundances of " <>
      topThreeAbun <>
      "%, respectively. " <>
      saturationFAs
  where
    topThree = nMostAbundantFAs 3 fas
    n = pack . show $ length topThree
    topThreeFA = pack . formatListToString $ identifiedFA <$> topThree
    topThreeAbun = pack . formatListToString' $ formatFloatN 2 . identifiedFAAbund <$> topThree
    saturationFAs =
      let sat = mostAbundantSaturatedFA fas
          mono = mostAbundantMonounsaturatedFA fas
          poly = mostAbundantFAPolyunsaturated fas
      in case (sat, mono, poly) of
        (Nothing, Nothing, Nothing) ->
          ""
        (Just fa, Nothing, Nothing) ->
          "The most abundant saturated fatty acids was " <>
          pack (show (identifiedFA fa)) <>
          " with an abundance of " <>
          pack (formatFloatN 2 (identifiedFAAbund fa)) <>
           "%. "
        (Nothing, Just fa, Nothing) ->
          "The most abundant monounsaturated fatty acids was " <>
          pack (show (identifiedFA fa)) <>
          " with an abundance of " <>
          pack (formatFloatN 2 (identifiedFAAbund fa)) <>
          "%. "
        (Nothing, Nothing, Just fa) ->
          "The most abundant polyunsaturated fatty acids was " <>
          pack (show (identifiedFA fa)) <>
          " with an abundance of " <>
          pack (formatFloatN 2 (identifiedFAAbund fa)) <>
          "%. "
        (Just sat', Just mono', Nothing) ->
          "The most abundant monounsaturated and saturated fatty acids were " <>
          pack (formatListToString (identifiedFA <$> [mono', sat'])) <>
          " with abundances of " <>
          pack (formatListToString' (formatFloatN 2 . identifiedFAAbund <$> [mono', sat'])) <>
          "%, respectively. "
        (Just sat', Nothing, Just poly') ->
          "The most abundant polyunsaturated, and saturated fatty acids were " <>
          pack (formatListToString (identifiedFA <$> [poly', sat'])) <>
          " with abundances of " <>
          pack (formatListToString' (formatFloatN 2 . identifiedFAAbund <$> [poly', sat'])) <>
          "%, respectively. "
        (Nothing, Just mono', Just poly') ->
          "The most abundant polyunsaturated, monounsaturated fatty acids were " <>
          pack (formatListToString (identifiedFA <$> [poly', mono'])) <>
          " with abundances of " <>
          pack (formatListToString' (formatFloatN 2 . identifiedFAAbund <$> [poly', mono'])) <>
          "%, respectively. "
        (Just sat', Just mono', Just poly') ->
          "The most abundant polyunsaturated, monounsaturated and saturated fatty acids were " <>
          pack (formatListToString (identifiedFA <$> [poly', mono', sat'])) <>
          " with abundances of " <>
          pack (formatListToString' (formatFloatN 2 . identifiedFAAbund <$> [poly', mono', sat'])) <>
          "%, respectively. "

identifedTAGInfo :: [TAGSummary] -> ByteString
identifedTAGInfo tags =
  "Triacylglycerols were identified for " <>
  pack (show (length (spectraWithIdentifiedTAGs tags))) <>
  " of the " <>
  pack (show (length tags)) <>
  " collision-induced dissociation spectra (" <>
  pack (formatFloatN 2 (normalisedAbundanceTotalTAGs tags)) <>
  "% by normalised abundance). " <>
  case mostAbundantTAG tags of
    Nothing -> "No triacylglycerols were identified."
    Just tag ->
      case condensedTriacylglycerols tag of
        [] -> "No triacylglycerols were identified."
        [tag'] -> "The most abundant triacylglycerol identified at the \
          \species level was " <>
          pack (show tag') <>
          " with a normalised abundance of " <>
          pack (formatFloatN 2 (tagSummayNormalisedAbundance tag)) <>
          "%. This species level triacylglycerol corresponds to the following \
          \possible alkyl bond level triacylglycerols: " <>
          pack (formatListToString (triacylglycerols tag)) <>
          "."
        tags' -> "The most abundant triacylglycerols identified at the \
          \species level were " <>
          pack (formatListToString tags') <>
          " with a combined normalised abundance of " <>
          pack ((formatFloatN 2 . tagSummayNormalisedAbundance) tag) <>
          "% (Note: Using low-resolution mass spectrometry, these species level \
          \triacylglycerols have the same m/z ratio). These species level \
          \triacylglycerols correspond to the following possible alkyl bond \
          \level triacylglycerols: " <>
          pack (formatListToString (triacylglycerols tag)) <>
          "."

report :: NormalisedAbundance -> [IdentifiedFA] -> [TAGSummary] -> ByteString
report a fas tags =
  "Fatty acids were identified from neutral losses observed in collision-induced \
  \dissociation spectra of protonated triacylglycerol ions. " <>
  abundanceAssignedFAs a <>
  faDescription fas <>
  "\n\nTriacylglycerols were identified by analysing the neutral losses observed in \
  \collision-induced dissociation spectra. " <>
  identifedTAGInfo tags
