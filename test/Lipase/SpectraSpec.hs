module Lipase.SpectraSpec (spec) where

import Spectra
import FattyAcid
import Triacylglycerol
import Data.Maybe
import Data.List
import Control.Lens
import Test.Hspec

tg_161_181_201_spectrum :: [SpectrumCsv]
tg_161_181_201_spectrum =
  [ SpectrumCsv 575.5 33
  , SpectrumCsv 580.0 1
  , SpectrumCsv 603.5 33
  , SpectrumCsv 631.6 33
  , SpectrumCsv 885.8 22
  ]

testMsSpectrum :: MSSpectrum Mz
testMsSpectrum = mkMSSpectrum
    [ SpectrumCsv 623.5 33
    , SpectrumCsv 732.6 33
    , SpectrumCsv 885.8 1100
    , SpectrumCsv 887.8 1300
    , SpectrumCsv 689.9 1200
    ]

filteredSpectrum :: MS2Spectrum Mz Mz
filteredSpectrum = filterByRelativeAbundance (> 5)
                   (mkMS2SpectrumRemovePrecursor 885.5 tg_161_181_201_spectrum)

lookup_575 :: Maybe (SpectrumRow Mz)
lookup_575 = lookupIon 0.1 575.5 filteredSpectrum

spec :: Spec
spec = do
  describe "removePrecursorIon" .
    it "removePrecursorIon should remove the precursor ion" $
      all (\x -> x^.csvMz /= 885.8)
          (removePrecursorIon 885.8 tg_161_181_201_spectrum)

  describe "mkMS2SpectrumRemovePrecursor" $ do
    it "precursor removed" .
      isNothing $ lookupIon 0.1 885.8
          (mkMS2SpectrumRemovePrecursor 885.5 tg_161_181_201_spectrum)
    it "filterByRelativeAbundance - m/z 580.0 ion removed" .
      isNothing . lookupIon 0.1 580.0 $ filterByRelativeAbundance (> 5)
          (mkMS2SpectrumRemovePrecursor 885.5 tg_161_181_201_spectrum)
    it "relativeAbundance of m/z 575.5 ion is 100%" $
      withinTolerance 0.1 100 (fromMaybe 0 (lookup_575^?_Just.relativeAbundance))
    it "normalisedAbundance of m/z 575.5 ion is 33%" $
      withinTolerance 0.1 33 (fromMaybe 0 (lookup_575^?_Just.normalisedAbundance))
    it "intensity of m/z 575.5 ion is 33" $
      withinTolerance 0.1 33 (fromMaybe 0 (lookup_575^?_Just.intensity))

  describe "calNeutralLosses" $
    it "310.3 Da neutral loss present in the spectrum of TG 16:1_18:1_20:1"
      (isJust . lookupIon 0.4 310.3 $ calNeutralLosses filteredSpectrum)

  describe "assignFAsFromNeutralLoss" .
    it "Should identify the fatty acids that make up TG 16:1_18:1_20:1" $
      sort ((assignFAsFromNeutralLoss . calNeutralLosses) filteredSpectrum
        ^..ms2Spectrum.traverse.ion._Just) `shouldBe`
        [ FattyAcid (FattyAcyl 16 1)
        , FattyAcid (FattyAcyl 18 1)
        , FattyAcid (FattyAcyl 20 1)
        ]

  describe "findPossibleTAGs" .
    it "TG 16:1_18:1_20:1 spectrum should identify TG 16:1_18:1_20:1 and TG 18:1_18:1_18:1" $
      (findPossibleTAGs testMsSpectrum . assignFAsFromNeutralLoss . calNeutralLosses)
        filteredSpectrum^.precursorIon `shouldBe`
        (885.5, 1100,
        [ Triacylglycerol (FattyAcyl 16 1)
                           (FattyAcyl 18 1)
                           (FattyAcyl 20 1)
         , Triacylglycerol (FattyAcyl 18 1)
                           (FattyAcyl 18 1)
                           (FattyAcyl 18 1)
         ])
