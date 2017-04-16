module Lipase.TriacylglycerolSpec (spec) where

import Triacylglycerol
import FattyAcid
import Spectra
import Isotope
import Data.List
import Test.Hspec
import Test.QuickCheck

tg_181_181_181 = Triacylglycerol (FattyAcyl 18 1)
                                 (FattyAcyl 18 1)
                                 (FattyAcyl 18 1)

tg_161_181_201 = Triacylglycerol (FattyAcyl 16 1)
                                 (FattyAcyl 18 1)
                                 (FattyAcyl 20 1)
spec :: Spec
spec = do
  describe "Elemental composition for Triacylglycerol" $ do
    it "Trioleoylglycerol elemental composition should be C57H104O6" $
      toMolecularFormula tg_181_181_181
      `shouldBe` mkMolecularFormula [(C,57), (H, 104), (O, 6)]

  describe "possibleTAGs" $ do
    it "Trioleoylglycerol elemental composition should be C57H104O6" $
      sort (possibleTAGs 0.3 885.8 [ FattyAcid (FattyAcyl 16 1)
                                   , FattyAcid (FattyAcyl 18 1)
                                   , FattyAcid (FattyAcyl 20 1)
                                   ])
      `shouldBe` sort [ tg_181_181_181
                      , tg_161_181_201
                      ]

  describe "protonated triacylglycerol mz" $ do
    it "Trioleoylglycerol m/z should equal 885.8" $
      withinTolerance 0.1 885.8 (mz (Protonated tg_181_181_181))
