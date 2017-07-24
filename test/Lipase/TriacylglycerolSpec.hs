module Lipase.TriacylglycerolSpec (spec) where

import Triacylglycerol
import FattyAcid
import Spectra
import Isotope
import Isotope.Ion
import Data.List (sort)
import Test.Hspec

tg_180_180_180 :: Triacylglycerol
tg_180_180_180 = Triacylglycerol (FattyAcyl 18 0)
                                 (FattyAcyl 18 0)
                                 (FattyAcyl 18 0)

tg_180_180_181 :: Triacylglycerol
tg_180_180_181 = Triacylglycerol (FattyAcyl 18 0)
                                 (FattyAcyl 18 0)
                                 (FattyAcyl 18 1)

tg_181_181_181 :: Triacylglycerol
tg_181_181_181 = Triacylglycerol (FattyAcyl 18 1)
                                 (FattyAcyl 18 1)
                                 (FattyAcyl 18 1)

tg_161_181_201 :: Triacylglycerol
tg_161_181_201 = Triacylglycerol (FattyAcyl 16 1)
                                 (FattyAcyl 18 1)
                                 (FattyAcyl 20 1)
spec :: Spec
spec = do
  describe "Elemental composition for Triacylglycerol" .
    it "Trioleoylglycerol elemental composition should be C57H104O6" $
      toElementalComposition tg_181_181_181
      `shouldBe` mkElementalComposition [(C,57), (H, 104), (O, 6)]

  describe "possibleTAGs" .
    it "possibleTAGs 0.3 885.8 with 16:1, 18:1 and 20:1" $
      sort (possibleTAGs 0.3 885.8 [ FattyAcid (FattyAcyl 16 1)
                                   , FattyAcid (FattyAcyl 18 1)
                                   , FattyAcid (FattyAcyl 20 1)
                                   ])
      `shouldBe` sort [ tg_181_181_181
                      , tg_161_181_201
                      ]

  describe "protonated triacylglycerol mz" .
    it "Trioleoylglycerol m/z should equal 885.8" $
      withinTolerance 0.1 885.8 (mz (Protonated tg_181_181_181))

  describe "isSaturated" $ do
    it "TG 18:0_18:0_18:0 is saturated" $
      isSaturated tg_180_180_180
      `shouldBe` True
    it "TG 18:0_18:0_18:1 is not saturated" $
      isSaturated tg_180_180_181
      `shouldBe` False
    it "TG 18:1_18:1_18:1 is not saturated" $
      isSaturated tg_181_181_181
      `shouldBe` False

  describe "isMonounsaturated" $ do
    it "TG 18:0_18:0_18:0 is monounsaturated" $
      isMonounsaturated tg_180_180_180
      `shouldBe` False
    it "TG 18:0_18:0_18:1 is not monounsaturated" $
      isMonounsaturated tg_180_180_181
      `shouldBe` True
    it "TG 18:1_18:1_18:1 is not monounsaturated" $
      isMonounsaturated tg_181_181_181
      `shouldBe` False

  describe "isPolyunsaturated" $ do
    it "TG 18:0_18:0_18:0 is polyunsaturated" $
      isPolyunsaturated tg_180_180_180
      `shouldBe` False
    it "TG 18:0_18:0_18:1 is not polyunsaturated" $
      isPolyunsaturated tg_180_180_181
      `shouldBe` False
    it "TG 18:1_18:1_18:1 is not polyunsaturated" $
      isPolyunsaturated tg_181_181_181
      `shouldBe` True
