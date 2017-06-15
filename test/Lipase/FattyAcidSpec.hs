module Lipase.FattyAcidSpec (spec) where

import Isotope
import FattyAcid
import Test.Hspec

spec :: Spec
spec = do
  describe "Elemental composition for FattyAcyl" $ do
    it "Stearic acid acyl group elemental composition should be C18H35O2" $
      toElementalComposition (FattyAcyl 18 0)
      `shouldBe` mkElementalComposition [(C, 18), (H, 35), (O, 2)]
    it "Oleic acid acyl group elemental composition should be C18H33O2" $
      toElementalComposition (FattyAcyl 18 1)
      `shouldBe` mkElementalComposition [(C, 18),(H, 33), (O, 2)]
    it "DHA elemental acyl group composition should be C22H31O2" $
      toElementalComposition (FattyAcyl 22 6)
      `shouldBe` mkElementalComposition [(C, 22),(H, 31), (O, 2)]

  describe "Elemental composition for FattyAcid" $ do
    it "Stearic acid elemental composition should be C18H36O2" $
      toElementalComposition (FattyAcid (FattyAcyl 18 0))
      `shouldBe` mkElementalComposition [(C, 18), (H, 36), (O, 2)]
    it "Oleic acid elemental composition should be C18H34O2" $
      toElementalComposition (FattyAcid (FattyAcyl 18 1))
      `shouldBe` mkElementalComposition [(C, 18),(H, 34), (O, 2)]
    it "DHA elemental composition should be C22H32O2" $
      toElementalComposition (FattyAcid (FattyAcyl 22 6))
      `shouldBe` mkElementalComposition [(C, 22),(H, 32), (O, 2)]

  describe "assignFA" $ do
    it "284.3 should be stearic acid" $
      assignFA 0.3 (MonoisotopicMass 284.3) fattyAcids
      `shouldBe` Just (FattyAcid (FattyAcyl 18 0))
    it "282.3 should be oleic acid" $
      assignFA 0.3 (MonoisotopicMass 282.3) fattyAcids
      `shouldBe` Just (FattyAcid (FattyAcyl 18 1))
    it "328.2 should be DHA" $
      assignFA 0.3 (MonoisotopicMass 328.2) fattyAcids
      `shouldBe` Just (FattyAcid (FattyAcyl 22 6))
