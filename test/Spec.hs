module Main where

import Test.Hspec
import Lipase.FattyAcidSpec
import Lipase.SpectraSpec
import Lipase.TriacylglycerolSpec

main :: IO ()
main = hspec $ do
  describe "FattyAcid" Lipase.FattyAcidSpec.spec
  describe "Spectra" Lipase.SpectraSpec.spec
  describe "Triacylglycerol" Lipase.TriacylglycerolSpec.spec
