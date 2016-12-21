
module Triacylglycerol where

import FattyAcid
import Spectra
import Isotope
import Data.List
import Control.Monad

data Triacylglycerol = Triacylglycerol {
    fa1 :: FattyAcyl
  , fa2 :: FattyAcyl
  , fa3 :: FattyAcyl
}

instance Eq Triacylglycerol where
  Triacylglycerol a1 b1 c1 == Triacylglycerol a2 b2 c2 =
    sort [a1, b1, c1] == sort [a2, b2, c2]

instance Ord Triacylglycerol where
  Triacylglycerol a1 b1 c1 <= Triacylglycerol a2 b2 c2 =
    sort [a1, b1, c1] <= sort [a2, b2, c2]

instance Show Triacylglycerol where
  show (Triacylglycerol a b c) =
    "TG(" ++ intercalate "_" (show <$> sort [a, b, c]) ++ ")"

instance ToMolecularFormula Triacylglycerol where
  toMolecularFormula (Triacylglycerol a b c) =
    mkMolecularFormula [(C, 3), (H, 5)] |+| foldMap toMolecularFormula [a, b, c]

instance ToElementalComposition Triacylglycerol where
  toElementalComposition = toElementalComposition . toMolecularFormula

-- Brute force function to find possible TAGs. Write a more efficient
-- implementation if there are performance issues.
possibleTAGs :: Double -> MonoisotopicMass -> [FattyAcyl] -> [Triacylglycerol]
possibleTAGs n prec fas = nub $ do
  fa1 <- fas
  fa2 <- fas
  fa3 <- fas
  guard (withinTolerance n prec (monoisotopicMass (Triacylglycerol fa1 fa2 fa3)))
  return $ Triacylglycerol fa1 fa2 fa3
