{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Triacylglycerol where

import Isotope
import Data.List

newtype NumCarbons = NumCarbons { getNumCarbons :: Int }
  deriving (Eq, Ord, Num)

instance Show NumCarbons where
  show (NumCarbons n) = show n

newtype NumDoubleBonds = NumDoubleBonds { getNumDoubleBonds :: Int }
  deriving (Eq, Ord, Num)

instance Show NumDoubleBonds where
  show (NumDoubleBonds n) = show n

data FattyAcyl = FattyAcyl {
    numCarbons     :: NumCarbons
  , numDoubleBonds :: NumDoubleBonds
} deriving (Eq, Ord)

instance Show FattyAcyl where
  show (FattyAcyl cs dbs) = show cs ++ ":" ++ show dbs

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

instance ToMolecularFormula FattyAcyl where
  toMolecularFormula (FattyAcyl cs dbs) = mkMolecularFormula
    [ (C, getNumCarbons cs)
    , (H, getNumCarbons cs * 2 + 1 - getNumDoubleBonds dbs * 2)
    , (O, 2)
    ]

instance ToElementalComposition FattyAcyl where
  toElementalComposition = toElementalComposition . toMolecularFormula

instance ToMolecularFormula Triacylglycerol where
  toMolecularFormula (Triacylglycerol a b c) =
    mkMolecularFormula [(C, 3), (H, 5)] |+| foldMap toMolecularFormula [a, b, c]

instance ToElementalComposition Triacylglycerol where
  toElementalComposition = toElementalComposition . toMolecularFormula
