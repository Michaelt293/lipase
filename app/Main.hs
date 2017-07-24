module Main where

import Spectra
import FattyAcid
import Triacylglycerol
import Results
import Report
import Isotope.Ion (Mz(..))
import System.Environment (getArgs)
import System.Directory (getDirectoryContents)
import System.Exit (die)
import Data.Foldable (traverse_)
import Data.List (isSuffixOf, isInfixOf, sort)
import Data.Char (isDigit)
import Data.Monoid ((<>))
import Data.Csv (encodeDefaultOrderedByName)
import qualified Data.ByteString.Lazy as ByteString
import qualified Data.ByteString as ByteStringStrict
import Control.Monad (void)


findPrecursorIonMz :: String -> Mz
findPrecursorIonMz fileName =
      (Mz . read . takeWhile isDigit $
       dropWhile (not . isDigit) fileName) + 0.5 -- maybe don't hard code 0.5

main :: IO ()
main = do
  (dir:_) <- getArgs
  csvFiles <- filter (isSuffixOf ".csv") <$> getDirectoryContents dir
  msSpectrumFile <-
        case filter (isInfixOf "pos_MS_spectrum") csvFiles of
          [spec] -> pure spec
          _ : _  -> die "more than one MS spectrum"
          _      -> die "no MS spectrum"
  let ms2SpectraFiles = filter (not . isInfixOf "pos_MS_spectrum") csvFiles
      precursorIonsMz = findPrecursorIonMz <$> ms2SpectraFiles
      readCsv = process . (\csvFile -> dir <> "/" <> csvFile)
  unprocessedMSSpectrum <- readCsv msSpectrumFile
  unprocessedMS2Spectra <- traverse readCsv ms2SpectraFiles
  let msSpectrum' = mkMSSpectrum unprocessedMSSpectrum
      ms2Spectra = sort $ filterByRelativeAbundance (> 5) <$>
                          zipWith
                            mkMS2SpectrumRemovePrecursor
                            precursorIonsMz
                            unprocessedMS2Spectra
      neutralLossSpectra  = calNeutralLosses <$> ms2Spectra
      fasFromNeutralLoss  = assignFAsFromNeutralLoss <$> neutralLossSpectra
      tagsFromAssignedFas = findPossibleTAGs msSpectrum' <$> fasFromNeutralLoss
      finalResults        = toFinalResults tagsFromAssignedFas
      -- Results
      (ionsIdentified, identifiedFAs') = identifiedFAs finalResults
      tagSummary'                      = tagSummary finalResults
  putStrLn "Total tentatively assigned fatty acids with normalised abundances:"
  traverse_ print identifiedFAs'
  print ionsIdentified
  void $ ByteString.writeFile "identifiedFAs.csv" (encodeDefaultOrderedByName identifiedFAs')
  putStrLn "Total assigned triacylglycerols:"
  traverse_ print tagSummary'
  void $ ByteString.writeFile "tagSummary.csv" (encodeDefaultOrderedByName tagSummary')
  void $ ByteStringStrict.writeFile "report.txt" (report ionsIdentified identifiedFAs' tagSummary')
