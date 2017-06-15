module Main where

import Spectra
import FattyAcid
import Triacylglycerol
import Isotope.Ion (Mz(..))
import System.Environment (getArgs)
import System.Directory (getDirectoryContents)
import Data.Foldable (traverse_)
import Data.List (isSuffixOf, isInfixOf, sort)
import Data.Char (isDigit)
import Data.Monoid ((<>))

findPrecursorIonMz :: String -> Mz
findPrecursorIonMz fileName =
      (Mz . read . takeWhile isDigit $
       dropWhile (not . isDigit) fileName) + 0.5 -- maybe don't hard code 0.5

main :: IO ()
main = do
  (dir:_) <- getArgs
  csvFiles <- filter (isSuffixOf ".csv") <$> getDirectoryContents dir
  let msSpectrumFile =
        case filter (isInfixOf "pos_MS_spectrum") csvFiles of
          [spec] -> spec
          _ : _  -> error "more than one MS spectrum"
          _      -> error "no MS spectrum"
      ms2SpectraFiles = filter (not . isInfixOf "pos_MS_spectrum") csvFiles
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
      neutralLossSpectra = calNeutralLosses <$> ms2Spectra
      fasFromNeutralLoss = assignFAsFromNeutralLoss <$> neutralLossSpectra
      tagsFromAssignedFas = findPossibleTAGs msSpectrum' <$> fasFromNeutralLoss
      finalResults = toFinalResults tagsFromAssignedFas
      -- Results
      identifiedFAs' = identifiedFAs finalResults
      identifiedTags' = identifiedTAGSummary finalResults
      allTagFAs' = allTagFAs finalResults
      identifiedCondensedTags' = identifiedCondensedTags finalResults
  putStrLn "Total tentatively assigned fatty acids with normalised abundances:"
  traverse_ putStrLn $ renderPairNormalisedAbundance <$> identifiedFAs'
  putStrLn "Total assigned triacylglycerols:"
  traverse_ putStrLn $
    (\(x, y) -> unwords (show <$> x) <> " : " <> unwords (show <$> y)) <$>
    sort identifiedTags'
  putStrLn "Total triacylglycerol fatty acids:"
  putStrLn . unwords $ show <$> allTagFAs'
  putStrLn "Total condensed triacylglycerols with normalised abundances:"
  traverse_ putStrLn $ renderPairNormalisedAbundance <$> identifiedCondensedTags'
