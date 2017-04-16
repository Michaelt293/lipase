module Main where

import Isotope
import Spectra
import FattyAcid
import Triacylglycerol
import System.Environment
import Data.List
import Data.Char
import Data.Maybe
import Data.Monoid
import System.Directory
import Data.List.Split

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
  let msSpectrum = mkMSSpectrum unprocessedMSSpectrum
      ms2Spectra = sort $ filterByRelativeAbundance (> 5) <$>
                          zipWith
                            mkMS2SpectrumRemovePrecursor
                            precursorIonsMz
                            unprocessedMS2Spectra
      neutralLossSpectra = calNeutralLosses <$> ms2Spectra
      fasFromNeutralLoss = assignFAsFromNeutralLoss <$> neutralLossSpectra
      tagsFromAssignedFas = findPossibleTAGs msSpectrum <$> fasFromNeutralLoss
      finalResults = toFinalResults tagsFromAssignedFas
      -- Results
      identifiedFAs' = identifiedFAs finalResults
      identifiedTags' = identifiedTags finalResults
      allTagFAs' = allTagFAs finalResults
      identifiedCondensedTags' = identifiedCondensedTags finalResults
  putStrLn "Total tentatively assigned fatty acids with normalised abundances:"
  putStrLn . intercalate "\n" $ renderPairNormalisedAbundance <$> identifiedFAs'
  putStrLn "Total assigned triacylglycerols:"
  putStrLn . unwords $ show <$> identifiedTags'
  putStrLn "Total triacylglycerol fatty acids:"
  putStrLn . unwords $ show <$> allTagFAs'
  putStrLn "Total condensed triacylglycerols with normalised abundances:"
  putStrLn . intercalate "\n" $ renderPairNormalisedAbundance <$> identifiedCondensedTags'
  --    finalResult' = finalResults finalResult
  -- putStrLn "Fatty acids in the search list"
  -- print $ sort fattyAcyls
  -- --putStrLn "Total tentatively assigned fatty acids"
  -- --print $ allTentativelyAssignedFAs finalResult
  -- mapM_ putStrLn $ tagMzNormalisedAbundances finalResult'
  -- --print $ sumShouldEqual1 finalResult'
  -- mapM_ (putStrLn . renderFattyAcylNormalisedAbundance) $ accumulateNormalisedAbundance finalResult'


  --putStrLn "Total assigned triacylglycerols"
  --print $ allAssignedTAGs finalResult
  --putStrLn "Total triacylglycerol fatty acids"
  --print $ allTagFAs finalResult
  --putStrLn . showList' $ formatNormalisedAbundanceFAsIndentified <$> finalResult
