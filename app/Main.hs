module Main where

import Spectra
import System.Environment
import Data.List
import Data.Char
import System.Directory
import Data.List.Split

main :: IO ()
main = do
  (dir:_) <- getArgs
  csvFiles <- filter (isSuffixOf ".csv") <$> getDirectoryContents dir
  let msSpectrumFile =
        case filter (isInfixOf "pos_MS_spectrum") csvFiles of
          [spec] -> spec
          _ : _  -> error "more than one MS spectrum"
          _      -> error "no MS spectrum"
  let ms2SpectraFiles = filter (not . isInfixOf "pos_MS_spectrum") csvFiles
  let precursorIonsMz = findPrecursorIonMz <$> ms2SpectraFiles
        where findPrecursorIonMz fileName = Mz $
                (read . takeWhile isDigit $
                 dropWhile (not . isDigit) fileName) + 0.5 -- maybe don't hard code this value
  let readCsv = readFile . (\x -> dir ++ "/" ++ x)
  unprocessMsSpectrum <- readCsv msSpectrumFile
  unprocessMs2Spectrum <- traverse readCsv ms2SpectraFiles
  let processCsv csv =  splitOn "," <$> drop 2 (lines csv)
  let msSpectrum = MSSpectrum . insertAbundances . readSpectrum . processCsv $ unprocessMsSpectrum
  let ms2Spectra = zipWith
                     MS2Spectrum
                     precursorIonsMz
                     (insertAbundances . readSpectrum . processCsv <$> unprocessMs2Spectrum)
  let filteredMS2Spectra = removePrecursorIon . filterByRelativeAbundance (RelativeAbundance 5) <$> ms2Spectra
  let neutralLossSpectra = toNeutralLossSpectrum <$> filteredMS2Spectra
  print neutralLossSpectra
