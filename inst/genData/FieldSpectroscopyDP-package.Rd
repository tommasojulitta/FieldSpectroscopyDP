\name{FieldSpectroscopyDP-package}

\alias{FieldSpectroscopyDP-package}
\alias{FieldSpectroscopyDP}

\title{Field Spectroscopy Data Processing}

\description{
Data processing of spectral measurements collected in the field. Conversion of raw data acquired in digital counts in physical unit (Radiance). Calculation of reflectance and reflectance based indices. Spectral convolution. Retrieval of Sun Induced Chlorophyll Fluorescence.
}

\details{

Basic Functions are 
\itemize{
\item{ Convert Digital number to Radiance: \code{\link{GetRadiance}}  }
\item{ Compute the Reflectance Factor from upwelling and downwelling radiance or counts: \code{\link{GetReflectance}}  }
\item{ Compute the linear interpolation of each wavelength of two measuremnts at a specific time occurred in between the two: \code{\link{SigTimeInterpolation}}  }
\item{ Compute the mathematic combination of n bands according to a specified formula and to a specified band convlution (average or gaussian): \code{\link{BandMath}}  }

}%itemize

Fluorescence Functions are 
\itemize{
\item{ Compute the single Fraunhofer Line Discriminator for fluorescence retrieval: \code{\link{sFLD}}  }
\item{ Compute the 3 Fraunhofer Line Discriminator for fluorescence retrieval: \code{\link{FLD3}}  }
\item{ Compute the improved Fraunhofer Line Discriminator for fluorescence retrieval: \code{\link{iFLD}}  }
\item{ Compute spectral fitting methods for fluorescence retrieval: \code{\link{SFM}}  }
}%itemize

Convolution Functions are 
\itemize{
\item{ Compute the statistics (mean, median, and standard deviation) of a spectrum given a certain spectral range: \code{\link{fConvol}}  }
\item{ Compute the gaussian convolution on a specific spectral range: \code{\link{fGaussianConvol}}  }
\item{ Compute the gaussian convolution on a specific spectral range: \code{\link{fWeightAvgConvol}}  }
\item{ Compute the standard deviationof gaussian convolution on a specific spectral range: \code{\link{fGaussianConvolSD}}  }
\item{ Compute the gaussian convolution on multiple spectral range: \code{\link{fGaussianConvolMultiWl}}  }
\item{ Resample based on gaussian convolution: \code{\link{fSpecRes}}  }
\item{ Write the Rdata object with the configuration of the satellites: \code{\link{fCreateConfSatellite}}  }
}%itemize

Input Output functions are 
\itemize{
\item{ Read Spectra from ASCII files and list them : \code{\link{ReadSpectra}}  }
\item{ Read Integration time from ASCII files acquired by Ocean View software: \code{\link{ReadSpecITOceanView}}  }
\item{ Read Data from csv files saved by FloX and list them: \code{\link{ReadFlOX}}  }
\item{ Read Data from csv files saved by DFloX and list them: \code{\link{ReadDFlOX}}  }
}%itemize

Further misc functionality: 
\itemize{
\item{ Convert date to julian day of the year and day fraction: \code{\link{DateToDOY}}  }
\item{ Convert FloX date format in as.POSIXct class: \code{\link{DateTimeFloX}}  }
\item{ Compute the spectral sum of reflected radiance and fluorescence: \code{\link{SpecFit}}  }
\item{ Calculate the first guess of fluorescence and true reflectance based on gaussian function and spline fitting on apparent reflectance: \code{\link{FirstGuess}}  }
\item{ Calculate the optimized fluorescence and true reflectance: \code{\link{SFMResults}}  }
\item{ Round to the next order fo magnitude: \code{\link{log10_ceiling}}  }
}%itemize

}%details

\author{Thomas Wutzler}



\keyword{ package }





