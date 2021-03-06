% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Convert_MICs.R
\name{Convert_MICs}
\alias{Convert_MICs}
\title{Convert_MICs - Processing raw AST data}
\usage{
Convert_MICs(x, drugs, mic_text = NULL, sign_text = "SIGN",
  append_orig = F, sep = "_", mic_format = NULL)
}
\arguments{
\item{x}{A matrix-like object containing MIC data.}

\item{drugs}{A character vector of column names corresponding to the drug names/abbreviations to be converted.}

\item{mic_text}{The character string used to identify columns containing MICs.}

\item{sign_text}{The character string used to identify columns containing the comparison signs. Ignored where mic_format = 'Standard'.}

\item{append_orig}{If true, formatted columns are appended to original data. If false, only}

\item{sep}{The character used to seperate the drug names from the mic_text and sign_text in column names.}

\item{mic_format}{A character string indicating the format of the data. Currently accepts}
}
\value{
a data.frame containing numeric MIC data, and metadata if selected.
}
\description{
Antimicrobial susceptibility testing (AST) data from surveillance programs is typically determined via serial microdilution and reported as interval censored minimum inhibitory concentrations (MIC). This function converts MIC values reported as a combination of a comparison sign characters (e.g., > or <=) and numerical characters (e.g., 0.5 or 8) to numerical values prior to analysis. For example, an entry of '= 2' indicates the bacteria was inhibited by 2 ug/mL of drug, but was not inhibited by 1 ug/mL. In most cases, the comparison character can be ignored, but entries with '>' indicate the value exceeds the highest tested drug concentration. To represent the increased resistance of these isolates, the reported concentration must be adjusted upwards by one concentration, even though that concentration was not tested. For example the entry '> 4' (growth was not inhibited at 4 ug/mL, the highest tested concentration) should be adjusted to 8 to represent these results indicated greater resistance than '= 4' (growth was inhibited at 4 ug/mL, but not at 2 ug/mL). This function removes the comparison signs and adjusts the values as neccessary, returning a more conveniently labeled data.frame containing numeric values and selected metadata. Two formats are currently supported: 1) The standard format where the sign is held in the same field as the concentration, and 2) the format used by the NARMS program where comparison sign is stored in a seperate field.
}
