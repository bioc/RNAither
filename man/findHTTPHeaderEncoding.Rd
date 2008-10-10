\name{findHTTPHeaderEncoding}
\alias{findHTTPHeaderEncoding}

\title{ Find the encoding of the HTTP response from the HTTP header }
\description{
  This function is currently made available so it can be called from C code to find the charset from the HTTP header in the response from an HTTP server. It maps this charset value to one of the known R encodings (UTF-8 or LATIN1) or returns the native encoding.
}
\usage{
findHTTPHeaderEncoding(str)
}

\arguments{
  \item{str}{ one or more lines from the HTTP header }
}

\value{
  NA or an integer value indicating the encoding to be used. This integer corresponds to the cetype\_t enumeration in Rinternals.h.
}

\references{ Writing R Extensions Manual and the section(s) on character encodings }

\examples{
findHTTPHeaderEncoding("Content-Type: text/html;charset=ISO-8859-1\r\n")
findHTTPHeaderEncoding("Content-Type: text/html; charset=utf-8\r\n")

}

\keyword{character}
