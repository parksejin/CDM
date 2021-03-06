\name{IRT.modelfit}
\alias{IRT.modelfit}
\alias{IRT.modelfit.din}
\alias{IRT.modelfit.gdina}
\alias{summary.IRT.modelfit.din}
\alias{summary.IRT.modelfit.gdina}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
S3 Methods for Assessing Model Fit
}

\description{
This S3 method asesses global (absolute) model fit using
the methods described in \code{\link{modelfit.cor.din}}.
}

\usage{
IRT.modelfit(object, \dots)

\method{IRT.modelfit}{din}(object, \dots)
\method{IRT.modelfit}{gdina}(object, \dots)

\method{summary}{IRT.modelfit.din}(object, \dots)
\method{summary}{IRT.modelfit.gdina}(object, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of classes \code{\link{din}} or \code{\link{gdina}}.
}
  \item{\dots}{
More arguments to be passed.
}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
See output of \code{\link{modelfit.cor.din}}.
}


%\references{
%% ~put references to the literature/web site here ~
%}


%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For extracting the individual likelihood or posterior see
\code{\link{IRT.likelihood}} or \code{\link{IRT.posterior}}.

The model fit of objects of class \code{\link{gdm}} can be obtained
by using the \code{tam.modelfit.IRT} function in the \pkg{TAM} package.
}
%%% TAM, sirt and mirt links

\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: Absolute model fit
#############################################################################
data(sim.dina)
data(sim.qmatrix)

#*** Model 1: DINA model for DINA simulated data
mod1 <- din(sim.dina, q.matr = sim.qmatrix, rule = "DINA" )
fmod1 <- IRT.modelfit( mod1 )
summary(fmod1)

#*** Model 2: GDINA model
mod2 <- gdina( sim.dina , q.matr = sim.qmatrix , rule="GDINA" )
fmod2 <- IRT.modelfit( mod2 )
summary(fmod2)	
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{Individual likelihood}
% \keyword{Individual posterior}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
