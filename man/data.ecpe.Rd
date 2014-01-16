\name{data.ecpe}
\alias{data.ecpe}
\docType{data}
\title{
Dataset ECPE
}
\description{
ECPE dataset from the Templin and Hoffman (2013) tutorial of
specifying cognitive diagnostic models in Mplus.
}

\usage{
data(data.ecpe)
}

\format{
The format of the data is a list containing the dichotomous item
response data \code{data} (2922 persons at 28 items) 
and the Q-matrix \code{q.matrix} (28 items and 3 skills):

\code{List of 2} \cr
\code{ $ data    :'data.frame':} \cr
\code{  ..$ id : int [1:2922] 1 2 3 4 5 6 7 8 9 10 ...} \cr
\code{  ..$ E1 : int [1:2922] 1 1 1 1 1 1 1 0 1 1 ...} \cr
\code{  ..$ E2 : int [1:2922] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ E3 : int [1:2922] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  ..$ E4 : int [1:2922] 0 1 1 1 1 1 1 1 1 1 ...} \cr
\code{  [...] } \cr
\code{  ..$ E27: int [1:2922] 1 1 1 1 1 1 1 0 1 1 ...} \cr
\code{  ..$ E28: int [1:2922] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ q.matrix:'data.frame':} \cr
\code{  ..$ skill1: int [1:28] 1 0 1 0 0 0 1 0 0 1 ...} \cr
\code{  ..$ skill2: int [1:28] 1 1 0 0 0 0 0 1 0 0 ...} \cr
\code{  ..$ skill3: int [1:28] 0 0 1 1 1 1 1 0 1 0 ...} \cr

The skills are

\code{skill1}: Morphosyntactic rules

\code{skill2}: Cohesive rules

\code{skill3}: Lexical rules.

}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
The dataset is used in Templin and Hoffman (2013).

The dataset was downloaded from
\url{http://psych.unl.edu/jtemplin/teaching/dcm/dcm12ncme/}.
}


\references{
Templin, J., & Hoffman, L. (2013). 
Obtaining diagnostic classification model estimates using Mplus. 
\emph{Educational Measurement: Issues and Practice}, \bold{32}, 37-50.
}

\section{Example Index}{
\code{\link{sim.din}} (Example 3),
\code{\link{sim.gdina}} (Examples 2, 3)
}

\examples{
\dontrun{
data(data.ecpe)

# Model 1: LCDM model 
mod1 <- gdina( data.ecpe$data[,-1], q.matrix= data.ecpe$q.matrix , link="logit")
summary(mod1)

# Model 2: DINA model 
mod2 <- gdina( data.ecpe$data[,-1], q.matrix= data.ecpe$q.matrix , rule="DINA")
summary(mod2)

# Model comparison using likelihood ratio test
anova(mod1,mod2)
##       Model   loglike Deviance Npars      AIC      BIC    Chisq df  p
##   2 Model 2 -42841.61 85683.23    63 85809.23 86185.97 206.0359 18  0
##   1 Model 1 -42738.60 85477.19    81 85639.19 86123.57       NA NA NA
}
}
\keyword{datasets}