\name{evaluate}
\alias{evaluate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 biomod2 modeling outputs evaluation
%%  ~~function to do ... ~~
}
\description{
 This function will evaluate biomod2 modeling output for given metrics (e.g 'TSS', 'ROC'...) for a given dataset.
}
\usage{
evaluate(model, data, stat, as.array=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ the model you want evaluate (either \code{"BIOMOD.models.out"}, \code{"BIOMOD.EnsembleModeling.out"} or \code{"biomod2_model"} )}
  \item{data}{ the \code{data.set} on which you want to perform analyses. Must be a dataset with first column containing the observed data for your species. The folowing columns must be the explanatory variables at observed points. Be sure that colnames of your dataset are the name of your species then the names of variables used for building models at previous steps.}
  \item{stat}{vector of statistic metrics names (e.g 'TSS','ROC') you want to perform. (see \code{\link[biomod2]{BIOMOD_Modeling}}) to get the list ofall available metrics) }
  \item{as.array}{logical, (FALSE by default) if FALSE a list of avaluation tables is returned (one item by models). If TRUE, the output will be return under 'classical' biomod2 array objects (see \code{\link[biomod2]{BIOMOD_Projection}})}
}

\details{
It's more or less base on the same principle than \code{\link[randomForest]{randomForest}} variables importance algorithm. The principle is to shuffle a single variable of the given data. Make model prediction with this 'shuffled' data.set. Then we compute a simple correlation (Pearson's by default) between references predictions and the 'shuffled' one. The return score is 1-cor(pred_ref,pred_shuffled). The highest the value, the more influence the variable has on the model. A value of this 0 assumes no influence of that variable on the model. Note that this technique does not account for interactions between the variables.
}

\value{
  a \code{list} or an \code{array} containing for each evaluation metric the score, the threshold considered to transform continuous data into binary ones (for all metrics exepted 'ROC') and associated sensibility and specificity. 
  
}

\author{
  Damien Georges
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{BIOMOD_EnsembleModeling}}, \code{\link[biomod2]{variables_importance}}
}

\examples{


# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                     models = c('SRE','CTA','RF'), 
                                     models.options = myBiomodOption, 
                                     NbRunEval=1, 
                                     DataSplit=80, 
                                     Yweights=NULL, 
                                     VarImport=3, 
                                     models.eval.meth = c('TSS'),
                                     SaveObj = TRUE,
                                     rescal.all.models = FALSE,
                                     do.full.models = FALSE,
                                     modeling.id='test')
                                     
# 4. Evaluate model over another dataset (here the full one)

## creation of suitable dataset
data <- cbind(GuloGulo=get_formal_data(myBiomodModelOut,'resp.var'), 
              get_formal_data(myBiomodModelOut,'expl.var'))

## evaluation
evaluate(myBiomodModelOut, data=data, stat=c('ROC','TSS'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ evaluation }
\keyword{ models }
\keyword{ score }% __ONLY ONE__ keyword per line
