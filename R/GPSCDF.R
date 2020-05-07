#'Generalized Propensity Score Cumulative Distribution Function (GPS-CDF)
#'
#'\code{GPSCDF} takes in a generalized propensity score (GPS) object with length
#'>2 and returns the GPS-CDF balancing score.
#'
#'The \code{GPSCDF} method is used to conduct propensity score matching and
#'stratification for both ordinal and multinomial treatments. The method
#'directly maps any GPS vector (with length >2) to a single scalar value that
#'can be used to produce either average treatment effect (ATE) or average
#'treatment effect among the treated (ATT) estimates. For the \code{K}
#'multinomial treatments setting, the balance achieved from each \code{K!}
#'ordering of the GPS should be assessed to find the optimal ordering of the GPS
#'vector (see Examples for more details).
#'
#'@importFrom dplyr ntile
#'@importFrom nbpMatching distancematrix nonbimatch
#'@importFrom nnet multinom
#'@importFrom MASS mvrnorm
#'@importFrom stats nls nls.control sd
#'@importFrom survival clogit coxph strata
#'@importFrom utils capture.output
#'
#'@param pscores The object containing the treatment ordered generalized
#'  propensity scores for each subject.
#'@param data An optional data frame to attach the calculated balancing score.
#'  The data frame will also be used in stratification and matching.
#'@param trt An optional object containing the treatment variable.
#'@param stratify Option to produce strata based on the power parameter
#'  (\code{ppar}). Default is \code{FALSE}.
#'@param nstrat An optional parameter for the number of strata to be created
#'  when \code{stratify} is set to \code{TRUE}. Default is \code{5} strata.
#'@param optimal Option to perform optimal matching of subjects based on the
#'  power parameter (\code{ppar}). Default is \code{FALSE}.
#'@param greedy Option to perform greedy matching of subjects based on the power
#'  parameter (\code{ppar}). Default is \code{FALSE}.
#'@param ordinal Specifies ordinal treatment groups for matching. Subjects are
#'  matched based on the ratio of the squared difference of power parameters for
#'  two subjects, \code{ppar_i} and \code{ppar_j}, in the numerator and the
#'  squared difference in observed treatment received, \code{trt_i} and
#'  \code{trt_j}, in the denominator: \code{(ppar_i-ppar_j)^2/(trt_i-trt_j)^2}.
#'  Default is \code{FALSE}.
#'@param multinomial Specifies multinomial treatment groups for matching.
#'  Subjects are matched based on the absolute difference of power parameters
#'  for two subjects, \code{ppar_i} and \code{ppar_j}, who received different
#'  treatments: \code{|ppar_i - ppar_j|}. Default is \code{FALSE}.
#'@param caliper An optional parameter for the caliper value used when
#'  performing greedy matching. Used when \code{greedy} is set to \code{TRUE}.
#'  Default is \code{.25*sd(ppar)}.
#'@return \item{ppar}{The power parameter scalar balancing score to be used in
#'  outcome analyses through stratification or matching.} \item{data}{The user
#'  defined dataset with power parameter (ppar), strata, and/or optimal matching
#'  variables attached.} \item{nstrat}{The number of strata used for
#'  stratification.} \item{strata}{The strata produced based on the calculated
#'  power parameter (\code{ppar}).} \item{optmatch}{The optimal matches produced
#'  based on the calculated power parameter (\code{ppar}).}
#'  \item{optdistance}{The average absolute total distance of power parameters
#'  (\code{ppars}) for optimally matched pairs.} \item{caliper}{The caliper
#'  value used for greedy matching.} \item{grddata}{The user defined dataset
#'  with greedy matching variable attached.} \item{grdmatch}{The greedy matches
#'  produced based on the calculated power parameter (\code{ppar}).}
#'  \item{grdydistance}{The average absolute total distance of power parameters
#'  (\code{ppars}) for greedy matched pairs.}
#'@author Derek W. Brown, Thomas J. Greene, Stacia M. DeSantis
#'@references Greene, TJ. (2017). Utilizing Propensity Score Methods for Ordinal
#'  Treatments and Prehospital Trauma Studies. Texas Medical Center
#'  Dissertations (via ProQuest).
#' @examples
#'
#'
#'### Example: Create data example
#' N<- 100
#'
#' set.seed(18201) # make sure data is repeatable
#' Sigma <- matrix(.2,4,4)
#' diag(Sigma) <- 1
#' data<-matrix(0, nrow=N, ncol=6,dimnames=list(c(1:N),
#'       c("Y","trt",paste("X",c(1:4),sep=""))))
#' data[,3:6]<-matrix(MASS::mvrnorm(N, mu=rep(0, 4), Sigma,
#'       empirical = FALSE) , nrow=N, ncol = 4)
#'
#' dat<-as.data.frame(data)
#'
#'
#' #Create Treatment Variable
#' tlogits<-matrix(0,nrow=N,ncol=2)
#' tprobs<-matrix(0,nrow=N,ncol=3)
#'
#' alphas<-c(0.25, 0.3)
#' strongbetas<-c(0.7, 0.4)
#' modbetas<-c(0.2, 0.3)
#'
#' for(j in 1:2){
#'   tlogits[,j]<- alphas[j] + strongbetas[j]*dat$X1 + strongbetas[j]*dat$X2+
#'                 modbetas[j]*dat$X3 + modbetas[j]*dat$X4
#' }
#'
#' for(j in 1:2){
#'   tprobs[,j]<- exp(tlogits[,j])/(1 + exp(tlogits[,1]) + exp(tlogits[,2]))
#'   tprobs[,3]<- 1/(1 + exp(tlogits[,1]) + exp(tlogits[,2]))
#' }
#'
#' set.seed(91187)
#' for(j in 1:N){
#'   data[j,2]<-sample(c(1:3),size=1,prob=tprobs[j,])
#' }
#'
#'
#' #Create Outcome Variable
#' ylogits<-matrix(0,nrow=N,ncol=1,dimnames=list(c(1:N),c("Logit(P(Y=1))")))
#' yprobs<-matrix(0,nrow=N,ncol=2,dimnames=list(c(1:N),c("P(Y=0)","P(Y=1)")))
#'
#' for(j in 1:N){
#'   ylogits[j,1]<- -1.1 + 0.7*data[j,2] + 0.6*dat$X1[j] + 0.6*dat$X2[j] +
#'                  0.4*dat$X3[j] + 0.4*dat$X4[j]
#'
#'   yprobs[j,2]<- 1/(1+exp(-ylogits[j,1]))
#'
#'   yprobs[j,1]<- 1-yprobs[j,2]
#' }
#'
#' set.seed(91187)
#' for(j in 1:N){
#'   data[j,1]<-sample(c(0,1),size=1,prob=yprobs[j,])
#' }
#'
#' dat<-as.data.frame(data)
#'
#'
#' ### Example: Using GPSCDF
#'
#' #Create the generalized propensity score (GPS) vector using any parametric or
#' #nonparametric model
#'
#' glm<- nnet::multinom(as.factor(trt)~ X1+ X2+ X3+ X4, data=dat)
#' probab<- round(predict(glm, newdata=dat, type="probs"),digits=8)
#' gps<-cbind(probab[,1],probab[,2],1-probab[,1]-probab[,2])
#'
#'
#' #Create scalar balancing power parameter
#' fit<-GPSCDF(pscores=gps)
#'
#' \dontrun{
#'   fit$ppar
#' }
#'
#'
#' #Attach scalar balancing power parameter to user defined data set
#' fit2<-GPSCDF(pscores=gps, data=dat)
#'
#' \dontrun{
#'   fit2$ppar
#'   fit2$data
#' }
#'
#'
#' ### Example: Ordinal Treatment
#'
#' #Stratification
#' fit3<-GPSCDF(pscores=gps, data=dat, stratify=TRUE, nstrat=5)
#'
#' \dontrun{
#'   fit3$ppar
#'   fit3$data
#'   fit3$nstrat
#'   fit3$strata
#'
#'   library(survival)
#'   model1<-survival::clogit(Y~as.factor(trt)+X1+X2+X3+X4+strata(strata),
#'                            data=fit3$data)
#'   summary(model1)
#' }
#'
#'
#' #Optimal Matching
#' fit4<- GPSCDF(pscores=gps, data=dat, trt=dat$trt, optimal=TRUE, ordinal=TRUE)
#'
#' \dontrun{
#'   fit4$ppar
#'   fit4$data
#'   fit4$optmatch
#'   fit4$optdistance
#'
#'   library(survival)
#'   model2<-survival::clogit(Y~as.factor(trt)+X1+X2+X3+X4+strata(optmatch),
#'                            data=fit4$data)
#'   summary(model2)
#' }
#'
#'
#' #Greedy Matching
#' fit5<- GPSCDF(pscores=gps, data=dat, trt=dat$trt, greedy=TRUE, ordinal=TRUE)
#'
#' \dontrun{
#'   fit5$ppar
#'   fit5$data
#'   fit5$caliper
#'   fit5$grddata
#'   fit5$grdmatch
#'   fit5$grdydistance
#'
#'   library(survival)
#'   model3<-survival::clogit(Y~as.factor(trt)+X1+X2+X3+X4+strata(grdmatch),
#'                            data=fit5$grddata)
#'   summary(model3)
#' }
#'
#'
#' ### Example: Multinomial Treatment
#'
#' #Create all K! orderings of the GPS vector
#' gps1<-cbind(gps[,1],gps[,2],gps[,3])
#' gps2<-cbind(gps[,1],gps[,3],gps[,2])
#' gps3<-cbind(gps[,2],gps[,1],gps[,3])
#' gps4<-cbind(gps[,2],gps[,3],gps[,1])
#' gps5<-cbind(gps[,3],gps[,1],gps[,2])
#' gps6<-cbind(gps[,3],gps[,2],gps[,1])
#'
#' gpsarry<-array(c(gps1, gps2, gps3, gps4, gps5, gps6), dim=c(N,3,6))
#'
#'
#' #Create scalar balancing power parameters for each ordering of the GPS vector
#' fit6<- matrix(0,nrow=N,ncol=6,dimnames=list(c(1:N),c("ppar1","ppar2","ppar3",
#'               "ppar4","ppar5","ppar6")))
#'
#' \dontrun{
#' for(i in 1:6){
#'   fit6[,i]<-GPSCDF(pscores=gpsarry[,,i])$ppar
#' }
#'
#'   fit6
#'
#' #Perform analyses (similar to ordinal examples) using each K! ordering of the
#' #GPS vector. Select ordering which achieves optimal covariate balance
#' #(i.e. minimal standardized mean difference).
#' }
#'
#'@export

GPSCDF<-function(pscores=NULL, data=NULL, trt=NULL, stratify=FALSE, nstrat=5, optimal=FALSE, greedy=FALSE, ordinal=FALSE, multinomial=FALSE, caliper=NULL){

  if(is.null(stratify)){
    stratify<- FALSE
  }

  if(is.null(optimal)){
    optimal <- FALSE
  }

  if(is.null(greedy)){
    greedy <- FALSE
  }

  N<-dim(pscores)[1] # Number of subjects
  size<-dim(pscores)[2] # Number of treatments

  cpscores<-t(apply(pscores[,], 1,cumsum))
  Z<- seq(1, dim(pscores)[2], by=1)

  if(sum(pscores)/N != 1){
    message("Ensure pscores sum to 1")
  }

    Znorm<-sort(unique(Z))/max(unique(Z))
    ppar<-rep(0,N)

    for( i in 1:N){
      y<-cpscores[i,]
      mod<-stats::nls(y~I(Znorm^exp(power)), control = stats::nls.control(maxiter = 150, tol = 1e-05, minFactor = 1/1024,printEval = FALSE, warnOnly = TRUE), start = list(power = 0),trace = F)
      parm<-summary(mod)$coefficients[1]

      mod<-stats::nls(y~I(Znorm^exp(power)), control = stats::nls.control(maxiter = 150, tol = 1e-05, minFactor = 1/1024,printEval = FALSE, warnOnly = TRUE), start =list(power =parm),trace = F)
      parm<-summary(mod)$coefficients[1]

      ppar[i] <- parm

    }

    if (!is.null(data)){
      data$a<- ppar
      data2<- data
    }

    #Set up Stratification
    if(stratify==TRUE){
      strata<-dplyr::ntile(ppar, n=nstrat)

      if (!is.null(data)){
        data$strata<-strata
      }

    } else{
      strata<- NULL
      nstrat<- NULL}


    #Set up Optimal Matching
    if(optimal==TRUE){

      if (is.null(data)){
        stop('Specify a dataframe to attach matches')
      } else{

        if (is.null(trt)){
          stop('Specify a treatment variable to proceed with matching')
        } else{

          if (ordinal==FALSE & multinomial==FALSE){
            stop('Specify Ordinal or Multinomial treatments')
          } else{


          # Set up matching score
          epsilon=1e-5
          deltamat<-matrix(0,nrow=N,ncol=N)
          deltamat2<-matrix(0,nrow=N,ncol=N)

          # Loops to set up delta matrix
          if(ordinal==TRUE){
            for(i in 1:N){
              deltamat[i,]<-((ppar[i]-ppar)^2+epsilon)/((trt[i]-trt)^2)
            }
          }

          if(multinomial==TRUE){
            for(i in 1:N){
              for(k in 1:N){
                if(trt[i]==trt[k]){deltamat[i,k]<-999} else{
                  deltamat[i,k]<- abs((ppar[i]-ppar[k]))}
              }
            }
          }


          for(i in 1:N){
            deltamat2[i,]<-abs((ppar[i]-ppar))
          }

          #Get rid of Inf and put in 999999
          deltamat[!is.finite(deltamat)]<-99
          diag(deltamat)<-99

          # Derigs algorithm only works with integers so we can multiply all distances by 10,000 to get accuracy to 4 decimal places
          deltamatint<-deltamat*100000
          # Use distancematrix function to reform so we can do NBP matching
          suppressWarnings(distmat<-nbpMatching::distancematrix(deltamatint))

          # Set up matches
          invisible(utils::capture.output(matchset<-nbpMatching::nonbimatch(distmat)))

          #Remove row if N is odd
          matches1<-matchset$halves[ grep("ghost", matchset$halves$Group2.ID, invert = TRUE) , ]
          matches<- matches1[ grep("ghost", matches1$Group1.ID, invert = TRUE) , ]

          #Find distances of matches
          matchmat<-matrix(NA, nrow=round(dim(matches)[1]), ncol=3)
          for(i in 1:dim(matches)[1]){
            pair<-matches[i,c(2,4)]
            value<- deltamat2[pair[1,1],pair[1,2]]

            matchmat[i,1]<-pair[1,1]
            matchmat[i,2]<-pair[1,2]
            matchmat[i,3]<-value
          }

          npairs<-dim(matchmat)[1]

          #Calculate Average Total Distance of Matched Pairs
          optdistance<- sum(matchmat[,3])/npairs

          data$optmatch<-0
          # Attach matches to data
          for(i in 1:dim(data)[1]){
            data$optmatch[matches[i,2]]<-i
            data$optmatch[matches[i,4]]<-i
          }
          optmatch<-data$optmatch

          }
        }
      }
    } else{
      optmatch<- NULL
      optdistance<- NULL}


    #Set up Greedy Matching
    if(greedy==TRUE){

      if (is.null(data)){
        stop('Specify a dataframe to attach matches')
      } else{

        if (is.null(trt)){
          stop('Specify a treatment variable to proceed with matching')
        } else{

          if (ordinal==FALSE & multinomial==FALSE){
            stop('Specify Ordinal or Multinomial treatments')
          } else{

          # Set up matching score
          epsilon=1e-5
          deltamat<-matrix(0,nrow=N,ncol=N)
          deltamat2<-matrix(0,nrow=N,ncol=N)

          #Set Caliper
          if (is.null(caliper)){
            caliper<-.25*stats::sd(ppar)
          }

          # Loops to set up delta matrix
          if(ordinal==TRUE){
            for(i in 1:N){
              deltamat[i,]<-((ppar[i]-ppar)^2+epsilon)/((trt[i]-trt)^2)
            }
          }

          if(multinomial==TRUE){
            for(i in 1:N){
              for(k in 1:N){
                if(trt[i]==trt[k]){deltamat[i,k]<-999} else{
                  deltamat[i,k]<- abs((ppar[i]-ppar[k]))}
              }
            }
          }


          for(i in 1:N){
            deltamat2[i,]<-abs((ppar[i]-ppar))
          }

          #Get rid of Inf and put in 999
          deltamat[!is.finite(deltamat)]<-999
          diag(deltamat)<-999

          # Use Greedy matching to get matches
          # Set up matches

          # Set up holding for matches
          matchmat<-matrix(NA, nrow=round(dim(deltamat)[1]/2), ncol=3)

          #Replace matched pairs with maximum of delta matrix so it wont be used again
          repnum<-max(deltamat)

          i<-0
          while(min(deltamat) < caliper){
            i<-i+1
            inds = which(deltamat== min(deltamat), arr.ind=TRUE)
            value= deltamat2[inds[1,1],inds[1,2]]

            pair<-inds[1,1:2]
            matchmat[i,1:2]<-pair
            matchmat[i,3]<-value

            deltamat[pair,]<-repnum
            deltamat[,pair]<-repnum

          }

          matchmat<-matchmat[is.na(matchmat[,1])==F,]

          npairs<-dim(matchmat)[1]

          #Calculate Average Total Distance of Matched Pairs
          grdydistance<- sum(matchmat[,3])/npairs

          # Attach matches to data
          data2<- data2[matchmat,]

          data2$grdmatch<-0
          for(j in 1:(dim(data2)[1]/2)){
            data2$grdmatch[j]<-j
            data2$grdmatch[j+dim(data2)[1]/2]<-j
          }
          grdmatch<-data2$grdmatch
          grddata<-data2

          }
        }
      }
    } else{
      caliper<- NULL
      grdmatch<- NULL
      grddata<- NULL
      grdydistance<-NULL}



  returnlist<-list(ppar=ppar, data=data, nstrat=nstrat, strata=strata, optmatch=optmatch, optdistance=optdistance, caliper=caliper, grddata=grddata, grdmatch=grdmatch, grdydistance=grdydistance, NULL=NULL)
  returnlistfinal<- returnlist[-which(sapply(returnlist, is.null))]
  return(returnlistfinal)

} #END
