#' Signed and Signed Rank Tests of \insertCite{dufour1995exact;textual}{Rdpack}
#'
#' This function provides the test statistic and the critical value for the nonparametric sign-based and signed-rank tests proposed by \insertCite{dufour1995exact;textual}{Rdpack}. 
#' These tests are valid in the presence of a single regressor and no nuisance parameters (e.g. intercept).
#' 
#' @param y the vector of dependent variables
#' @param x the \eqn{n\times 1} vector or \eqn{n\times k} matrix of regressors.
#' @param level the level of the test -i.e. \eqn{\alpha}. Default value is 0.05.
#' @param p is the success probability of the binomial distribution for each trial. Default value is 0.5.
#' @param W includes the Wilcoxon signed rank test variate when set to TRUE.
#' @param SR includes signed rank test variate proposed by \insertCite{dufour1995exact;textual}{Rdpack} when set to TRUE.
#' @keywords Nonparametric, inference, exact, signs
#' @export
#' @examples
#' CD_95(y,x)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt

CD_95<-function(y,x,level=0.05,p=0.5,W=FALSE,SR=FALSE){

	if(length(x)!=length(y)){

		"Lengths of x and y must be equal!"

	}else{



		n<-length(y)

		sgn_z<-1*((y*x)>=0)

		CDStat<-sum(sgn_z)

		print(CDStat)

		print(paste("Signed: CD_95 Sign-based test statistic: ",CDStat))

		invisible(CDStat)



		CritVal_U<-qbinom(1-(level/2),n,p)
		CritVal_L<-n-CritVal_U

		print(paste("Signed: Two-sided critical values at ",level," level are ", CritVal_L," and ",CritVal_U))

		invisible(CritVal_U)
		invisible(CritVal_L)


		if(CDStat>CritVal_U || CDStat<CritVal_L){

			print(paste("Signed: Rejected the null hypotehsis at ",level," level"))
			Sign<-1

		}else{

			print(paste("Signed: Failed to Reject the null hypothesis at ",level," level"))
			Sign<-0

		}



		if(W==TRUE){

			Z_1<-abs(y*x)

			W_1<-order(Z_1)

			CDRStat<-sum(sgn_z*W_1)

			print(CDRStat)

			print(paste("Wilcoxon-Signed rank: CD_95 test statistic: ",CDRStat))

			invisible(CDRStat)


			CritValR_U<-qsignrank(1-(level/2),n)		
			CritValR_L<-(n*(n+1)/2)-CritValR_U

			print(paste("Wilcoxon-Signed rank: TWo-sided critical values at the ",level," level are ", CritValR_L," and ",CritVal_U))

			invisible(CritValR_U)
			invisible(CritValR_L)


			if(CDRStat>CritValR_U || CDRStat<CritValR_L){

				print(paste("Wilcoxon-Signed rank: Rejected the null hypotehsis at ",level," level"))
				WilcoxonSign<-1

			}else{

				print(paste("Wilcoxon-Signed rank: Failed to Reject the null hypothesis at ",level," level"))
				WilcoxonSign<-0

			}

		}

		if(SR==TRUE){

			Z_2<-abs(y)

			W_2<-order(Z_2)

			CDWStat<-sum(sgn_z*W_2)

			print(CDWStat)

			print(paste("Signed rank: CD_95 test statistic: ",CDWStat))

			invisible(CDWStat)

			CritValW_U<-qsignrank(1-(level/2),n)
			CritValW_L<-(n*(n+1)/2)-CritValW_U
			

			print(paste("Signed rank: TWo-sided critical values at the ",level," level are ", CritValW_L," and ",CritValW_U))

			invisible(CritValW_L)
			invisible(CritValW_U)


			if(CDWStat>CritValW_U || CDWStat<CritValW_L){

				print(paste("Signed rank: Rejected the null hypotehsis at ",level," level"))
				SignRank<-1

			}else{

				print(paste("Signed rank: Failed to Reject the null hypothesis at ",level," level"))
				SignRank<-0

			}

		}


		if(W==TRUE && SR==TRUE){

			R_95<-c(Sign,WilcoxonSign,SignRank)

		}else if(W==TRUE && SR==FALSE){

			R_95<-c(Sign,WilcoxonSign)

		}else if(W==FALSE && SR==TRUE){

			R_95<-c(Sign,SignRank)

		}

		return(R_95)

	}

}
