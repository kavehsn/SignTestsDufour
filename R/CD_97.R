#' Signed and Signed Rank Tests of \insertCite{campbell1997exact;textual}{Rdpack}
#'
#' This function extends the nonparametric signed and signed rank tests of \insertCite{campbell1997exact;textual}{Rdpack} by allowing the presence of a nuisance parameter. 
#' The nuisance parameter problem is dealt with by employing bound-type procedures. These tests are valid in the presence of a single regressor and a nuisance parameter (e.g. intercept).
#' 
#' @param y the vector of dependent variables
#' @param x the \eqn{n\times 1} vector or \eqn{n\times k} matrix of regressors.
#' @param level the level of the test -i.e. \eqn{\alpha}. Default value is 0.05.
#' @param alpha_1 or \eqn{\alpha_1} for empolying the bound-type procedure. Default value is \eqn{0.14\times\alpha}. Note: \eqn{\alpha_1<\alpha}.
#' @param p the success probability of the binomial distribution for each trial. Default value is 0.5.
#' @param SR includes signed rank test variate proposed by Campbell and Dufour (1995) when set to TRUE.
#' @keywords Nonparametric, inference, exact, signs, Bonferroni, nuisance-parameter
#' @export
#' @examples
#' CD_97(y,x)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt

CD_97<-function(y,x,level=0.05,alpha_1=0.14*level,p=0.5,SRTest=FALSE){

	if(length(x)!=length(y)){

		"Lengths of x and y must be equal!"

	}else{

		n<-length(y)

		alpha_2<-level-alpha_1
		alpha_3<-level+alpha_1

		k<-qbinom(alpha_1/2,n,p)

		LowerQuantile<-k+1
		UpperQuantile<-n-k

		OrderStats<-y[order(y)]

		CI<-OrderStats[LowerQuantile:UpperQuantile]


		SR<-rep(0,times=length(CI))
		SL<-rep(0,times=length(CI))
		STT<-rep(0,times=length(CI))

		for(i in 1:length(SR)){

			sgn_z<-1*((y-CI[i])*x>=0)

			SR[i]<-sum(sgn_z)
			SL[i]<-n-SR[i]
			STT[i]<-max(SR[i],SL[i])

		}

		Q_L<-min(STT)
		Q_U<-max(STT)

		CritVal_2<-qbinom(1-(alpha_2/2),n,p)
		CritVal_3<-n-qbinom((alpha_3/2),n,p)


		if(Q_L>CritVal_2){

			print(paste("Signed: Rejected the null hypotehsis at ",level," level"))
			Sign<-1

		}else if(Q_U<=CritVal_3){

			print(paste("Signed: Failed to reject the null hypotehsis at ",level," level"))
			Sign<-0

		}else{

			print(paste("Signed: The test is inconclusive at ",level," level"))
			Sign<-0

		}

		if(SRTest!=FALSE){

			SRR<-rep(0,times=length(CI))
			SRL<-rep(0,times=length(CI))
			SRTT<-rep(0,times=length(CI))

			for(i in 1:length(SRR)){

				sgn_z<-1*((y-CI[i])*x>=0)

				Z_2<-abs(y-CI[i])

				W_2<-order(Z_2)


				SRR[i]<-sum(sgn_z*W_2)
				SRL[i]<-(n*(n+1)/2)-SRR[i]
				SRTT[i]<-max(SRR[i],SRL[i])

			}


			CritVal_2_R<-qsignrank(1-(alpha_2/2),n,p)
			CritVal_3_R<-(n*(n+1)/2)-qbinom((alpha_3/2),n,p)


			QR_L<-min(SRTT)
			QR_U<-max(SRTT)


			if(QR_L>CritVal_2_R){

				print(paste("SignedRank: Rejected the null hypotehsis at ",level," level"))
				SignRank<-1

			}else if(QR_U<=CritVal_3_R){

				print(paste("SignedRank: Failed to reject the null hypotehsis at ",level," level"))
				SignRank<-0

			}else{

				print(paste("SignedRank: The test is inconclusive at ",level," level"))
				SignRank<-0

			}


		}



		if(SRTest==TRUE){

			R<-c(Sign,SignRank)

		}else{

			R<-c(Sign)

		}


	}

}
