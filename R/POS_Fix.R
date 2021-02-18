#' Point-Optimal Sign-Based Tests of \insertCite{dufour2010exact;textual}{Rdpack}
#'
#' This function provides the test statistic and the critical values for the nonparametric point-optimal sign-based tests proposed by \insertCite{dufour2010exact;textual}{Rdpack}. 
#' The proposed tests are exact, distribution-free and valid in the presence of nonstandard distributions and heteroskedasticity of unknown form. Moreover, they have the highest power amongst commonly encountered tests that are intended to be
#' robust against heteroskedasticity.
#' 
#' @param y the vector of dependent variables
#' @param x the vector of regressors
#' @param null the null hypothesis
#' @param level is the level of the test. Default value is 0.05.
#' @param p is the success probability of the binomial distribution for each trial. Default value is 0.5.
#' @param B is the number iterations for simulating the 
#' @keywords Nonparametric point-optimal Sign-tests Dufour and Taamouti (2010)
#' @import MASS
#' @export
#' @examples
#' POS_Fix(y,x,null=c(0,0))
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt

POS_Fix<-function(y,x,null=c(0,0),level=0.05,p=0.5,B=10000,...){

	if(is.matrix(x)==TRUE){

		if(nrow(x)!=length(y)){

			"The number of rows of matrix X and the length of vector y must be equal!"

		}else{

			n<-length(y)

			z<-matrix(data=c(y,x),nrow=n,ncol=ncol(x)+1)

			smp_size<-floor(0.1*nrow(z))

			z_Alt<-head(z,smp_size)
			z_Test<-tail(z,n-smp_size)

			y_Alt<-z_Alt[,1]
			x_Alt<-z_Alt[,2:ncol(z_Alt)]

			y_Test<-z_Test[,1]
			x_Test<-z_Test[,2:ncol(z_Test)]


			X_Test<-matrix(data=c(rep(1,times=nrow(x_Test)),x_Test),nrow=nrow(x_Test),ncol=ncol(x)+1)

			sgn_y<-1*((y_Test-(X_Test%*%null))>=0)

			betahat<-rlm(y_Alt~x_Alt,...)


			a1<-(1/pnorm(X_Test%*%(betahat$coefficients-null)))-1

			POSStat<-sum(log(1/a1)*sgn_y)


			print(POSStat)
			print(paste("POS_Fix test statistic: ",POSStat))
			invisible(POSStat)

			POSStat_Sim <- rep(0, times = B)

			for(i in 1:B){

				sgn_y_sim<-rbinom(length(y_Test),1,p)
				POSStat_Sim[i]<-sum(log(1/a1)*sgn_y_sim)

			}


			CritVal<-quantile(POSStat_Sim,1-level)

			print(paste("The critical value at ",level," level is", CritVal))

			if(POSStat>CritVal){

				print(paste("Rejected the null hypotehsis at ",level," level"))

			}else{

				print(paste("Failed to Reject the null hypothesis at ",level," level"))

			}




		}
	}else{


		if(length(x)!=length(y)){

			"The lengths of the vectors x and y must be equal!"

		}else{

			n<-length(y)

			z<-matrix(data=c(y,x),nrow=n,ncol=2)

			smp_size<-floor(0.1*nrow(z))

			z_Alt<-head(z,smp_size)
			z_Test<-tail(z,n-smp_size)

			y_Alt<-z_Alt[,1]
			x_Alt<-z_Alt[,2:ncol(z_Alt)]

			y_Test<-z_Test[,1]
			x_Test<-z_Test[,2:ncol(z_Test)]

			X_Test<-matrix(data=c(rep(1,times=length(x_Test)),x_Test),nrow=length(x_Test),ncol=2)

			sgn_y<-1*((y_Test-(X_Test%*%null))>=0)

			betahat<-rlm(y_Alt~x_Alt,...)


			a1<-(1/pnorm(X_Test%*%(betahat$coefficients-null)))-1

			POSStat<-sum(log(1/a1)*sgn_y)


			print(POSStat)
			print(paste("POS_Fix test statistic: ",POSStat))
			invisible(POSStat)

			POSStat_Sim <- rep(0, times = B)

			for(i in 1:B){

				sgn_y_sim<-rbinom(length(y_Test),1,p)
				POSStat_Sim[i]<-sum(log(1/a1)*sgn_y_sim)

			}


			CritVal<-quantile(POSStat_Sim,1-level)

			print(paste("The critical value at ",level," level is", CritVal))

			if(POSStat>CritVal){

				print(paste("Rejected the null hypotehsis at ",level," level"))

			}else{

				print(paste("Failed to Reject the null hypothesis at ",level," level"))

			}



		}


	}
}
