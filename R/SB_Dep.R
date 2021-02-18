#' Signed and HAC-Corrected Signed Tests of \insertCite{coudin2009finite;textual}{Rdpack}
#' 
#'
#' This function provides the nonparametric signed and HAC-corrected signed statstics proposed by \insertCite{coudin2009finite;textual}{Rdpack}. 
#' The tests are distribution-free, exact and valid in the presence of multiple regressors, and errors exhibiting serial nonlinear dependence and heterogeneous volatility.
#' The distribution of the test statistic is simulated under the null hypothesis with sufficient number of replications and the critical values are derived respectively. 
#' Unlike \insertCite{coudin2009finite;textual}{Rdpack}, this function is only limited to finite sample distributions. As such, the assumption of strict exogeneity must hold for the tests to be valid in small samples.
#' The R code for estimating the HAC covariance matrix has been retrieved from the paper by \insertCite{heberle2017fast;textual}{Rdpack}.  
#' 
#' @param y the vector of dependent variables.
#' @param x the \eqn{n\times 1} vector or \eqn{n\times k} matrix of regressors.
#' @param null the null hypothesis. Default vector is \eqn{(0,0)'} for a bivariate regression.
#' @param level the level of the test -i.e. \eqn{\alpha}. Default value is 0.05.
#' @param p the success probability of the binomial distribution for each trial. Default value is 0.5.
#' @param Dpndnt includes HAC-corrected signed statistic for serial (nonlinear) dependent data when set to TRUE.
#' @param bwidth the bandwidth for the Bartlett Kernell estimator.
#' @param B number of iterations for simulating the distribution of the test statistic under the null hypothesis.
#' @keywords Nonparametric, inference, dependent, exact, signs, HAC.
#' @export
#' @import matlib tsapp
#' @examples
#' SB_Dep(y,x,null=c(0,0),bwdith=2)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt



SB_Dep<-function(y,x,null=c(0,0),level=0.05,p=0.5,Dpndnt=FALSE,B=10000,bwidth){

	if(is.matrix(x)==TRUE){

		if(nrow(x)!=length(y)){

			"The number of rows of matrix X and the length of vector y must be equal!"

		}else{


			n<-length(y)

			X<-matrix(data=c(rep(1,time=n),x),nrow=n,ncol=ncol(x)+1)


			sgn_u<-1*((y-(X%*%null))>=0)

			P<-X%*%inv((t(X)%*%X)/n)%*%t(X)

			SF<-t(sgn_u)%*%P%*%sgn_u


			for(i in 1:B){

				sgn_u_Sim<-rbinom(n,1,p)
				SF_Sim<-t(sgn_u_Sim)%*%P%*%sgn_u_Sim

			}



			CritVal<-quantile(SF_Sim,1-level)


			print(paste("SF: The critical value at ",level," level is", CritVal))

			if(SF>CritVal){

				print(paste("SF: Rejected the null hypotehsis at ",level," level"))
				return(1)

			}else{

				print(paste("SF: Failed to Reject the null hypothesis at ",level," level"))
				return(0)

			}


			if(Dpndnt==TRUE){


				HAC.zei <-function(mcond,method="bartlett",bw){
					Nlen <- dim (mcond)[1]
					ww<-kweightsHAC(kernel=method,Nlen,bw)
					LL<-crossprod (mcond,mcond)/Nlen
					for(i in 1:bw){
						GG<-(crossprod(mcond[(i+1):Nlen,],mcond[1:(Nlen-i),]))/Nlen
						LL<-LL+(ww[i]*(GG+t(GG)))
					}
					return(LL)
				}



				A_m<-0

				for(j in 1:n){

					A_m<-A_m+(sgn_u[j]*X[j,])%*%(t(sgn_u[j]*X[j,]))


				}	

				A_m<-A_m/n

				


				cvHAC<-HAC.zei(A_m,bw=bwidth)

				Omeg<-cvHAC/n

				P_HAC<-X%*%inv(Omeg)%*%t(X)

				SHAC<-t(sgn_u)%*%P_HAC%*%sgn_u


				for(i in 1:B){

					sgn_u_Sim<-rbinom(n,1,p)	
					SHAC_Sim<-t(sgn_u_Sim)%*%P_HAC%*%sgn_u_Sim

				}



				CritVal_HAC<-quantile(SHAC_Sim,1-level)


				print(paste("SHAC_Sim: The critical value at ",level," level is", CritVal_HAC))

				if(SHAC>CritVal_HAC){

					print(paste("SHAC: Rejected the null hypotehsis at ",level," level"))
					return(1)

				}else{

					print(paste("SHAC: Failed to Reject the null hypothesis at ",level," level"))
					return(0)

				}



			}



		}


	}else{


		if(length(x)!=length(y)){

			"The lengths of the vectors x and y must be equal!"

		}else{


			n<-length(y)

			X<-matrix(data=c(rep(1,time=n),x),nrow=n,ncol=2)


			sgn_u<-1*((y-(X%*%null))>=0)

			P<-X%*%inv((t(X)%*%X)/n)%*%t(X)

			SF<-t(sgn_u)%*%P%*%sgn_u


			for(i in 1:B){

				sgn_u_Sim<-rbinom(n,1,p)
				SF_Sim<-t(sgn_u_Sim)%*%P%*%sgn_u_Sim

			}



			CritVal<-quantile(SF_Sim,1-level)
			

			print(paste("SF: The critical value at ",level," level is", CritVal))

			if(SF>CritVal){

				print(paste("SF: Rejected the null hypotehsis at ",level," level"))
				return(1)

			}else{

				print(paste("SF: Failed to Reject the null hypothesis at ",level," level"))
				return(0)

			}


			if(Dpndnt==TRUE){
				

				HAC.zei <-function(mcond,method="bartlett",bw){
					Nlen <- dim (mcond)[1]
					ww<-kweightsHAC(kernel=method,Nlen,bw)
					LL<-crossprod (mcond,mcond)/Nlen
					for(i in 1:bw){
						GG<-(crossprod(mcond[(i+1):Nlen,],mcond[1:(Nlen-i),]))/Nlen
						LL<-LL+(ww[i]*(GG+t(GG)))
					}
					return(LL)
				}



				A_m<-0

				for(j in 1:n){

					A_m<-A_m+(sgn_u[j]*X[j,])%*%(t(sgn_u[j]*X[j,]))


				}	

				A_m<-A_m/n



				cvHAC<-HAC.zei(A_m,bw=1)

				Omeg<-cvHAC/n

				P_HAC<-X%*%inv(Omeg)%*%t(X)

				SHAC<-t(sgn_u)%*%P_HAC%*%sgn_u


				for(i in 1:B){

					sgn_u_Sim<-rbinom(n,1,p)	
					SHAC_Sim<-t(sgn_u_Sim)%*%P_HAC%*%sgn_u_Sim

				}



				CritVal_HAC<-quantile(SHAC_Sim,1-level)


				print(paste("SHAC_Sim: The critical value at ",level," level is", CritVal_HAC))

				if(SHAC>CritVal_HAC){

					print(paste("SHAC: Rejected the null hypotehsis at ",level," level"))
					return(1)

				}else{

					print(paste("SHAC: Failed to Reject the null hypothesis at ",level," level"))
					return(0)

				}



			}



		}


	}
}