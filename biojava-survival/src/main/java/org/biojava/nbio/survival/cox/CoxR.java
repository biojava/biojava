/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.survival.cox;

import org.biojava.nbio.survival.cox.matrix.Matrix;
import org.biojava.nbio.survival.cox.stats.ChiSq;
import org.biojava.nbio.survival.cox.stats.Cholesky2;
import org.biojava.nbio.survival.data.WorkSheet;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;

/**
 *   This is a port of the R survival code used for doing Cox Regression. The algorithm was a fairly easy port from C code to Java where the challenge was
 *   making the code a little more object friendly. In the R code everything is passed around as an array and a large portion of the code is spent extracting
 *   data from the array for use in different calculations. By organizing the data in a class for each data point was able to simplify much of the code.
 *   Not all variants of different methods that you can select for doing various statistical calculations are implemented. Wouldn't be difficult to go back in
 *   add them in if they are important.
 *
 *<p>In R you can pass in different paramaters to override defaults which requires parsing of the paramaters. In the Java code tried to be a little more exact
 *   in the code related to paramaters where using strata, weighting, robust and cluster are advance options. Additionaly code is implemented from Bob Gray
 *   to do variance correction when using weighted paramaters in a data set.
 *   /Users/Scooter/NetBeansProjects/biojava3-survival/docs/wtexamples.docx
 *
 *<p>The CoxHelper class is meant to hide some of the implementation details.
 *
 *<p>Issues
 *<ul>
 *<li>sign in CoxMart?
 *<li>double toler_chol = 1.818989e-12; Different value for some reason
 *<li>In robust linear_predictor set to 0 which implies score = 1 but previous score value doesn't get reset
 *</ul>
 *
 *  Cox regression fit, replacement for coxfit2 in order
 *    to be more frugal about memory: specificly that we
 *    don't make copies of the input data.
 *
 *
 *
 *
 *
 * <p>
 * the input parameters are
 *
 * <pre>
 *      maxiter      :number of iterations
 *      time(n)      :time of status or censoring for person i
 *      status(n)    :status for the ith person    1=dead , 0=censored
 *      covar(nv,n)  :covariates for person i.
 *                       Note that S sends this in column major order.
 *      strata(n)    :marks the strata.  Will be 1 if this person is the
 *                      last one in a strata.  If there are no strata, the
 *                      vector can be identically zero, since the nth person's
 *                      value is always assumed to be = to 1.
 *      offset(n)    :offset for the linear predictor
 *      weights(n)   :case weights
 *      init         :initial estimate for the coefficients
 *      eps          :tolerance for convergence.  Iteration continues until
 *                      the percent change in loglikelihood is <= eps.
 *      chol_tol     : tolerance for the Cholesky decompostion
 *      method       : 0=Breslow, 1=Efron
 *      doscale      : 0=don't scale the X matrix, 1=scale the X matrix
 * </pre>
 * returned parameters
 * <pre>
 *      means(nv)    : vector of column means of X
 *      beta(nv)     :the vector of answers (at start contains initial est)
 *      u(nv)        :score vector
 *      imat(nv,nv)  :the variance matrix at beta=final
 *                     (returned as a vector)
 *      loglik(2)    :loglik at beta=initial values, at beta=final
 *      sctest       :the score test at beta=initial
 *      flag         :success flag  1000  did not converge
 *                                  1 to nvar: rank of the solution
 *      iterations         :actual number of iterations used
 * </pre>
 * work arrays
 * <pre>
 *      mark(n)
 *      wtave(n)
 *      a(nvar), a2(nvar)
 *      cmat(nvar,nvar)       ragged array
 *      cmat2(nvar,nvar)
 *      newbeta(nvar)         always contains the "next iteration"
 * </pre>
 * calls functions:  cholesky2, chsolve2, chinv2
 * <p>
 * the data must be sorted by ascending time within strata
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxR {

	/**
	 *
	 * @param variables
	 * @param DataT
	 * @param useStrata
	 * @param useWeighted
	 * @param robust
	 * @param cluster
	 * @return
	 * @throws Exception
	 */
	public CoxInfo process(ArrayList<String> variables, ArrayList<SurvivalInfo> DataT, boolean useStrata, boolean useWeighted, boolean robust, boolean cluster) throws Exception {
		//from coxph.control.S
		int maxiter2 = 20;
		double eps2 = 1e-9;
		double toler2 = Math.pow(eps2, .75);
		int doscale2 = 1;
		//int method2 = 0;
		double[] beta = new double[variables.size()];
		return process(variables, DataT, maxiter2, CoxMethod.Efron, eps2, toler2, beta, doscale2, useStrata, useWeighted, robust, cluster);

	}

	/**
	 *
	 * @param variables
	 * @param data
	 * @param maxiter
	 * @param method
	 * @param eps
	 * @param toler
	 * @param beta
	 * @param doscale
	 * @param useStrata
	 * @param useWeighted
	 * @param robust
	 * @param cluster
	 * @return
	 * @throws Exception
	 */
	public CoxInfo process(ArrayList<String> variables, ArrayList<SurvivalInfo> data, int maxiter, CoxMethod method, double eps, double toler, double[] beta, int doscale, boolean useStrata, boolean useWeighted, boolean robust, boolean cluster) throws Exception {
		//make sure data is converted to numbers if labels are used
		SurvivalInfoHelper.categorizeData(data);
		//create variables if testing for interaction
		for (String variable : variables) {
			if (variable.indexOf(":") != -1) {
				String[] d = variable.split(":");
				SurvivalInfoHelper.addInteraction(d[0], d[1], data);
			}
		}

		Collections.sort(data);
		// Collections.reverse(data);
		CoxInfo coxInfo = new CoxInfo();
		coxInfo.setSurvivalInfoList(data);



		int i, j, k, person;
		boolean gotofinish = false;
		double[][] cmat, imat;  /*ragged arrays covar[][], */
		double wtave;
		double[] a, newbeta;
		double[] a2;
		double[][] cmat2;
		double[] scale;
		double denom = 0, zbeta, risk;
		double temp, temp2;
		int ndead;  /* actually, the sum of their weights */
		double newlk = 0;
		double dtime, d2;
		double deadwt;  /*sum of case weights for the deaths*/
		double efronwt; /* sum of weighted risk scores for the deaths*/
		int halving;    /*are we doing step halving at the moment? */
		@SuppressWarnings("unused")
		int nrisk = 0;   /* number of subjects in the current risk set */

		/* copies of scalar input arguments */
		int nused, nvar;




		/* vector inputs */
		//  double *time, *weights, *offset;
		//  int *status, *strata;

		/* returned objects */
		// double imat2[][];
		double[] u, loglik, means;


		double sctest;
		int flag = 0;
		int iter = 0;
		//SEXP rlist, rlistnames;
		//  int nprotect;  /* number of protect calls I have issued */

		/* get local copies of some input args */
		nused = data.size(); // LENGTH(offset2);
		nvar = variables.size(); // ncols(covar2);


		//       imat2 = new double[nvar][nvar];
//        nprotect++;
		imat = new double[nvar][nvar]; //dmatrix(REAL(imat2),  nvar, nvar);
		a = new double[nvar]; //(double *) R_alloc(2*nvar*nvar + 4*nvar, sizeof(double));
		newbeta = new double[nvar]; //a + nvar;
		a2 = new double[nvar]; //newbeta + nvar;
		scale = new double[nvar]; //a2 + nvar;
		cmat = new double[nvar][nvar]; //dmatrix(scale + nvar,   nvar, nvar);
		cmat2 = new double[nvar][nvar]; //dmatrix(scale + nvar +nvar*nvar, nvar, nvar);

		/*
		 ** create output variables
		 */
//    PROTECT(beta2 = duplicate(ibeta));
//    beta = REAL(beta2);
		//  beta = new double[nvar];
		// beta = beta2;
		//  PROTECT(means2 = allocVector(REALSXP, nvar));
		//  means = REAL(means2);
		means = new double[nvar];
		double[] sd = new double[nvar];
		//double[] se = new double[nvar];

		//   means = means2;
		//   PROTECT(u2 = allocVector(REALSXP, nvar));
		//   u = REAL(u2);
		u = new double[nvar];
		//   u = u2;
//    PROTECT(loglik2 = allocVector(REALSXP, 2));
//    loglik = REAL(loglik2);
		loglik = new double[2];
		//   loglik = loglik2;
//    PROTECT(sctest2 = allocVector(REALSXP, 1));
//    sctest = REAL(sctest2);
//        sctest = new double[1];
		//   sctest = sctest2;
//    PROTECT(flag2 = allocVector(INTSXP, 1));
//    flag = INTEGER(flag2);
//        flag = new int[1];
		//     flag = flag2;
//    PROTECT(iter2 = allocVector(INTSXP, 1));
//    iterations = INTEGER(iter2);
//        iterations = new int[1];
//        iterations = iter2;
		//       nprotect += 7;

		/*
		 ** Subtract the mean from each covar, as this makes the regression
		 **  much more stable.
		 */
		double[] time = new double[nused];
		int[] status = new int[nused];
		double[] offset = new double[nused];
		double[] weights = new double[nused];
		int[] strata = new int[nused];

		double[][] covar = new double[nvar][nused];
		ArrayList<String> clusterList = null;

		if(cluster){
			clusterList = new ArrayList<String>();
		}
		//copy data over to local arrays to minimuze changing code
		for (person = 0; person < nused; person++) {
			SurvivalInfo si = data.get(person);
			time[person] = si.getTime();
			status[person] = si.getStatus();
			offset[person] = si.getOffset();
			if(cluster){
				if(si.getClusterValue() == null && si.getClusterValue().length() == 0){
					throw new Exception("Cluster value is not valid for " + si.toString());
				}
				clusterList.add(si.getClusterValue());
			}
			if (useWeighted) {
				weights[person] = si.getWeight();
			} else {
				weights[person] = 1.0;
			}
			if (useStrata) {
				strata[person] = si.getStrata();
			} else {
				strata[person] = 0;
			}
			for (i = 0; i < variables.size(); i++) {
				String variable = variables.get(i);
				covar[i][person] = si.getVariable(variable);
			}
		}

		double tempsd = 0;
		i = 0;
		for (i = 0; i < nvar; i++) {

			temp = 0;
			tempsd = 0;
			//calculate the mean sd

			for (person = 0; person < nused; person++) {

				temp += covar[i][person]; // * weights[person];
				tempsd += (covar[i][person]) * (covar[i][person]); //*weights[person] * weights[person]
			}
			temp /= nused;
			//   temp /= weightCount;
			means[i] = temp;
			tempsd /= nused;
			//  tempsd /= weightCount;
			tempsd = Math.sqrt(tempsd - temp * temp);
			sd[i] = tempsd; //standard deviation
			//subtract the mean
			for (person = 0; person < nused; person++) {
				covar[i][person] -= temp;
			}
			if (doscale == 1) {  /* and also scale it */
				temp = 0;
				for (person = 0; person < nused; person++) {
					temp += Math.abs(covar[i][person]); //fabs
				}
				if (temp > 0) {
					temp = nused / temp;   /* scaling */
				} else {
					temp = 1.0; /* rare case of a constant covariate */
				}
				scale[i] = temp;
				for (person = 0; person < nused; person++) {
					covar[i][person] *= temp;
				}
			}
		}
		if (doscale == 1) {
			for (i = 0; i < nvar; i++) {
				beta[i] /= scale[i]; /*rescale initial betas */
			}
		} else {
			for (i = 0; i < nvar; i++) {
				scale[i] = 1.0;
			}
		}

		/*
		 ** do the initial iteration step
		 */
		strata[nused - 1] = 1;
		loglik[1] = 0;
		for (i = 0; i < nvar; i++) {
			u[i] = 0;  //u = s1
			a2[i] = 0; //a2 = a
			for (j = 0; j < nvar; j++) {
				imat[i][j] = 0;  //s2
				cmat2[i][j] = 0; //a
			}
		}

		for (person = nused - 1; person >= 0;) {
			if (strata[person] == 1) {
				nrisk = 0;
				denom = 0;
				for (i = 0; i < nvar; i++) {
					a[i] = 0;
					for (j = 0; j < nvar; j++) {
						cmat[i][j] = 0;
					}
				}
			}

			dtime = time[person];
			ndead = 0; /*number of deaths at this time point */
			deadwt = 0;  /* sum of weights for the deaths */
			efronwt = 0;  /* sum of weighted risks for the deaths */
			while (person >= 0 && time[person] == dtime) {
				/* walk through the this set of tied times */
				nrisk++;
				zbeta = offset[person];    /* form the term beta*z (vector mult) */
				for (i = 0; i < nvar; i++) {
					zbeta += beta[i] * covar[i][person]; //x
				}
				zbeta = coxsafe(zbeta);

				risk = Math.exp(zbeta) * weights[person]; //risk = v
				denom += risk;

				/* a is the vector of weighted sums of x, cmat sums of squares */
				for (i = 0; i < nvar; i++) {
					a[i] += risk * covar[i][person]; //a = s1
					for (j = 0; j <= i; j++) {
						cmat[i][j] += risk * covar[i][person] * covar[j][person]; //cmat = s2;
					}
				}

				if (status[person] == 1) {
					ndead++;
					deadwt += weights[person];
					efronwt += risk;
					loglik[1] += weights[person] * zbeta;

					for (i = 0; i < nvar; i++) {
						u[i] += weights[person] * covar[i][person];
					}
					if (method == CoxMethod.Efron) { /* Efron */
						for (i = 0; i < nvar; i++) {
							a2[i] += risk * covar[i][person];
							for (j = 0; j <= i; j++) {
								cmat2[i][j] += risk * covar[i][person] * covar[j][person];
							}
						}
					}
				}

				person--;
				if (person >= 0 && strata[person] == 1) { //added catch of person = 0 and person-- = -1
					break;  /*ties don't cross strata */
				}
			}


			if (ndead > 0) {  /* we need to add to the main terms */
				if (method == CoxMethod.Breslow) { /* Breslow */
					loglik[1] -= deadwt * Math.log(denom);

					for (i = 0; i < nvar; i++) {
						temp2 = a[i] / denom;  /* mean */
						u[i] -= deadwt * temp2;
						for (j = 0; j <= i; j++) {
							imat[j][i] += deadwt * (cmat[i][j] - temp2 * a[j]) / denom;
						}
					}
				} else { /* Efron */
					/*
					 ** If there are 3 deaths we have 3 terms: in the first the
					 **  three deaths are all in, in the second they are 2/3
					 **  in the sums, and in the last 1/3 in the sum.  Let k go
					 **  from 0 to (ndead -1), then we will sequentially use
					 **     denom - (k/ndead)*efronwt as the denominator
					 **     a - (k/ndead)*a2 as the "a" term
					 **     cmat - (k/ndead)*cmat2 as the "cmat" term
					 **  and reprise the equations just above.
					 */
					for (k = 0; k < ndead; k++) {
						temp = (double) k / ndead;
						wtave = deadwt / ndead;
						d2 = denom - temp * efronwt;
						loglik[1] -= wtave * Math.log(d2);
						for (i = 0; i < nvar; i++) {
							temp2 = (a[i] - temp * a2[i]) / d2;
							u[i] -= wtave * temp2;
							for (j = 0; j <= i; j++) {
								imat[j][i] += (wtave / d2)
										* ((cmat[i][j] - temp * cmat2[i][j])
										- temp2 * (a[j] - temp * a2[j]));
							}
						}
					}

					for (i = 0; i < nvar; i++) {
						a2[i] = 0;
						for (j = 0; j < nvar; j++) {
							cmat2[i][j] = 0;
						}
					}
				}
			}
		}   /* end  of accumulation loop */
		loglik[0] = loglik[1]; /* save the loglik for iterations 0 */

		/* am I done?
		 **   update the betas and test for convergence
		 */
		for (i = 0; i < nvar; i++) /*use 'a' as a temp to save u0, for the score test*/ {
			a[i] = u[i];
		}

		flag = Cholesky2.process(imat, nvar, toler);
		chsolve2(imat, nvar, a);        /* a replaced by  a *inverse(i) */

		temp = 0;
		for (i = 0; i < nvar; i++) {
			temp += u[i] * a[i];
		}
		sctest = temp;  /* score test */

		/*
		 **  Never, never complain about convergence on the first step.  That way,
		 **  if someone HAS to they can force one iterations at a time.
		 */
		for (i = 0; i < nvar; i++) {
			newbeta[i] = beta[i] + a[i];
		}
		if (maxiter == 0) {
			chinv2(imat, nvar);
			for (i = 0; i < nvar; i++) {
				beta[i] *= scale[i];  /*return to original scale */
				u[i] /= scale[i];
				imat[i][i] *= scale[i] * scale[i];
				for (j = 0; j < i; j++) {
					imat[j][i] *= scale[i] * scale[j];
					imat[i][j] = imat[j][i];
				}
			}
			// goto finish;
			gotofinish = true;

		}

		/*
		 ** here is the main loop
		 */
		if (!gotofinish) {
			halving = 0;             /* =1 when in the midst of "step halving" */
			for (iter = 1; iter <= maxiter; iter++) {
				newlk = 0;
				for (i = 0; i < nvar; i++) {
					u[i] = 0;
					for (j = 0; j < nvar; j++) {
						imat[i][j] = 0;
					}
				}

				/*
				 ** The data is sorted from smallest time to largest
				 ** Start at the largest time, accumulating the risk set 1 by 1
				 */
				for (person = nused - 1; person >= 0;) {
					if (strata[person] == 1) { /* rezero temps for each strata */
						denom = 0;
						nrisk = 0;
						for (i = 0; i < nvar; i++) {
							a[i] = 0;
							for (j = 0; j < nvar; j++) {
								cmat[i][j] = 0;
							}
						}
					}

					dtime = time[person];
					deadwt = 0;
					ndead = 0;
					efronwt = 0;
					while (person >= 0 && time[person] == dtime) {
						nrisk++;
						zbeta = offset[person];
						for (i = 0; i < nvar; i++) {
							zbeta += newbeta[i] * covar[i][person];
						}
						zbeta = coxsafe(zbeta);


						risk = Math.exp(zbeta) * weights[person];
						denom += risk;

						for (i = 0; i < nvar; i++) {
							a[i] += risk * covar[i][person];
							for (j = 0; j <= i; j++) {
								cmat[i][j] += risk * covar[i][person] * covar[j][person];
							}
						}

						if (status[person] == 1) {
							ndead++;
							deadwt += weights[person];
							newlk += weights[person] * zbeta;
							for (i = 0; i < nvar; i++) {
								u[i] += weights[person] * covar[i][person];
							}
							if (method == CoxMethod.Efron) { /* Efron */
								efronwt += risk;
								for (i = 0; i < nvar; i++) {
									a2[i] += risk * covar[i][person];
									for (j = 0; j <= i; j++) {
										cmat2[i][j] += risk * covar[i][person] * covar[j][person];
									}
								}
							}
						}

						person--;
						if (person >= 0 && strata[person] == 1) { //added catch of person = 0 and person-- = -1
							break;  /*ties don't cross strata */
						}
					}

					if (ndead > 0) {  /* add up terms*/
						if (method == CoxMethod.Breslow) { /* Breslow */
							newlk -= deadwt * Math.log(denom);
							for (i = 0; i < nvar; i++) {
								temp2 = a[i] / denom;  /* mean */
								u[i] -= deadwt * temp2;
								for (j = 0; j <= i; j++) {
									imat[j][i] += (deadwt / denom)
											* (cmat[i][j] - temp2 * a[j]);
								}
							}
						} else { /* Efron */
							for (k = 0; k < ndead; k++) {
								temp = (double) k / ndead;
								wtave = deadwt / ndead;
								d2 = denom - temp * efronwt;
								newlk -= wtave * Math.log(d2);
								for (i = 0; i < nvar; i++) {
									temp2 = (a[i] - temp * a2[i]) / d2;
									u[i] -= wtave * temp2;
									for (j = 0; j <= i; j++) {
										imat[j][i] += (wtave / d2)
												* ((cmat[i][j] - temp * cmat2[i][j])
												- temp2 * (a[j] - temp * a2[j]));
									}
								}
							}

							for (i = 0; i < nvar; i++) { /*in anticipation */
								a2[i] = 0;
								for (j = 0; j < nvar; j++) {
									cmat2[i][j] = 0;
								}
							}
						}
					}
				}   /* end  of accumulation loop  */

				/* am I done?
				 **   update the betas and test for convergence
				 */
				flag = Cholesky2.process(imat, nvar, toler);

				if (Math.abs(1 - (loglik[1] / newlk)) <= eps && halving == 0) { /* all done */
					loglik[1] = newlk;
					chinv2(imat, nvar);     /* invert the information matrix */
					for (i = 0; i < nvar; i++) {
						beta[i] = newbeta[i] * scale[i];
						u[i] /= scale[i];
						imat[i][i] *= scale[i] * scale[i];
						for (j = 0; j < i; j++) {
							imat[j][i] *= scale[i] * scale[j];
							imat[i][j] = imat[j][i];
						}
					}
					//  goto finish;
					gotofinish = true;
					break;
				}

				if (iter == maxiter) {
					break;  /*skip the step halving calc*/
				}

				if (newlk < loglik[1]) {    /*it is not converging ! */
					halving = 1;
					for (i = 0; i < nvar; i++) {
						newbeta[i] = (newbeta[i] + beta[i]) / 2; /*half of old increment */
					}
				} else {
					halving = 0;
					loglik[1] = newlk;
					chsolve2(imat, nvar, u);
					j = 0;
					for (i = 0; i < nvar; i++) {
						beta[i] = newbeta[i];
						newbeta[i] = newbeta[i] + u[i];
					}
				}
			}   /* return for another iteration */
		}

		if (!gotofinish) {
			/*
			 ** We end up here only if we ran out of iterations
			 */
			loglik[1] = newlk;
			chinv2(imat, nvar);
			for (i = 0; i < nvar; i++) {
				beta[i] = newbeta[i] * scale[i];
				u[i] /= scale[i];
				imat[i][i] *= scale[i] * scale[i];
				for (j = 0; j < i; j++) {
					imat[j][i] *= scale[i] * scale[j];
					imat[i][j] = imat[j][i];
				}
			}
			flag = 1000;
		}

//finish:
		/*
		 for (j = 0; j < numCovariates; j++) {
		 b[j] = b[j] / SD[j];
		 * ix = j * (numCovariates + 1) + j
		 SE[j] = Math.sqrt(a[ix(j, j, numCovariates + 1)]) / SD[j];
		 //            o = o + ("   " + variables.get(j) + "    " + Fmt(b[j]) + Fmt(SE[j]) + Fmt(Math.exp(b[j])) + Fmt(Norm(Math.abs(b[j] / SE[j]))) + Fmt(Math.exp(b[j] - 1.95 * SE[j])) + Fmt(Math.exp(b[j] + 1.95 * SE[j])) + NL);
		 CoxCoefficient coe = coxInfo.getCoefficient(variables.get(j));
		 coe.coeff = b[j];
		 coe.stdError = SE[j];
		 coe.hazardRatio = Math.exp(b[j]);
		 coe.pvalue = Norm(Math.abs(b[j] / SE[j]));
		 coe.hazardRatioLoCI = Math.exp(b[j] - 1.95 * SE[j]);
		 coe.hazardRatioHiCI = Math.exp(b[j] + 1.95 * SE[j]);
		 }

		 */

		coxInfo.setScoreLogrankTest(sctest);
		coxInfo.setDegreeFreedom(beta.length);
		coxInfo.setScoreLogrankTestpvalue(ChiSq.chiSq(coxInfo.getScoreLogrankTest(), beta.length));
		coxInfo.setVariance(imat);
		coxInfo.u = u;

		//     for (int n = 0; n < beta.length; n++) {
		//         se[n] = Math.sqrt(imat[n][n]); // / sd[n];
		//     }


		//       System.out.println("coef,se, means,u");
		for (int n = 0; n < beta.length; n++) {
			CoxCoefficient coe = new CoxCoefficient();
			coe.name = variables.get(n);
			coe.mean = means[n];
			coe.standardDeviation = sd[n];
			coe.coeff = beta[n];
			coe.stdError = Math.sqrt(imat[n][n]);
			coe.hazardRatio = Math.exp(coe.getCoeff());
			coe.z = coe.getCoeff() / coe.getStdError();
			coe.pvalue = ChiSq.norm(Math.abs(coe.getCoeff() / coe.getStdError()));
			double z = 1.959964;
			coe.hazardRatioLoCI = Math.exp(coe.getCoeff() - z * coe.getStdError());
			coe.hazardRatioHiCI = Math.exp(coe.getCoeff() + z * coe.getStdError());

			coxInfo.setCoefficient(coe.getName(), coe);
			// System.out.println(beta[n] + "," + se[n] + "," + means[n] + "," + sd[n] + "," + u[n]); //+ "," + imat[n] "," + loglik[n] + "," + sctest[n] + "," + iterations[n] + "," + flag[n]

		}

		coxInfo.maxIterations = maxiter;
		coxInfo.eps = eps;
		coxInfo.toler = toler;

		coxInfo.iterations = iter;
		coxInfo.flag = flag;
		coxInfo.loglikInit = loglik[0];
		coxInfo.loglikFinal = loglik[1];
		coxInfo.method = method;

		//    System.out.println("loglik[0]=" + loglik[0]);
		//    System.out.println("loglik[1]=" + loglik[1]);

		//    System.out.println("chisq? sctest[0]=" + sctest[0]);
		//    System.out.println("?overall model p-value=" + chiSq(sctest[0], beta.length));


		//      System.out.println();
		//       for (int n = 0; n < covar[0].length; n++) {
		//           System.out.print(n);
		//           for (int variable = 0; variable < covar.length; variable++) {
		//               System.out.print("\t" + covar[variable][n]);

		//           }
		//           System.out.println();
		//       }
		//      for (SurvivalInfo si : data) {
		//          System.out.println(si.order + " " + si.getScore());
		//      }
//        coxInfo.dump();


		coxphfitSCleanup(coxInfo, useWeighted, robust,clusterList);
		return coxInfo;
	}

	/**
	 *
	 * @param ci
	 * @param useWeighted
	 * @param robust
	 * @param cluster
	 * @throws Exception
	 */
	public void coxphfitSCleanup(CoxInfo ci, boolean useWeighted,boolean robust, ArrayList<String> cluster) throws Exception {
		//Do cleanup found after coxfit6 is called in coxph.fit.S
		//infs <- abs(coxfit$u %*% var)
		//[ a1 b1] * [a1 b1]
		//           [a2 b2]
		double[][] du = new double[1][ci.u.length];
		du[0] = ci.u;
		double[] infs = Matrix.abs(Matrix.multiply(ci.u, ci.getVariance()));
//        StdArrayIO.print(infs);

		ArrayList<CoxCoefficient> coxCoefficients = new ArrayList<CoxCoefficient>(ci.getCoefficientsList().values());

		for (int i = 0; i < infs.length; i++) {
			double inf = infs[i];
			double coe = coxCoefficients.get(i).getCoeff();
			if (inf > ci.eps && inf > (ci.toler * Math.abs(coe))) {
				ci.message = "Loglik converged before variable ";
			}
		}

		//sum(coef*coxfit$means)
		double sumcoefmeans = 0;
		for (CoxCoefficient cc : coxCoefficients) {
			sumcoefmeans = sumcoefmeans + cc.getCoeff() * cc.getMean();
		}

		// coxph.fit.S line 107
		//lp <- c(x %*% coef) + offset - sum(coef*coxfit$means)
		for (SurvivalInfo si : ci.survivalInfoList) {
			double offset = si.getOffset();
			double lp = 0;
			for (CoxCoefficient cc : coxCoefficients) {
				String name = cc.getName();
				double coef = cc.getCoeff();
				double value = si.getVariable(name);
				lp = lp + value * coef;
			}
			lp = lp + offset - sumcoefmeans;
			si.setLinearPredictor(lp);
			si.setScore(Math.exp(lp));

//           System.out.println("lp score " + si.order + " " + si.time + " " + si.getWeight() + " " + si.getClusterValue() + " " + lp + " " + Math.exp(lp));
		}
//       ci.dump();
		//begin code after call to coxfit6 in coxph.fit.S
		//Compute the martingale residual for a Cox model
		// appears to be C syntax error for = - vs -=
		//(if (nullmodel) in coxph.fit
		double[] res = CoxMart.process(ci.method, ci.survivalInfoList, false);

		for(int i = 0; i < ci.survivalInfoList.size(); i++){
			SurvivalInfo si = ci.survivalInfoList.get(i);
			si.setResidual(res[i]);
		}

		//this represents the end of coxph.fit.S code and we pickup
		//after call to fit <- fitter(X, Y, strats ....) in coxph.R

		if (robust) {
			ci.setNaiveVariance(ci.getVariance());
			double[][] temp;
			double[][] temp0;

			if (cluster != null) {

				temp = ResidualsCoxph.process(ci, ResidualsCoxph.Type.dfbeta, useWeighted, cluster);
				//# get score for null model
				//    if (is.null(init))
				//          fit2$linear.predictors <- 0*fit$linear.predictors
				//    else
				//          fit2$linear.predictors <- c(X %*% init)
				//Set score to 1

				double[] templp = new double[ci.survivalInfoList.size()];
				double[] tempscore = new double[ci.survivalInfoList.size()];
				int i = 0;
				for (SurvivalInfo si : ci.survivalInfoList) {
					templp[i] = si.getLinearPredictor();
					tempscore[i] = si.getScore();
					si.setLinearPredictor(0);
					si.setScore(1.0); //this erases stored value which isn't how the R code does it
					i++;
				}

				temp0 = ResidualsCoxph.process(ci, ResidualsCoxph.Type.score, useWeighted, cluster);

				i = 0;
				for (SurvivalInfo si : ci.survivalInfoList) {
					si.setLinearPredictor(templp[i]);
					si.setScore(tempscore[i]); //this erases stored value which isn't how the R code does it
					i++;
				}


			} else {
				temp = ResidualsCoxph.process(ci, ResidualsCoxph.Type.dfbeta, useWeighted, null);
				//     fit2$linear.predictors <- 0*fit$linear.predictors
				double[] templp = new double[ci.survivalInfoList.size()];
				double[] tempscore = new double[ci.survivalInfoList.size()];
				int i = 0;
				for (SurvivalInfo si : ci.survivalInfoList) {
					templp[i] = si.getLinearPredictor();
					tempscore[i] = si.getScore();
					si.setLinearPredictor(0);
					si.setScore(1.0);
				}
				temp0 = ResidualsCoxph.process(ci, ResidualsCoxph.Type.score, useWeighted, null);

				i = 0;
				for (SurvivalInfo si : ci.survivalInfoList) {
					si.setLinearPredictor(templp[i]);
					si.setScore(tempscore[i]); //this erases stored value which isn't how the R code does it
					i++;
				}
			}
			//fit$var<- t(temp) % * % temp
			double[][] ttemp = Matrix.transpose(temp);
			double[][] var = Matrix.multiply(ttemp, temp);
			ci.setVariance(var);
			//u<- apply(as.matrix(temp0), 2, sum)
			double[] u = new double[temp0[0].length];
			for (int i = 0; i < temp0[0].length; i++) {
				for (int j = 0; j < temp0.length; j++) {
					u[i] = u[i] + temp0[j][i];
				}
			}
			//fit$rscore <- coxph.wtest(t(temp0)%*%temp0, u, control$toler.chol)$test
			double[][] wtemp = Matrix.multiply(Matrix.transpose(temp0),temp0);
			double toler_chol = 1.818989e-12;
		  //  toler_chol = ci.toler;
			WaldTestInfo wti = WaldTest.process(wtemp,u,toler_chol);
			//not giving the correct value
			ci.setRscore(wti.getTest());
		}
		calculateWaldTestInfo(ci);




	}

	static public void calculateWaldTestInfo(CoxInfo ci){
		if(ci.getNumberCoefficients() > 0){
			double toler_chol = 1.818989e-12;
		  //  toler_chol = ci.toler;
			double[][] b = new double[1][ci.getNumberCoefficients()];
			int i = 0;
			for(CoxCoefficient coe : ci.getCoefficientsList().values()){
				b[0][i] = coe.getCoeff();
				i++;
			}
			ci.setWaldTestInfo(WaldTest.process(ci.getVariance(), b, toler_chol));
		}
	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
		CoxR coxr = new CoxR();


		if (true) {
			try {
			   InputStream is = coxr.getClass().getClassLoader().getResourceAsStream("uis-complete.txt");




				WorkSheet worksheet = WorkSheet.readCSV(is, '\t');
				ArrayList<SurvivalInfo> survivalInfoList = new ArrayList<SurvivalInfo>();
				int i = 0;
				for (String row : worksheet.getRows()) {
					double time = worksheet.getCellDouble(row, "TIME");
					double age = worksheet.getCellDouble(row, "AGE");
					double treat = worksheet.getCellDouble(row, "TREAT");
					double c = worksheet.getCellDouble(row, "CENSOR");
					int censor = (int) c;

					SurvivalInfo si = new SurvivalInfo(time, censor);
					si.setOrder(i);
					si.addContinuousVariable("AGE", age);
					si.addContinuousVariable("TREAT", treat);

					survivalInfoList.add(si);
					i++;
				}

				CoxR cox = new CoxR();
				ArrayList<String> variables = new ArrayList<String>();
				//               variables.add("AGE");

				variables.add("AGE");
				variables.add("TREAT");

				//       variables.add("TREAT:AGE");
			  //  ArrayList<Integer> cluster = new ArrayList<Integer>();
				CoxInfo ci = cox.process(variables, survivalInfoList, false, true,false, false);
				System.out.println(ci);
			} catch (Exception e) {
				e.printStackTrace();
			}


		}

//        if (false) {
//
//            try {
//
//
//                WorkSheet worksheet = WorkSheet.readCSV("/Users/Scooter/NetBeansProjects/AssayWorkbench/src/edu/scripps/assayworkbench/cox/uis-complete.txt", '\t');
//                ArrayList<String> rows = worksheet.getRows();
//                ArrayList<String> variables = new ArrayList<String>();
//                variables.add("AGE");
//                variables.add("TREAT");
//                double[] time2 = new double[rows.size()];
//                int[] status2 = new int[rows.size()];
//                double[][] covar2 = new double[variables.size()][rows.size()];
//                double[] offset2 = new double[rows.size()];
//                double[] weights2 = new double[rows.size()];
//                int[] strata2 = new int[rows.size()];
//
//
//                for (int i = 0; i < rows.size(); i++) {
//                    String row = rows.get(i);
//                    double time = worksheet.getCellDouble(row, "TIME");
//                    //      double age = worksheet.getCellDouble(row, "AGE");
//                    //      double treat = worksheet.getCellDouble(row, "TREAT");
//                    double c = worksheet.getCellDouble(row, "CENSOR");
//                    int censor = (int) c;
//
//                    time2[i] = time;
//                    status2[i] = censor;
//                    offset2[i] = 0;
//                    weights2[i] = 1;
//                    strata2[i] = 0;
//
//                    for (int j = 0; j < variables.size(); j++) {
//                        String variable = variables.get(j);
//                        double v = worksheet.getCellDouble(row, variable);
//                        covar2[j][i] = v;
//                    }
//
//
//                }
//                //from coxph.control.S
//                int maxiter2 = 20;
//                double eps2 = 1e-9;
//                double toler2 = Math.pow(eps2, .75);
//                int doscale2 = 1;
//                int method2 = 0;
//                //toler.chol = eps ^ .75
//                //toler.inf=sqrt(eps)
//                //outer.max=10
//
//                CoxR cox = new CoxR();
//                //        cox.coxfit6(maxiter2, time2, status2, covar2, offset2, weights2, strata2, method2, eps2, toler2, time2, doscale2);
//
//
//
//
//
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        }

	}

	/* $Id: chinv2.c 11357 2009-09-04 15:22:46Z therneau $
	 **
	 ** matrix inversion, given the FDF' cholesky decomposition
	 **
	 ** input  **matrix, which contains the chol decomp of an n by n
	 **   matrix in its lower triangle.
	 **
	 ** returned: the upper triangle + diagonal contain (FDF')^{-1}
	 **            below the diagonal will be F inverse
	 **
	 **  Terry Therneau
	 */
	void chinv2(double[][] matrix, int n) {
		double temp;
		int i, j, k;

		/*
		 ** invert the cholesky in the lower triangle
		 **   take full advantage of the cholesky's diagonal of 1's
		 */
		for (i = 0; i < n; i++) {
			if (matrix[i][i] > 0) {
				matrix[i][i] = 1 / matrix[i][i];   /*this line inverts D */
				for (j = (i + 1); j < n; j++) {
					matrix[j][i] = -matrix[j][i];
					for (k = 0; k < i; k++) /*sweep operator */ {
						matrix[j][k] += matrix[j][i] * matrix[i][k];
					}
				}
			}
		}

		/*
		 ** lower triangle now contains inverse of cholesky
		 ** calculate F'DF (inverse of cholesky decomp process) to get inverse
		 **   of original matrix
		 */
		for (i = 0; i < n; i++) {
			if (matrix[i][i] == 0) {  /* singular row */
				for (j = 0; j < i; j++) {
					matrix[j][i] = 0;
				}
				for (j = i; j < n; j++) {
					matrix[i][j] = 0;
				}
			} else {
				for (j = (i + 1); j < n; j++) {
					temp = matrix[j][i] * matrix[j][j];
					if (j != i) {
						matrix[i][j] = temp;
					}
					for (k = i; k < j; k++) {
						matrix[i][k] += temp * matrix[j][k];
					}
				}
			}
		}
	}

	/*  $Id: chsolve2.c 11376 2009-12-14 22:53:57Z therneau $
	 **
	 ** Solve the equation Ab = y, where the cholesky decomposition of A and y
	 **   are the inputs.
	 **
	 ** Input  **matrix, which contains the chol decomp of an n by n
	 **   matrix in its lower triangle.
	 **        y[n] contains the right hand side
	 **
	 **  y is overwriten with b
	 **
	 **  Terry Therneau
	 */
	void chsolve2(double[][] matrix, int n, double[] y) {
		int i, j;
		double temp;

		/*
		 ** solve Fb =y
		 */
		for (i = 0; i < n; i++) {
			temp = y[i];
			for (j = 0; j < i; j++) {
				temp -= y[j] * matrix[i][j];
			}
			y[i] = temp;
		}
		/*
		 ** solve DF'z =b
		 */
		for (i = (n - 1); i >= 0; i--) {
			if (matrix[i][i] == 0) {
				y[i] = 0;
			} else {
				temp = y[i] / matrix[i][i];
				for (j = i + 1; j < n; j++) {
					temp -= y[j] * matrix[j][i];
				}
				y[i] = temp;
			}
		}
	}



	/**
	 *
	 * @param x
	 * @return
	 */
	public double coxsafe(double x) {
		if (x < -200) {
			return -200;
		}
		if (x > 22) {
			return 22;
		}
		return x;
	}
}
