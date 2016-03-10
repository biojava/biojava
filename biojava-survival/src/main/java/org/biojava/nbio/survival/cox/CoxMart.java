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
/*
 * /* Ported from $Id: coxmart.c 11166 2008-11-24 22:10:34Z therneau $
 *
 ** Compute the martingale residual for a Cox model
 **
 ** Input
 **      n       number of subjects
 **      method  will be ==1 for the Efron method
 **      time    vector of times
 **      status  vector of status values
 **      score   the vector of subject scores, i.e., exp(beta*z)
 **      strata  is =1 for the last obs of a strata
 **      mark    carried forward from the coxfit routine
 **
 ** Output
 **      expected the expected number of events for the subject
 **
 ** The martingale residual is more of a nuisance for the Efron method
 **
 */
package org.biojava.nbio.survival.cox;

import java.util.ArrayList;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxMart {

	/**
	 *
	 * @param method
	 * @param survivalInfoList
	 * @param useStrata
	 * @return
	 */
	static public double[] process(CoxMethod method, ArrayList<SurvivalInfo> survivalInfoList, boolean useStrata) {
		int i, j;
		int lastone;
		int n = survivalInfoList.size();
		double deaths, denom = 0, e_denom = 0;
		double hazard;
		double temp, wtsum;
		double downwt;
		double[] time = new double[n];
		double[] status = new double[n];
		double[] strata = new double[n];
		double[] wt = new double[n];
		double[] score = new double[n];

		double[] expect = new double[survivalInfoList.size()];

		for (int p = 0; p < n; p++) {
			SurvivalInfo si = survivalInfoList.get(p);
			time[p] = si.getTime();
			status[p] = si.getStatus();
			if (useStrata) {
				strata[p] = si.getStrata();
			} else {
				strata[p] = 0;
			}
			wt[p] = si.getWeight();
			score[p] = si.getScore();
		}

		strata[n - 1] = 1;  /*failsafe */
		/* Pass 1-- store the risk denominator in 'expect' */
		for (i = n - 1; i >= 0; i--) {  // Error because of no bounds checking in C it is an error on the get i - 1
		 //   SurvivalInfo si = survivalInfoList.get(i);

			if (strata[i] == 1) {
				denom = 0; //strata[i]
			}
			denom += score[i] * wt[i];      //score[i]*wt[i];
			if (i == 0 || strata[i - 1] == 1 || time[i - 1] != time[i]) //strata[i-1]==1 ||  time[i-1]!=time[i]
			{
			 //   si.setResidual(denom);
				expect[i] = denom;
			} else {
			 //   si.setResidual(0); //expect[i] =0;
				expect[i] = 0;
			}
		}

		/* Pass 2-- now do the work */
		deaths = 0;
		wtsum = 0;
		e_denom = 0;
		hazard = 0;
		lastone = 0;
		for (i = 0; i < n; i++) {
	//         SurvivalInfo si = survivalInfoList.get(i);
	//         SurvivalInfo sip1 = null;
	//         if (i + 1 < n) {
	//             sip1 = survivalInfoList.get(i + 1);
	//         }
	//         if (si.getResidual() != 0) {
	//             denom = si.getResidual();
	//         }
			if(expect[i] != 0)
				denom = expect[i];
	//         si.setResidual(status[i]);//expect[i] = status[i];
			expect[i] = status[i];

			deaths += status[i]; //status[i];
			wtsum += status[i] * wt[i]; // status[i]*wt[i];
			e_denom += score[i] * status[i] * wt[i];//score[i]*status[i] *wt[i];
			if ( strata[i] == 1 || time[i + 1] != time[i]) { //strata[i]==1 ||  time[i+1]!=time[i]
		/*last subject of a set of tied times */
				if (deaths < 2 || method == CoxMethod.Breslow) { //*method==0
					hazard += wtsum / denom;
					for (j = lastone; j <= i; j++) {
					//     SurvivalInfo sj = survivalInfoList.get(j);
					//     double res =  sj.getResidual() - score[j] * hazard;
						expect[j] -= score[j] * hazard;
					//     sj.setResidual(res); //expect[j] -= score[j]*hazard;

					}
				} else {
					temp = hazard;
					wtsum /= deaths;
					for (j = 0; j < deaths; j++) {
						downwt = j / deaths;
						hazard += wtsum / (denom - e_denom * downwt);
						temp += wtsum * (1 - downwt) / (denom - e_denom * downwt);
					}
					for (j = lastone; j <= i; j++) {
					 //   SurvivalInfo sj = survivalInfoList.get(j);
						if (status[j] == 0) {
							//this appears to be an error for = - versus -=
						 //   double res = -score[j] * hazard;
							expect[j] = -score[j] * hazard;
						 //   sj.setResidual(res);//expect[j] = -score[j]*hazard; This appears to be an error of -score vs -=
						} else {
							// double res = sj.getResidual() - score[j] * temp;
							expect[j] -= score[j] * temp;
							//expect[j] -=  score[j]* temp;
						}
					}
				}
				lastone = i + 1;
				deaths = 0;
				wtsum = 0;
				e_denom = 0;
			}
			if (strata[i] == 1) {
				hazard = 0;
			}
		}


		for (j = lastone; j < n; j++) {
			expect[j] -= score[j] * hazard;
		}

		return expect;
	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
	}
}
