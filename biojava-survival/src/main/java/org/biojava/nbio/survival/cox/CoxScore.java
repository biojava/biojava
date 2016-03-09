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

import java.util.ArrayList;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxScore {

	/**
	 *
	 * @param method
	 * @param survivalInfoList
	 * @param coxInfo
	 * @param useStrata
	 * @return
	 */
	public static double[][] process(CoxMethod method, ArrayList<SurvivalInfo> survivalInfoList, CoxInfo coxInfo, boolean useStrata) {
		int i, j, k;
		double temp;
		int n = survivalInfoList.size();

		ArrayList<String> variables = new ArrayList<String>(coxInfo.getCoefficientsList().keySet());
		int nvar = variables.size();

		double deaths;
		int dd;
		double[] time = new double[n];
		double[] status = new double[n];
		double[] strata = new double[n];
		double[] weights = new double[n];
		double[] score = new double[n];
		double[] a = new double[nvar];
		double[] a2 = new double[nvar];
		double denom = 0, e_denom;
		double risk;
		double[][] covar = new double[nvar][n];
		double[][] resid = new double[nvar][n];
		double hazard, meanwt;
		double downwt, temp2;
		double mean;

		//  n = *nx;
		//  nvar  = *nvarx;
		for (int p = 0; p < n; p++) {
			SurvivalInfo si = survivalInfoList.get(p);
			time[p] = si.getTime();
			status[p] = si.getStatus();
			if (useStrata) {
				strata[p] = si.getStrata();
			} else {
				strata[p] = 0;
			}
			weights[p] = si.getWeight();
			score[p] = si.getScore();

			for(int v = 0; v < variables.size(); v++){
				String variable = variables.get(v);
				Double value = si.getVariable(variable);
				covar[v][p] = value;
			}

		}



		//  a = scratch;
		//  a2 = a+nvar;
	/*
		 **  Set up the ragged array
		 */
		//   covar=  dmatrix(covar2, n, nvar);
		//   resid=  dmatrix(resid2, n, nvar);

		e_denom = 0;
		deaths = 0;
		meanwt = 0;
		for (i = 0; i < nvar; i++) {
			a2[i] = 0;
		}
		strata[n - 1] = 1;  /*failsafe */
		for (i = n - 1; i >= 0; i--) {
			if (strata[i] == 1) {
				denom = 0;
				for (j = 0; j < nvar; j++) {
					a[j] = 0;
				}
			}

			risk = score[i] * weights[i];
			denom += risk;
			if (status[i] == 1) {
				deaths++;
				e_denom += risk;
				meanwt += weights[i];
				for (j = 0; j < nvar; j++) {
					a2[j] += risk * covar[j][i];
				}
			}
			for (j = 0; j < nvar; j++) {
				a[j] += risk * covar[j][i];
				resid[j][i] = 0;
			}

			if (deaths > 0 && (i == 0 || strata[i - 1] == 1 || time[i] != time[i - 1])) {
				/* last obs of a set of tied death times */
				if (deaths < 2 || method == CoxMethod.Breslow) {
					hazard = meanwt / denom;
					for (j = 0; j < nvar; j++) {
						temp = (a[j] / denom);     /* xbar */
						for (k = i; k < n; k++) {
							temp2 = covar[j][k] - temp;
							if (time[k] == time[i] && status[k] == 1) {
								resid[j][k] += temp2;
							}
							resid[j][k] -= temp2 * score[k] * hazard;
							if (strata[k] == 1) {
								break;
							}
						}
					}
				} else {  /* the harder case */
					meanwt /= deaths;
					for (dd = 0; dd < deaths; dd++) {
						downwt = dd / deaths;
						temp = denom - downwt * e_denom;
						hazard = meanwt / temp;
						for (j = 0; j < nvar; j++) {
							mean = (a[j] - downwt * a2[j]) / temp;
							for (k = i; k < n; k++) {
								temp2 = covar[j][k] - mean;
								if (time[k] == time[i] && status[k] == 1) {
									resid[j][k] += temp2 / deaths;
									resid[j][k] -= temp2 * score[k] * hazard
											* (1 - downwt);
								} else {
									resid[j][k] -= temp2 * score[k] * hazard;
								}
								if (strata[k] == 1) {
									break;
								}
							}
						}
					}
				}
				e_denom = 0;
				deaths = 0;
				meanwt = 0;
				for (j = 0; j < nvar; j++) {
					a2[j] = 0;
				}
			}
		}

		for (int p = 0; p < n; p++) {
			SurvivalInfo si = survivalInfoList.get(p);
			for (int v = 0; v < variables.size(); v++) {
				si.setResidualVariable(variables.get(v), resid[v][p]);
			}

		}

		//appears to be backward internally
		double[][] flipresid = new double[n][nvar];

		for(int s = 0; s < resid.length; s++){
			for(int t = 0; t  < resid[0].length; t++){
				flipresid[t][s] = resid[s][t];
			}
		}

		return flipresid;

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
	}
}
