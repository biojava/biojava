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
package org.biojava.nbio.survival.cox.stats;

import org.biojava.nbio.survival.cox.CoxInfo;
import org.biojava.nbio.survival.cox.CoxMethod;
import org.biojava.nbio.survival.cox.SurvivalInfo;

import java.util.ArrayList;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class AgScore {

	/**
	 *
	 * @param method
	 * @param survivalInfoList
	 * @param coxInfo
	 * @param useStrata
	 * @return
	 */
	public static double[][] process(CoxMethod method, ArrayList<SurvivalInfo> survivalInfoList, CoxInfo coxInfo, boolean useStrata) {
		int i, k;
		//double temp;
		int n = survivalInfoList.size();

		ArrayList<String> variables = new ArrayList<String>(coxInfo.getCoefficientsList().keySet());
		int nvar = variables.size();


		int dd;

		double[] event = new double[n];
		double[] start = new double[n];
		double[] stop = new double[n];

		double[] strata = new double[n];
		double[] weights = new double[n];
		double[] score = new double[n];

		double[] a = new double[nvar];
		double[] a2 = new double[nvar];
		double[] mean = new double[nvar];
		double[] mh1 = new double[nvar];
		double[] mh2 = new double[nvar];
		double[] mh3 = new double[nvar];

		double denom = 0;
		double time = 0;
		double e_denom = 0;
		double meanwt = 0;
		double deaths = 0;
		double risk;
		double[][] covar = new double[nvar][n];
		double[][] resid = new double[nvar][n];
		double hazard;
		double downwt, temp1, temp2, d2;


		int person = 0;

		//  n = *nx;
		//  nvar  = *nvarx;
		for (int p = 0; p < n; p++) {
			SurvivalInfo si = survivalInfoList.get(p);
			stop[p] = si.getTime();
			event[p] = si.getStatus();
			if (useStrata) {
				strata[p] = si.getStrata();
			} else {
				strata[p] = 0;
			}
			weights[p] = si.getWeight();
			score[p] = si.getScore();

			for (int v = 0; v < variables.size(); v++) {
				String variable = variables.get(v);
				Double value = si.getVariable(variable);
				covar[v][p] = value;
			}

		}

		for (person = 0; person < n;) {
			if (event[person] == 0) {
				person++;
			} else {
				/*
				 ** compute the mean over the risk set, also hazard at this time
				 */
				denom = 0;
				e_denom = 0;
				meanwt = 0;
				deaths = 0;
				for (i = 0; i < nvar; i++) {
					a[i] = 0;
					a2[i] = 0;
				}
				time = stop[person];
				for (k = person; k < n; k++) {
					if (start[k] < time) {
						risk = score[k] * weights[k];
						denom += risk;
						for (i = 0; i < nvar; i++) {
							a[i] = a[i] + risk * covar[i][k];
						}
						if (stop[k] == time && event[k] == 1) {
							deaths++;
							e_denom += risk;
							meanwt += weights[k];
							for (i = 0; i < nvar; i++) {
								a2[i] = a2[i] + risk * covar[i][k];
							}
						}
					}
					if (strata[k] == 1) {
						break;
					}
				}

				/* add things in for everyone in the risk set*/
				if (deaths < 2 || method == CoxMethod.Breslow) {
					/* easier case */
					hazard = meanwt / denom;
					for (i = 0; i < nvar; i++) {
						mean[i] = a[i] / denom;
					}
					for (k = person; k < n; k++) {
						if (start[k] < time) {
							risk = score[k];
							for (i = 0; i < nvar; i++) {
								resid[i][k] -= (covar[i][k] - mean[i]) * risk * hazard;
							}
							if (stop[k] == time) {
								person++;
								if (event[k] == 1) {
									for (i = 0; i < nvar; i++) {
										resid[i][k] += (covar[i][k] - mean[i]);
									}
								}
							}
						}
						if (strata[k] == 1) {
							break;
						}
					}
				} else {
					/*
					 ** If there are 3 deaths, let m1, m2, m3 be the three
					 **   weighted means,  h1, h2, h3 be the three hazard jumps.
					 ** Then temp1 = h1 + h2 + h3
					 **      temp2 = h1 + (2/3)h2 + (1/3)h3
					 **      mh1   = m1*h1 + m2*h2 + m3*h3
					 **      mh2   = m1*h1 + (2/3)m2*h2 + (1/3)m3*h3
					 **      mh3   = (1/3)*(m1+m2+m3)
					 */
					temp1 = 0;
					temp2 = 0;
					for (i = 0; i < nvar; i++) {
						mh1[i] = 0;
						mh2[i] = 0;
						mh3[i] = 0;
					}
					meanwt /= deaths;
					for (dd = 0; dd < deaths; dd++) {
						downwt = dd / deaths;
						d2 = denom - downwt * e_denom;
						hazard = meanwt / d2;
						temp1 += hazard;
						temp2 += (1 - downwt) * hazard;
						for (i = 0; i < nvar; i++) {
							mean[i] = (a[i] - downwt * a2[i]) / d2;
							mh1[i] += mean[i] * hazard;
							mh2[i] += mean[i] * (1 - downwt) * hazard;
							mh3[i] += mean[i] / deaths;
						}
					}
					for (k = person; k < n; k++) {
						if (start[k] < time) {
							risk = score[k];
							if (stop[k] == time && event[k] == 1) {
								for (i = 0; i < nvar; i++) {
									resid[i][k] += covar[i][k] - mh3[i];
									resid[i][k] -= risk * covar[i][k] * temp2;
									resid[i][k] += risk * mh2[i];
								}
							} else {
								for (i = 0; i < nvar; i++) {
									resid[i][k] -= risk * (covar[i][k] * temp1 - mh1[i]);
								}
							}
						}
						if (strata[k] == 1) {
							break;
						}
					}
					for (; stop[person] == time; person++) {
						if (strata[person] == 1) {
							break;
						}
					}
				}
			}
		}


		//appears to be backward internally
		double[][] flipresid = new double[n][nvar];

		for (int s = 0; s < resid.length; s++) {
			for (int t = 0; t < resid[0].length; t++) {
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
