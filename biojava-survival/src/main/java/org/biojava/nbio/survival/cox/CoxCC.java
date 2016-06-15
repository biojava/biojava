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

import org.apache.commons.math.stat.correlation.Covariance;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.biojava.nbio.survival.cox.matrix.Matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxCC {

	/**
	 *
	 * @param ci
	 * @throws Exception
	 */
	static public void process(CoxInfo ci) throws Exception {
		ArrayList<SurvivalInfo> survivalInfoList = ci.survivalInfoList;
		//r
		ArrayList<String> variables = new ArrayList<String>(ci.getCoefficientsList().keySet());

		ArrayList<Integer> strataClass = new ArrayList<Integer>(survivalInfoList.size());
		double[] wt = new double[survivalInfoList.size()];
		for (int i = 0; i < survivalInfoList.size(); i++) {
			SurvivalInfo si = survivalInfoList.get(i);
			strataClass.add(si.getStrata());
			wt[i] = si.getWeight();
		}


		double[][] r = ResidualsCoxph.process(ci, ResidualsCoxph.Type.score, false, null); // dn not use weighted

		// ArrayList<String> variables = ci.survivalInfoList.get(0).getDataVariables();
//        if (false) {
//            for (int i = 0; i < survivalInfoList.size(); i++) {
//                SurvivalInfo si = survivalInfoList.get(i);
//                System.out.print("Cox cc " + si.getOrder());
//                for (int j = 0; j < variables.size(); j++) {
//                    System.out.print(" " + r[i][j]);
//                }
//                System.out.println();
//            }
//        }

		double[][] rvar = null;

		if (ci.getNaiveVariance() != null) {
			rvar = ci.getNaiveVariance();
		} else {
			rvar = ci.getVariance();
		}
		//nj
		LinkedHashMap<Integer, Double> nj = new LinkedHashMap<Integer, Double>();
		Collections.sort(strataClass);
		for (Integer value : strataClass) {
			Double count = nj.get(value);
			if (count == null) {
				count = 0.0;
			}
			count++;
			nj.put(value, count);
		}
		//Nj
		LinkedHashMap<Integer, Double> Nj = new LinkedHashMap<Integer, Double>();
		//N = N + Nj[key];
		double N = 0;
		for (int i = 0; i < survivalInfoList.size(); i++) {
			SurvivalInfo si = survivalInfoList.get(i);
			Integer strata = si.getStrata();
			Double weight = si.getWeight();
			Double sum = Nj.get(strata);
			if (sum == null) {
				sum = 0.0;
			}
			sum = sum + weight;
			Nj.put(strata, sum);

		}

		for(Double value : Nj.values()){
			N = N + value;
		}

		LinkedHashMap<Integer, Double> k1j = new LinkedHashMap<Integer, Double>();
		for (Integer key : nj.keySet()) {
			double _nj = (nj.get(key)); //trying to copy what R is doing on precision
			double _Nj = (Nj.get(key));
			//         System.out.println("nj=" + _nj + " Nj=" + _Nj);
			k1j.put(key, _Nj * ((_Nj / _nj) - 1));
		}

		double[][] V = new double[variables.size()][variables.size()];

		for (Integer i : k1j.keySet()) {
			//          System.out.println("Strata=" + i + " " + k1j.get(i) + " " + Nj.get(i) + " " + nj.get(i));
			if (nj.get(i) > 1) {
				LinkedHashMap<String, DescriptiveStatistics> variableStatsMap = new LinkedHashMap<String, DescriptiveStatistics>();

				for (int p = 0; p < survivalInfoList.size(); p++) {
					SurvivalInfo si = survivalInfoList.get(p);
					if (si.getStrata() != i) {
						continue;
					}
					//              System.out.print(si.order + " ");
					for (int col = 0; col < variables.size(); col++) {
						String v = variables.get(col);
						DescriptiveStatistics ds = variableStatsMap.get(v);
						if (ds == null) {
							ds = new DescriptiveStatistics();
							variableStatsMap.put(v, ds);
						}
						ds.addValue(r[p][col]);
						//                  System.out.print(si.getResidualVariable(v) + "  ");
					}
					//              System.out.println();
				}
				//calculate variance covariance matrix var(r[class==levels(class)[i],],use='comp')
				double[][] var_covar = new double[variables.size()][variables.size()];
				for (int m = 0; m < variables.size(); m++) {
					String var_m = variables.get(m);
					for (int n = 0; n < variables.size(); n++) {
						String var_n = variables.get(n);
						if (m == n) {
							DescriptiveStatistics ds = variableStatsMap.get(var_m);
							var_covar[m][n] = ds.getVariance();
						} else {
							DescriptiveStatistics ds_m = variableStatsMap.get(var_m);
							DescriptiveStatistics ds_n = variableStatsMap.get(var_n);
							Covariance cv = new Covariance();
							double covar = cv.covariance(ds_m.getValues(), ds_n.getValues(), true);
							var_covar[m][n] = covar;
						}
					}
				}
		 //              System.out.println();
		 //              System.out.println("sstrat=" + i);
		 //              StdArrayIO.print(var_covar);

					   V = Matrix.add(V, Matrix.scale(var_covar, k1j.get(i))  );

		 //       for (int m = 0; m < V.length; m++) {
		 //           for (int n = 0; n < V.length; n++) {
		 //               V[m][n] = V[m][n] + (k1j.get(i) * var_covar[m][n]);
		  //
		 //           }
		  //      }
			}
		}
		//     System.out.println("V");
		//     StdArrayIO.print(V);
		//     System.out.println();
		//z$var <- rvar + rvar %*% V %*% rvar # replace variance in z
		double[][] imat1 = Matrix.multiply(rvar, V);
		imat1 = Matrix.multiply(imat1, rvar);
		imat1 = Matrix.add(rvar, imat1);
		//  System.out.println("New var");
		//  StdArrayIO.print(imat1);
		ci.setVariance(imat1);

		//need to update walsh stats for overall model
		CoxR.calculateWaldTestInfo(ci);
		//per Bob/Kathryn email on 4/23/2014 in a weighted model LogRank p-value is no longer valid so should erase it
		ci.setScoreLogrankTest(Double.NaN);
		ci.setScoreLogrankTestpvalue(Double.NaN);
	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
	}
}
