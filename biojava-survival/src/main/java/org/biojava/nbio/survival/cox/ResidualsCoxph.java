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

import java.util.ArrayList;
import java.util.LinkedHashMap;

//import java.util.Collections;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ResidualsCoxph {

	/**
	 *
	 */
	public enum Type {

		/**
		 *
		 */
		dfbeta,
		/**
		 *
		 */
		dfbetas,
		/**
		 *
		 */
		score;
	}

	/**
	 *
	 * @param ci
	 * @param type
	 * @param useWeighted
	 * @param cluster
	 * @return
	 * @throws Exception
	 */
	public static double[][] process(CoxInfo ci, Type type, boolean useWeighted, ArrayList<String> cluster) throws Exception {
		Type otype = type;
		if (type == Type.dfbeta || type == Type.dfbetas) {
			type = Type.score;
			//if missing weighted is a required so never missing
		} //64 2 625 310

		double[][] rr = null;
		if (type == Type.score) {
			rr = CoxScore.process(ci.method, ci.survivalInfoList, ci, false);
		}

		//debug
//        if (false) {
//            for (int i = 0; i < ci.survivalInfoList.size(); i++) {
//                SurvivalInfo si = ci.survivalInfoList.get(i);
//                System.out.print("residuals " + si.getOrder() + " " + si.getClusterValue());
//                for (int j = 0; j < 2; j++) {
//                    System.out.print(" " + rr[i][j]);
//                }
//                System.out.println();
//            }
//        }


		double[][] vv = null;
		if (ci.getNaiveVariance() != null) {
			vv = ci.getNaiveVariance();
		} else {
			vv = ci.getVariance();
		}
		if (otype == Type.dfbeta) {
			//rr <- rr %*% vv
			rr = Matrix.multiply(rr, vv);
		} else if (otype == Type.dfbetas) {
			//rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
			double[][] d1 = Matrix.multiply(rr, vv);
			double[][] d2 = Matrix.diag(Matrix.sqrt(Matrix.oneDivide(Matrix.diag(vv))));
			rr = Matrix.multiply(d1, d2);
		}



		if (useWeighted) {
			double[] weighted = ci.getWeighted();
			rr = Matrix.scale(rr, weighted);
		}
		if (cluster != null && cluster.size() > 0) {
			rr = rowsum(rr, cluster);
		}


		return rr;
	}

	/**
	 * From R in residuals.coxph.S rowsum(rr, collapse)
	 *
	 * @param rr
	 * @param sets
	 * @return
	 */
	private static double[][] rowsum(double[][] rr, ArrayList<String> sets) throws Exception {
		LinkedHashMap<String, Double> sumMap = new LinkedHashMap<String, Double>();
		if (rr.length != sets.size()) {
			throw new Exception("Cluster value for each sample are not of equal length n=" + rr.length + " cluster length=" + sets.size());
		}
		double[][] sum = null;
		for (int j = 0; j < rr[0].length; j++) {
			for (int i = 0; i < sets.size(); i++) {
				String s = sets.get(i);
				Double v = sumMap.get(s); //get in order
				if (v == null) {
					v = 0.0;
				}
				v = v + rr[i][j];
				sumMap.put(s, v);

			}
			if (sum == null) {
				sum = new double[sumMap.size()][rr[0].length];
			}

			ArrayList<String> index = new ArrayList<String>(sumMap.keySet());
			//sorting does seem to make a difference in test cases at the .0000000001
	   //     ArrayList<Integer> in = new ArrayList<Integer>();
	   //     for (String s : index) {
	   //         in.add(Integer.parseInt(s));
	   //     }
	   //     Collections.sort(index);

			for (int m = 0; m < index.size(); m++) {
				String key = index.get(m);
				sum[m][j] = sumMap.get(key);
			}

			sumMap.clear();
		}

		return sum;

	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
	}
}
