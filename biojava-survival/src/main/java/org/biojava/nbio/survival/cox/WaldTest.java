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

import org.biojava.nbio.survival.cox.stats.Cholesky2;
import org.biojava.nbio.survival.cox.stats.Chsolve2;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class WaldTest {
//coxph_wtest, df=as.integer(nvar),as.integer(ntest),as.double(var),tests= as.double(b),solve= double(nvar*ntest),as.double(toler.chol))
	//coxph_wtest(Sint *nvar2, Sint *ntest, double *var, double *b,double *solve, double *tolerch)

	/**
	 *
	 * @param var
	 * @param b
	 * @param toler_chol
	 * @return
	 */
	public static WaldTestInfo process(double[][] var, double[] b, double toler_chol) {
		double[][] b_ = new double[1][b.length];

		for(int i = 0; i < b.length; i++){
			b_[0][i] = b[i];
		}

		return process(var,b_,toler_chol);

	}

	/**
	 *
	 * @param var
	 * @param b
	 * @param toler_chol
	 * @return
	 */
	public static WaldTestInfo process(double[][] var, double[][] b, double toler_chol) {


		int i = 0;

	//      if(ci.coefficientsList.size() == 1){
	//          double b_ = b[0][i];
	//          double t = (b_ * b_) / var[0][0];
	//          return;
	//      }


		//  double toler_chol = ci.toler;
		int ntest = 1;
		int nvar = b[0].length;
		double sum = 0;
		double[][] solve = new double[ntest][nvar];
		double[] bsum = new double[ntest];

		Cholesky2.process(var, nvar, toler_chol);

		int df = 0;
		for (i = 0; i < nvar; i++) {
			if (var[i][i] > 0) {
				df++;  /* count up the df */
			}
		}

		for (i = 0; i < ntest; i++) {
			for (int j = 0; j < nvar; j++) {
				solve[i][j] = b[i][j];
			}
			Chsolve2.process(var, nvar, solve, i);   /*solve now has b* var-inverse */

			sum = 0;
			for (int j = 0; j < nvar; j++) {
				sum += b[i][j] * solve[i][j];
			}
			bsum[i] = sum;                     /* save the result */
			//b += nvar;    /*move to next column of b */
			// solve += nvar;
		}
		//* nvar2 = df;
		WaldTestInfo waldTestInfo = new WaldTestInfo();

		waldTestInfo.setDf(df);
		waldTestInfo.solve = solve;
		waldTestInfo.bsum = bsum;

		return waldTestInfo;
	}
}

