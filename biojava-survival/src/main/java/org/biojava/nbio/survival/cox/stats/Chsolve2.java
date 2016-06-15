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

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class Chsolve2 {

	/**
	 *
	 * @param matrix
	 * @param n
	 * @param y
	 * @param test
	 */
	static public void process(double[][] matrix, int n, double[][] y, int test) {
		int i, j;
		double temp;

		/*
		 ** solve Fb =y
		 */
		for (i = 0; i < n; i++) {
			temp = y[test][i];
			for (j = 0; j < i; j++) {
				temp -= y[test][j] * matrix[i][j];
			}
			y[test][i] = temp;
		}
		/*
		 ** solve DF'z =b
		 */
		for (i = (n - 1); i >= 0; i--) {
			if (matrix[i][i] == 0) {
				y[test][i] = 0;
			} else {
				temp = y[test][i] / matrix[i][i];
				for (j = i + 1; j < n; j++) {
					temp -= y[test][j] * matrix[j][i];
				}
				y[test][i] = temp;
			}
		}

	}
}
