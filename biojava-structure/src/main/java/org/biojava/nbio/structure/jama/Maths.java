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
package org.biojava.nbio.structure.jama;

public class Maths {

	/** sqrt(a^2 + b^2) without under/overflow.
	 * @param a a double
	 * @param b a double
	 * @return a double
	 * */

	public static double hypot(double a, double b) {
		double r;
		if (Math.abs(a) > Math.abs(b)) {
			r = b/a;
			r = Math.abs(a)*Math.sqrt(1+r*r);
		} else if (b != 0) {
			r = a/b;
			r = Math.abs(b)*Math.sqrt(1+r*r);
		} else {
			r = 0.0;
		}
		return r;
	}
}
