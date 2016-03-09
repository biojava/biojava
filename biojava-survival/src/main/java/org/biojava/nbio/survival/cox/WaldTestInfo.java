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

import org.biojava.nbio.survival.cox.stats.ChiSq;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class WaldTestInfo {

	private int df = 0;
	double[][] solve;
	double[] bsum;

	/**
	 *
	 * @return
	 */
	public double getTest() {
		return bsum[0];
	}

	@Override
	public String toString() {
		return "Wald test=" + getTest() + "df=" + df + " p-value=" + getPvalue(); //To change body of generated methods, choose Tools | Templates.
	}

	/**
	 *
	 * @return
	 */
	public double getPvalue() {
		double pvalue = ChiSq.chiSq(getTest(), df);
		return pvalue;
	}

	/**
	 * @return the df
	 */
	public int getDf() {
		return df;
	}

	/**
	 * @param df the df to set
	 */
	public void setDf(int df) {
		this.df = df;
	}
}
