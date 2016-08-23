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
package org.biojava.nbio.survival.cox.comparators;

import org.biojava.nbio.survival.cox.CoxInfo;
import org.biojava.nbio.survival.cox.CoxVariables;

import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class MeanModelComparator implements Comparator<CoxVariables>, Serializable {
    private static final long serialVersionUID = 1;

	String variable = "";

	/**
	 *
	 * @param variable
	 */
	public MeanModelComparator(String variable) {
	  this.variable = variable;
	}

	@Override
	public int compare(CoxVariables coxVariables1, CoxVariables coxVariables2) {
		CoxInfo ci1LTMean = coxVariables1.getCoxInfo("<MEAN");
		CoxInfo ci1GTMean = coxVariables1.getCoxInfo(">MEAN");

		if(ci1LTMean == null || ci1GTMean == null)
			return 0;

		double c1LTpvalue = ci1LTMean.getCoefficient(variable).getPvalue();
		double c1GTpvalue = ci1GTMean.getCoefficient(variable).getPvalue();

		double c1ratio = Math.min(c1LTpvalue, c1GTpvalue) / Math.max(c1LTpvalue, c1GTpvalue);

		CoxInfo ci2LTMean = coxVariables2.getCoxInfo("<MEAN");
		CoxInfo ci2GTMean = coxVariables2.getCoxInfo(">MEAN");

		double c2LTpvalue = ci2LTMean.getCoefficient(variable).getPvalue();
		double c2GTpvalue = ci2GTMean.getCoefficient(variable).getPvalue();

	   double c2ratio = Math.min(c2LTpvalue, c2GTpvalue) / Math.max(c2LTpvalue, c2GTpvalue);

		if (c1ratio > c2ratio) {
			return 1;
		} else if (c1ratio < c2ratio) {
			return -1;
		} else {
			return 0;
		}
		//ascending order
		// return coxVariables1.compareTo(coxVariables2);
	}
}
