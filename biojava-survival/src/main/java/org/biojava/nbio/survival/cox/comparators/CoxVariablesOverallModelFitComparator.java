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
public class CoxVariablesOverallModelFitComparator implements Comparator<CoxVariables>, Serializable {
    private static final long serialVersionUID = 1;

	String variables = "";

	/**
	 * Variables are stored as a string representation of an ArrayList
	 * [META_GENE] or [trtg, META_GENE] add variables used in cox regression to an array and then do toString.
	 * @param variables
	 */
	public CoxVariablesOverallModelFitComparator(String variables) {
		this.variables = variables;
	}

	@Override
	public int compare(CoxVariables coxVariables1, CoxVariables coxVariables2) {
		CoxInfo ci1 = coxVariables1.getCoxInfo(variables);
		CoxInfo ci2 = coxVariables2.getCoxInfo(variables);

		if (ci1.getWaldTestInfo().getPvalue() < ci2.getWaldTestInfo().getPvalue()) {
			return -1;
		} else if (ci1.getWaldTestInfo().getPvalue() > ci2.getWaldTestInfo().getPvalue()) {
			return 1;
		} else {
			return 0;
		}
		//ascending order
		// return coxVariables1.compareTo(coxVariables2);
	}
}
