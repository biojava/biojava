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

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxVariablesVariableComparator implements CoxComparatorInterface {

	String variables = "";
	String variable = "";

	/**
	 *
	 * @param variables
	 * @param variable
	 */
	public CoxVariablesVariableComparator(String variables, String variable) {
		this.variables = variables;
		this.variable = variable;
		description = "Signatures ranked by model " + variables + " and variable " + variable + " p-value";
	}

	@Override
	public int compare(CoxVariables coxVariables1, CoxVariables coxVariables2) {
		if(coxVariables1.equals(coxVariables2))
			return 0;
		CoxInfo ci1 = coxVariables1.getCoxInfo(variables);
		CoxInfo ci2 = coxVariables2.getCoxInfo(variables);
		if(ci1 == null && ci2 == null)
			return 0;
		if(ci1 == null && ci2 != null)
			return 1;
		if(ci2 == null && ci1 != null)
			return -1;
		if (ci1.getCoefficientsList().get(variable).getPvalue() < ci2.getCoefficientsList().get(variable).getPvalue()) {
			return -1;
		} else if (ci1.getCoefficientsList().get(variable).getPvalue() > ci2.getCoefficientsList().get(variable).getPvalue()) {
			return 1;
		} else {
			return 0;
		}
		//ascending order
		// return coxVariables1.compareTo(coxVariables2);
	}

	@Override
	public String getDescription() {
	   return description;
	}
	String description = "Signatures ranked by p-value ";

	@Override
	public void setDescription(String description) {
	   this.description = description;
	}

	@Override
	public String getModelVariables() {
		return variables;
	}

	@Override
	public String getSortVariable() {
		return variable;
	}



}
