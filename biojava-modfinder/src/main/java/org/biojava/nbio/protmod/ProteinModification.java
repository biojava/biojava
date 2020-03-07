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
 * Created on May 27, 2010
 * Author: Jianjiong Gao
 *
 */

package org.biojava.nbio.protmod;

import java.util.Set;

/**
 * This interface defines information about a specific protein
 * modification.
 *
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ProteinModification {

	/**
	 *
	 * @return modification id.
	 */
    String getId();

	/**
	 *
	 * @return Protein Data Bank Chemical Component ID.
	 */
    String getPdbccId();

	/**
	 *
	 * @return Protein Data Bank Chemical Component name.
	 */
    String getPdbccName();

	/**
	 *
	 * @return RESID ID.
	 */
    String getResidId();

	/**
	 *
	 * @return RESID name.
	 */
    String getResidName();

	/**
	 *
	 * @return PSI-MOD ID.
	 */
    String getPsimodId();

	/**
	 *
	 * @return PSI-MOD name.
	 */
    String getPsimodName();

	/**
	 *
	 * @return Systematic name.
	 */
    String getSystematicName();

	/**
	 *
	 * @return Description.
	 */
    String getDescription();

	/**
	 *
	 * @return a set of keywords.
	 */
    Set<String> getKeywords();

	/**
	 *
	 * @return {@link ModificationCondition}
	 */
    ModificationCondition getCondition();

	/**
	 *
	 * @return formula of the modified residue.
	 */
    String getFormula();

	/**
	 *
	 * @return the modification category.
	 */
    ModificationCategory getCategory();

	/**
	 *
	 * @return the modification occurrence type.
	 */
    ModificationOccurrenceType getOccurrenceType();

}
