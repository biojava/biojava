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
package org.biojava.nbio.core.alignment.matrices;

import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;


public interface AAIndexProvider {

	 /**
	  * Gets a substitution matrix by its name. The matrices are defined in
	    {@code}src/main/resources/matrices/AAINDEX.txt{@code}
	  * @param matrixName
	  * @return The @{code}SubstitutionMatrix{@code} or null if not exists
	  */
	 SubstitutionMatrix<AminoAcidCompound> getMatrix(String matrixName);
}
