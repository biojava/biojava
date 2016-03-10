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
package org.biojava.nbio.structure.io.mmcif;

import org.biojava.nbio.structure.io.mmcif.model.ChemComp;

/** Interface that is implemented by all classes that can provide {@link ChemComp} definitions.
 *
 * @author Andreas Prlic
 * @since 3.0
 */
public interface ChemCompProvider {

	/** Returns a new instance of a chemical component definition.
	 *
	 * @param recordName the ID of the {@link ChemComp}
	 * @return a new {@link ChemComp} definition.
	 */
	ChemComp getChemComp(String recordName) ;

}
