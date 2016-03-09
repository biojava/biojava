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
 * created at May 25, 2008
 */
package org.biojava.nbio.structure.gui.util;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;

/** To be implemented by JPanels that are part of the GUI to trigger structure aligmnents.
 *
 *
 * @author Andreas Prlic
 * @since 1.7
 *
 */
public interface StructurePairSelector {

	public Structure getStructure1() throws StructureException, IOException;
	public Structure getStructure2() throws StructureException, IOException;

}
