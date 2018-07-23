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
package org.biojava.nbio.structure.align.gui;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;

public class StructureAlignmentDisplay {

	/** Display an AFPChain alignment
	 *
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return a StructureAlignmentJmol instance
	 * @throws StructureException
	 */
	public static StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException {

		if ( ca1.length < 1 || ca2.length < 1){
			throw new StructureException("length of atoms arrays is too short! " + ca1.length + "," + ca2.length);
		}

		Group[] twistedGroups = AlignmentTools.prepareGroupsForDisplay(afpChain, ca1, ca2);

		List<Group> hetatms  = StructureTools.getUnalignedGroups(ca1);
		List<Group> hetatms2 = StructureTools.getUnalignedGroups(ca2);

		return DisplayAFP.display(afpChain, twistedGroups, ca1, ca2, hetatms, hetatms2);

	}

}
