/*
 *                  BioJava development code
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
 * Created on Jun 8, 2007
 *
 */
package org.biojava.bio.structure.align.util;

import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;

public class AlignmentToolsTest extends TestCase {
	
	public void testIsSequential() throws StructureException, IOException {
		AtomCache cache = new AtomCache();
		
		String name1, name2;
		Atom[] ca1, ca2;
		AFPChain afpChain;
		StructureAlignment ce;
		
		
		// CP case
		name1="1QDM.A"; // swaposin
		name2="1NKL"; // saposin
		
		ca1=cache.getAtoms(name1);
		ca2=cache.getAtoms(name2);
		
		ce = StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
		afpChain = ce.align(ca1,ca2);
		
		assertFalse("CeCPMain should give non-sequential alignments (between blocks).",AlignmentTools.isSequentialAlignment(afpChain,false));
		assertFalse("CeCPMain should give non-sequential alignments (within blocks).",AlignmentTools.isSequentialAlignment(afpChain,true));

		// linear case		
		ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
		afpChain = ce.align(ca1,ca2);
		
		assertTrue("CeMain should give sequential alignments (within blocks).",AlignmentTools.isSequentialAlignment(afpChain,true));
		assertTrue("CeMain should give sequential alignments (between blocks).",AlignmentTools.isSequentialAlignment(afpChain,false));

		// now change the block interior a bit
		
		int[][][] optAln = afpChain.getOptAln();
		int tmp;
		tmp = optAln[0][0][2];
		optAln[0][0][2] = optAln[0][0][1];
		optAln[0][0][1] = tmp;
		tmp = optAln[0][1][2];
		optAln[0][1][2] = optAln[0][1][1];
		optAln[0][1][1] = tmp;
		
		assertTrue("Modifying block interior shouldn't effect block sequence.",AlignmentTools.isSequentialAlignment(afpChain,false));
		assertFalse("Modifying block interior should be not sequential.",AlignmentTools.isSequentialAlignment(afpChain,true));

	}
}
