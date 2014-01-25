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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure.align.symm.protodomain;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.symm.CeSymm;
import org.biojava.bio.structure.align.symm.protodomain.ResourceList.NameProvider;
import org.junit.Test;

/**
 * A unit test for {@link ResourceList}.
 * @author dmyerstu
 */
public class ResourceListTest {

	@Test
	public void testSequentially() {
		NameProvider provider = NameProvider.defaultNameProvider();
		ResourceList.set(provider, ResourceList.DEFAULT_PDB_DIR);
		StructureAlignment ceSymm = new CeSymm();
		StructureAlignment ce = new CeMain();
		StructureAlignment ceCp = new CeCPMain();
		testMe(ceSymm, "d1b8aa2", "d1b8aa2");
		testMe(ce, "d1md6a_", "d1b8aa2");
		testMe(ceCp, "2a9u.A_39-139", "d1b8aa2");
	}
	
	private AFPChainAndAtoms testMe(StructureAlignment alg, String name1, String name2) {
		NameProvider provider = ResourceList.get().getNameProvider();
		
		// first, we delete the file
		provider.getAlignmentFile(alg.getAlgorithmName(), name1, name2).delete();
		
		// next, we load the alignment
		AFPChainAndAtoms acaa = ResourceList.get().load(alg, name1, name2);
		
		assertEquals("Alignment result name 1 is wrong", acaa.getName1(), name1);
		assertEquals("Alignment result name 2 is wrong", acaa.getName2(), name2);
		assertEquals("Alignment AFPChain name 1 is wrong", acaa.getAfpChain().getName1(), name1);
		assertEquals("Alignment AFPChain name 2 is wrong", acaa.getAfpChain().getName2(), name2);
		
		// the load should have also created the file
		assertTrue("Alignment result file doesn't exist", provider.alignmentExists(alg.getAlgorithmName(), name1, name2));
		// now we load the file again
		// this time it should not save, so check the mod timestamps
		long modBefore = provider.getAlignmentFile(alg.getAlgorithmName(), name1, name2).lastModified();
		ResourceList.get().load(alg, name1, name2);
		long modAfter = provider.getAlignmentFile(alg.getAlgorithmName(), name1, name2).lastModified();
		assertEquals("Existing alignment was recreated", modBefore, modAfter);
		
		return acaa;
	}
	
}
