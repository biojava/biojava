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
 * Created on Jun 2, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import junit.framework.TestCase;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ProteinModificationRegistryTest extends TestCase {
	
	public void testRegisterModification() {
		List<Component> components = new ArrayList<Component>(2);
		components.add(Component.of("COMP1", ComponentType.AMINOACID));
		components.add(Component.of("COMP2", ComponentType.AMINOACID, true, false));
		
		ModificationLinkage linkage = new ModificationLinkage(components, 0, "ATOM1", 1, "ATOM2");
		
		ModificationCondition condition = new ModificationConditionImpl(components, 
				Collections.singletonList(linkage));
		
		ProteinModification.register("TEST", 
				ModificationCategory.CROSS_LINK_2,
				ModificationOccurrenceType.NATURAL,
				condition)
				.setDescription("TEST")
				.setFormula("TEST")
				.setResidId("TEST")
				.setResidName("TEST")
				.setPsimodId("TEST")
				.setPsimodName("TEST")
				.setSystematicName("TEST");
		assertNotNull(ProteinModification.getById("TEST"));
	}
	
	/**
	 * Test the initialization registry of common protein modifications. 
	 */
	public void testRegisterCommonModification() {		
		Set<ProteinModification> mods = ProteinModification.allModifications();
		assertTrue(mods!=null && !mods.isEmpty());
		
//		System.out.println("There are totally "+mods.size()
//				+" protein modifications registered.");
//		
//		printModifications(mods);
	}
	
	public void testGetBy() {
		ProteinModification mod;
		Set<ProteinModification> mods;

		mod = ProteinModification.getById("0001");
		assertNotNull(mod);

		mods = ProteinModification.getByPdbccId("SEP");
		assertNotNull(mods);

		mods = ProteinModification.getByResidId("AA0076");
		assertNotNull(mods);

		mods = ProteinModification.getByPsimodId("MOD:00110");
		assertNotNull(mods);

		mods = ProteinModification.getByComponent(Component.of("FAD", ComponentType.LIGAND));
		assertNotNull(mods);

		mods = ProteinModification.getByCategory(ModificationCategory.ATTACHMENT);
		assertNotNull(mods);

		mods = ProteinModification.getByOccurrenceType(ModificationOccurrenceType.NATURAL);
		assertNotNull(mods);
	}
	
//	/**
//	 * Print modifications.
//	 * @param mods {@link ProteinModification}s.
//	 */
//	private void printModifications(Set<ProteinModification> mods) {
//		for (ProteinModification mod:mods) {
//			System.out.println(mod);	
//		}
//	}
	
}
