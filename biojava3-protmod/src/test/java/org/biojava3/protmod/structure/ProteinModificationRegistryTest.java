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

package org.biojava3.protmod.structure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ComponentType;
import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModificationConditionImpl;
import org.biojava3.protmod.ModificationLinkage;
import org.biojava3.protmod.ModificationOccurrenceType;
import org.biojava3.protmod.ProteinModification;
import org.biojava3.protmod.ProteinModificationImpl;
import org.biojava3.protmod.ProteinModificationRegistry;

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
		
		ProteinModification mod = new ProteinModificationImpl.Builder("TEST", 
				ModificationCategory.CROSS_LINK_2,
				ModificationOccurrenceType.NATURAL,
				condition)
			.setDescription("TEST")
			.setFormula("TEST")
			.setResidId("TEST")
			.setResidName("TEST")
			.setPsimodId("TEST")
			.setPsimodName("TEST")
			.setSystematicName("TEST")
			.build();
		
		
		ProteinModificationRegistry.register(mod);
		assertNotNull(ProteinModificationRegistry.getById("TEST"));
	}
	
	/**
	 * Test the initialization registry of common protein modifications. 
	 */
	public void testRegisterCommonModification() {		
		Set<ProteinModification> mods = ProteinModificationRegistry.allModifications();
		assertTrue(mods!=null && !mods.isEmpty());
		
//		System.out.println("There are totally "+mods.size()
//				+" protein modifications registered.");
//		
//		printModifications(mods);
	}
	
	public void testGetBy() {
		ProteinModification mod;
		Set<ProteinModification> mods;

		mod = ProteinModificationRegistry.getById("0001");
		assertNotNull(mod);

		mods = ProteinModificationRegistry.getByPdbccId("SEP");
		assertNotNull(mods);

		mods = ProteinModificationRegistry.getByResidId("AA0076");
		assertNotNull(mods);

		mods = ProteinModificationRegistry.getByPsimodId("MOD:00110");
		assertNotNull(mods);

		mods = ProteinModificationRegistry.getByComponent(Component.of("FAD", ComponentType.LIGAND));
		assertNotNull(mods);

		mods = ProteinModificationRegistry.getByCategory(ModificationCategory.ATTACHMENT);
		assertNotNull(mods);

		mods = ProteinModificationRegistry.getByOccurrenceType(ModificationOccurrenceType.NATURAL);
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
