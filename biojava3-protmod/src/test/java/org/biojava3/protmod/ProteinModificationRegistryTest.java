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

import java.util.Set;

import org.biojava3.protmod.ProteinModification;

import junit.framework.TestCase;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ProteinModificationRegistryTest extends TestCase {
	
	public void testRegisterModification() {
		ModificationCondition condition
			= new ModificationConditionImpl.Builder(
					new Component[] {
						Component.register("COMP1", ComponentType.AMINOACID),
						Component.register("COMP2", ComponentType.AMINOACID, true, false)
					})
			.addLinkage(0, 1, "ATOM1", "ATOM2")
			.build();
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
		Set<ProteinModification> mods = ProteinModification.getProteinModifications();
		assertTrue(mods!=null && !mods.isEmpty());
		
		System.out.println("There are totally "+mods.size()
				+" protein modifications registered.");
		
		printModifications(mods);
	}
	
	public void testGetBy() {
		ProteinModification mod;
		Set<ProteinModification> mods;
		
		System.out.println("getById");
		mod = ProteinModification.getById("0001");
		assertNotNull(mod);
		System.out.println(mod);

		System.out.println("getByPdbccId");
		mods = ProteinModification.getByPdbccId("SEP");
		assertNotNull(mods);
		System.out.println(mods);

		System.out.println("getByResidId");
		mods = ProteinModification.getByResidId("AA0076");
		assertNotNull(mods);
		System.out.println(mods);

		System.out.println("getByPsimodId");
		mods = ProteinModification.getByPsimodId("MOD:00110");
		assertNotNull(mods);
		System.out.println(mods);

		System.out.println("getByComponent");
		mods = ProteinModification.getByComponent(Component.of("FAD"));
		assertNotNull(mods);
		printModifications(mods);

		System.out.println("getByCategory");
		mods = ProteinModification.getByCategory(ModificationCategory.ATTACHMENT);
		assertNotNull(mods);
		printModifications(mods);

		System.out.println("getByOccurrenceType");
		mods = ProteinModification.getByOccurrenceType(ModificationOccurrenceType.NATURAL);
		assertNotNull(mods);
		printModifications(mods);
	}
	
	/**
	 * Print modifications.
	 * @param mods {@link ProteinModification}s.
	 */
	private void printModifications(Set<ProteinModification> mods) {
		for (ProteinModification mod:mods) {
			System.out.println(mod);	
		}
	}
	
}
