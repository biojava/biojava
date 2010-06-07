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

public class ProteinModificationRegistryTest extends TestCase {
	/**
	 * Test the initialization registry of common protein modifications. 
	 */
	public void testRegisteryInit() {
		System.out.println("===Begin Testing on initialzation of the registry===");
		
		Set<ProteinModification> mods = ProteinModification.getProteinModifications();
		assertTrue(mods!=null && !mods.isEmpty());
		
		System.out.println("There are totally "+mods.size()
				+" protein modifications registered.");
		
		printModifications(mods);
		
		System.out.println("===End Testing on initialzation of the registry===");
	}
	
	public void testGetBy() {
		ProteinModification mod;
		Set<ProteinModification> mods;
		
		System.out.println("===Begin Testing getBy... methods===");
		
		System.out.println("getById");
		mod = ProteinModification.getById("0001");
		assertNotNull(mod);
		printModification(mod);

		System.out.println("getByPdbccId");
		mod = ProteinModification.getByPdbccId("SEP");
		assertNotNull(mod);
		printModification(mod);

		System.out.println("getByResidId");
		mod = ProteinModification.getByResidId("AA0076");
		assertNotNull(mod);
		printModification(mod);

		System.out.println("getByPsimodId");
		mod = ProteinModification.getByPsimodId("MOD:00110");
		assertNotNull(mod);
		printModification(mod);

		System.out.println("getByComponent");
		mods = ProteinModification.getByComponent(Component.of("ACE"));
		assertNotNull(mods);
		printModifications(mods);

		System.out.println("getByComponent");
		mods = ProteinModification.getByComponent(Component.of("ACE"), Component.of("MET",true,false));
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
		
		System.out.println("===End Testing getBy... methods===");
	}
	
	/**
	 * Print modifications.
	 * @param mods {@link ProteinModification}s.
	 */
	private void printModifications(Set<ProteinModification> mods) {
		for (ProteinModification mod:mods) {
			printModification(mod);	
		}
	}
	
	/**
	 * Print a modification.
	 * @param mod a {@link ProteinModification}.
	 */
	private void printModification(ProteinModification mod) {
		System.out.print(mod.getId());
		System.out.print("\t"+mod.getPdbccId());
		System.out.print("\t"+mod.getPdbccName());
		System.out.print("\t"+mod.getResidId());
		System.out.print("\t"+mod.getResidName());
		System.out.print("\t"+mod.getPsimodId());
		System.out.print("\t"+mod.getPsimodName());
		System.out.print("\t"+mod.getDescription());
		System.out.print("\t"+mod.getSystematicName());
		System.out.print("\t"+mod.getCategory().label());
		System.out.print("\t"+mod.getOccurrenceType().label());
		
		ModificationCondition condition = mod.getCondition();
		
		Component[] comps = condition.getComponents();
		int sizeComp = comps.length;
		assertTrue(comps!=null&&sizeComp>0);
		System.out.print("\t");
		for (int i=0; i<sizeComp; i++) {
			Component comp = comps[i];
			String str = comp.getPdbccId();
			str += "["+comp.getType().label()+"]";
			if (comp.isCTerminal()) {
				str += "(C)";
			} else if (comp.isNTerminal()) {
				str += "(N)";
			}				
			System.out.print(str+";");
		}
		
		AtomBond[] bonds = condition.getBonds();
		System.out.print("\t");
		if (bonds!=null) {
			for (AtomBond bond:bonds) {
				String str = bond.getComponent1().getPdbccId();
				str += "("+bond.getAtom1()+")<=>";
				str += bond.getComponent2().getPdbccId();
				str += "("+bond.getAtom2()+")";
				System.out.print(str+";");
			}
		}		
		System.out.println();
	}
}
