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
		System.out.println("Testing on initialzation of the registry.");
		Set<ProteinModification> mods = ProteinModification.getProteinModifications();
		assertTrue(mods!=null && !mods.isEmpty());
		
		System.out.println("There are totally "+mods.size()
				+" protein modifications registered.");
		
		for (ProteinModification mod:mods) {
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
}
