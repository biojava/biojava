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

package org.biojava3.ptm;

import java.util.Set;

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
			System.out.print(mod.id());
			System.out.print("\t"+mod.pdbccId());
			System.out.print("\t"+mod.pdbccName());
			System.out.print("\t"+mod.residId());
			System.out.print("\t"+mod.residName());
			System.out.print("\t"+mod.psimodId());
			System.out.print("\t"+mod.psimodName());
			System.out.print("\t"+mod.description());
			System.out.print("\t"+mod.systematicName());
			System.out.print("\t"+mod.category().label());
			System.out.print("\t"+mod.occurrenceType().label());
			
			String[] comps = mod.components();
			int sizeComp = comps.length;
			assertTrue(comps!=null&&sizeComp>0);
			System.out.print("\t");
			for (int i=0; i<sizeComp; i++) {
				System.out.print(comps[i]+";");
			}
			
			String[][] atoms = mod.atoms();
			assertTrue((atoms==null) || (atoms.length==sizeComp
					&& atoms[0].length==sizeComp));
			System.out.print("\t");
			if (atoms!=null) {
				for (int i=0; i<sizeComp-1; i++) {
					for (int j=i+1; j<sizeComp; j++) {
						if (atoms[i][j]!=null) {
							assertTrue(atoms[j][i]!=null);
							System.out.print(comps[i]+"("+atoms[i][j]+")"
									+"<=>"+comps[j]+"("+atoms[j][i]+");");
						} else {
							assertTrue(atoms[j][i]==null);
						}
					}
				}
			}
			
			System.out.println();
		}
	}
}
