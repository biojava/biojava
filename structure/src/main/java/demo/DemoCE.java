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
 * Created on Jan 21, 2010
 *
 */
package demo;


import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

public class DemoCE {

	public static void main(String[] args){
		
		String pdbFilePath="/Users/andreas/WORK/PDB/";
		
		boolean isSplit = true;
		
		String name1 = "4hhb.A";
		String name2 = "4hhb.B";
		
		StructureAlignment algorithm  = new CeMain();
		
		AtomCache cache = new AtomCache(pdbFilePath, isSplit);
				
		Structure structure1 = null;
		Structure structure2 = null;

		try {

			structure1 = cache.getStructure(name1);
			structure2 = cache.getStructure(name2);
			
			Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
			Atom[] ca2 = StructureTools.getAtomCAArray(structure2);
			
			AFPChain afpChain = algorithm.align(ca1,ca2);
			//AFPChain afpChain = fatCat.alignFlexible(ca1,ca2,params);

			afpChain.setName1(name1);
			afpChain.setName2(name2);

			// flexible original results:
			System.out.println(afpChain.toFatcat(ca1,ca2));

						
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
}
