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
 * Created on May 17, 2010
 * Author: Andreas Prlic 
 *
 */

package demo;

import java.util.List;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.StructureIOFile;

/** An example of how to read MMcif files
 * 
 * @author Andreas Prlic
 * 
 */
public class DemoMMCIFReader
{

	public static void main(String[] args){
		String pdbId = "1A4W";

		StructureIOFile pdbreader = new MMCIFFileReader();

		try {
			pdbreader.setAutoFetch(true);
			Structure s = pdbreader.getStructureById(pdbId);
			
			Chain h = s.getChainByPDB("H");

			List<Group> ligands = h.getAtomLigands();
			
			System.out.println("These ligands have been found in chain " + h.getChainID());
			
			for (Group l:ligands){
				System.out.println(l);
			}
			
			System.out.println("Accessing QWE directly: ");
			Group qwe = h.getGroupByPDB(new ResidueNumber("H",373,null));

			System.out.println(qwe.getChemComp());
			
			System.out.println(h.getSeqResSequence());
			System.out.println(h.getAtomSequence());
			System.out.println(h.getAtomGroups(GroupType.HETATM));
			
			System.out.println("Compounds: " + s.getCompounds());
			
		} catch (Exception e) {
			e.printStackTrace();
		}


	}
}
