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
 * Created on Jun 8, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.parser;

import java.io.IOException;

import java.net.URL;

import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.PDBFileReader;

import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModifiedCompound;
import org.biojava3.protmod.ProteinModification;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class AttachmentParserTest extends TestCase {

	public void testParser() throws IOException, StructureException {
		System.out.println("===Begin testing on ModifiedResidueParser");
		
		String server = "http://www.pdb.org/pdb/files/";
		URL url;
		String name;
		
//		name = "3HN3"; // NAG
//		System.out.println(name);
//		url = new URL(server+name+".pdb");
//		assertNotNull(url);
//		parserTest(url);
//		
//		name = "1CPO"; // XYS
//		System.out.println(name);
//		url = new URL(server+name+".pdb");
//		assertNotNull(url);
//		parserTest(url);
		
		name = "1AL2"; // MYR
		System.out.println(name);
		url = new URL(server+name+".pdb");
		assertNotNull(url);
		parserTest(url);
		
		System.out.println("===End testing on ModifiedResidueParser");
	}
	
	private void parserTest(URL pdbUrl) throws IOException, StructureException {
		PDBFileReader pdbReader = new PDBFileReader();
		Structure struc = pdbReader.getStructure(pdbUrl);
		
		AttachmentParser parser = new AttachmentParser(0.4);
		
		int nrmodel = struc.nrModels();
		for (int modelnr=0; modelnr<nrmodel; modelnr++) {
			System.out.println("Model "+(modelnr+1));
			
			List<ModifiedCompound> mcs = parser.parse(struc, 
					ProteinModification.getByCategory(ModificationCategory.ATTACHMENT),
					modelnr);
			
			int i=0;
			for (ModifiedCompound mc : mcs) {
				System.out.println("Attachment #"+(++i)+":");
				
				Atom[] atoms = mc.getAtomBonds().get(0);
				
				Group residue = mc.getProteinResidues().get(0);
				Chain chain = residue.getParent();
				System.out.println("\t"+residue.getPDBName()+"\t"+chain.getName()+"\t"
						+residue.getPDBCode()+"\t"+atoms[0].getName());
				
				Group group = mc.getOtherGroups().get(0);
				assertEquals(chain, group.getParent());
				System.out.println("\t"+group.getPDBName()+"\t"+chain.getName()+"\t"
						+group.getPDBCode()+"\t"+atoms[1].getName());
				
				System.out.println("\t"+Calc.getDistance(atoms[0], atoms[1]));
			}
		}
	}
}
