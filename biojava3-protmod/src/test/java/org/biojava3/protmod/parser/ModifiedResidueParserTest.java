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
 * Created on Jun 7, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.parser;

import java.io.IOException;

import java.net.URL;

import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModifiedResidue;
import org.biojava3.protmod.ProteinModification;

public class ModifiedResidueParserTest extends TestCase {

	public void testParser() throws IOException {
		System.out.println("===Begin testing on ModifiedResidueParser");
		
		URL fileUrl = ModifiedResidueParserTest.class.getResource("3MVJ.pdb");
		assertNotNull(fileUrl);
		
		PDBFileReader pdbReader = new PDBFileReader();
		Structure struc = pdbReader.getStructure(fileUrl);
		
		ModifiedResidueParser parser = new ModifiedResidueParser();
		List<ModifiedResidue> residues = parser.parse(struc, 
				ProteinModification.getByCategory(
				ModificationCategory.CHEMICAL_MODIFICATION));
		
		for (ModifiedResidue residue : residues) {
			ProteinModification mod = residue.getModification();
			AminoAcid aa = residue.getModifiedAminoAcid();
			Chain chain = aa.getParent();
			System.out.println(mod.getPdbccId()+"\t"+chain.getName()+"\t");
		}

		System.out.println("===End testing on ModifiedResidueParser");
	}
}
