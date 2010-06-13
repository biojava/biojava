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

import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModifiedCompound;
import org.biojava3.protmod.ProteinModification;
import org.biojava3.protmod.TmpAtomCache;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ModifiedResidueParserTest extends TestCase {

	public void testParser() {
		String[] names = new String[] {
			"3MVJ", // SEP, TPO
			"1KZU", // FME
			"1AA6", // CSE
			"1NT0", // AHB
			"1ERM", // BHD
			"1QGW", // LYZ
			"2G66", // HY3, HYP
			"1A39", // PCA
			"1AG7", // CUG, HYP
			"1D5W", // PHD
			"1H9C", // CSP
			"1EUD", // NEP
			"1NSQ", // HIP
			"3LXN", // PTR
			"1ZM2", // DDE
			"1E0Z", // ALY
			"1DM3", // SCY
			"2NPP", // MAA
			"1GK8", // MME, HYP
			"1DOJ", // MEA, TYS
			"1G42", // 2MR
			"2B2U", // DA2, M3L
			"1ALL", // MEN
			"3FMY", // MEQ
			"1E6Y", // MHS, AGM
			"1IV8", // MLY, MLZ
			"1ZTO", // AAR
			"1D7T", // CY3, HYP
			"1D5M", // CLE
			"1XAE", // NFA, C-terminal modification, but occurs in non-terminal residue in 1XAE
			"2H9E", // LPD
			"2BF9", // TYC, error reading PDB file
			"1YYL", // VLM
			"1AEX", // SCH
			"1OMW", // CMT
			"2C0J", // P1L
			"1AA1", // KCX
			"1O5K", // MCL
		};
		
		for (String name : names) {
			System.out.println("===\n"+name);
			try {
				parserTest(name);
			} catch (Exception e){
				e.printStackTrace();
				fail(e.getMessage());
			}
		}
	}
	
	private void parserTest(String pdbId) throws IOException, StructureException {		
		Structure struc = TmpAtomCache.cache.getStructure(pdbId);
		
		ModifiedResidueParser parser = new ModifiedResidueParser();
		
		int nrmodel = struc.nrModels();
		for (int modelnr=0; modelnr<nrmodel; modelnr++) {
			List<ModifiedCompound> residues = parser.parse(struc, 
					ProteinModification.getByCategory(ModificationCategory.CHEMICAL_MODIFICATION),
					modelnr);
			
			System.out.println("Model #"+(modelnr+1));
			
			int i=0;
			for (ModifiedCompound residue : residues) {
				Group g = residue.getProteinResidues().get(0);
				Chain chain = g.getParent();
				System.out.println("Modified Residue #"+(++i));
				System.out.println("\t"+g.getPDBName()+"\t"+chain.getName()+"\t"+g.getPDBCode());
			}
		}
	}
}
