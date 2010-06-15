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
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
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
public class ProteinModificationParserTest extends TestCase {


	public void testMultiParser(){
		String[] names = new String[] {
				// Attachments
				"3HN3", // NAG
				"1CPO", // XYS
				"1AL2", // MYR
				"1L9H", // PLM
				"1BDO", // BTN
				//"2IQD", // no successful test case for LPA
				"1AT9", // RET
				"1DJP", // DO2, (bond length error 3.0)
				"1ALL", // CYC
				"1B8D", // PEB
				
				// Modified resdiues
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
				//"2BF9", // TYC, error reading PDB file
				"1YYL", // VLM
				"1AEX", // SCH
				"1OMW", // CMT
				"2C0J", // P1L
				"1AA1", // KCX
				"1O5K", // MCL
				"1A8I", // LLP
				"2J4Y", // LYR
				//PVL not exist in PDB

				// Cross link
				"3M6S", // Disulfide bond
				"1A6L", // F3S
				"1A70", // FES
				"1RPB", // Disulfide bond, and isopeptide (Cys - ASP)
				"3B2M", // isopeptide (Lys - Asn)
		};
		for ( String name : names){
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

		DefaultProteinModificationParser parser = new DefaultProteinModificationParser();
//		parser.setbondLengthTolerance(20);

		int nrmodel = struc.nrModels();
		for (int modelnr=0; modelnr<nrmodel; modelnr++) {
			System.out.println("Model "+(modelnr+1));

			List<ModifiedCompound> mcs = parser.parse(struc, 
					ProteinModification.getProteinModifications(),
//					ProteinModification.getByCategory(ModificationCategory.ATTACHMENT),
//					ProteinModification.getByCategory(ModificationCategory.CHEMICAL_MODIFICATION),
					modelnr);

			int i=0;
			for (ModifiedCompound mc : mcs) {
				System.out.println("Modification #"+(++i)+":");
				printModification(mc);
			}
		}
	}
	
	private void printModification(ModifiedCompound mc){
		ProteinModification mod = mc.getModification();
		ModificationCategory cat = mod.getCategory();
		System.out.println(cat.label()+": "+mod.getId());
		if (cat == ModificationCategory.CHEMICAL_MODIFICATION) {
			Group g = mc.getGroups().get(0);
			Chain chain = g.getParent();
			System.out.println("\t"+g.getPDBName()+"\t"+chain.getName()+"\t"+g.getPDBCode());
		} else {
			for (Atom[] atoms : mc.getAtomLinkages()) {
				Group group = atoms[0].getParent();
				Chain chain = group.getParent();
				System.out.println("\t"+group.getPDBName()+"\t"+chain.getName()+"\t"
						+group.getPDBCode()+"\t"+atoms[0].getName());

				group = atoms[1].getParent();
				assertEquals(chain, group.getParent());
				System.out.println("\t"+group.getPDBName()+"\t"+chain.getName()+"\t"
						+group.getPDBCode()+"\t"+atoms[1].getName());

				try {
					System.out.println("\t"+Calc.getDistance(atoms[0], atoms[1]));
				} catch (StructureException e) {
					e.printStackTrace();
				}
			}
		}
	}
}
