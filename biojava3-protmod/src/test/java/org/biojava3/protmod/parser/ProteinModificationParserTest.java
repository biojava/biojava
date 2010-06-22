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
import java.util.Set;

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
				"1OGP", // MTQ
				"1EL5", // FAD on CYS
				"1W1O", // FAD on HIS
				"1DII", // FAD on TYR
				"2KJS", "1LK9", // PNS
				"1D7E", // HC4
				"2TMD", // FMN
				"1VAO", // FAD on HIS
				"1PDA", // DPM
				"2J96", // PVN
				"2HIL", // OPE
				
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
				"1A2V", // TPQ
				"1JJU", // TRQ
				"1WCT", // GTH
				"1A2C", // TYS
				"1WCT", // BTR
				"1AUK", // FGL
				"148L", // DAL
				"1C4B", // DIL
				"1T5M", // DSG
				"1CZQ", // DTR
				"2JUE", // DTH
				"1A7Y", // DVA
				"1CXP", // CSO
				"1F8W", // CSX
				"1FFV", // ARO
				"1CKN", // GPL
				"1BUW", // SNC
				"1CZI", // SMC
				"1E93", // OMT
				"1ACD", // CSD
				"1C0T", // CSW
				"1E6Y", // GL3
				"1BI0", // CSS
				"1E6Y", // AGM
				"1HBM", // MGN
				"1FFU", // CSZ
				"3H5R", // SNN, note: SNN is not at C-terminal in some structures, e.g. 3I4W

				// Cross link
				"3M6S", // Disulfide bond
				"1A6L", // F3S
				"1A70", // FES
				"1RPB", // Disulfide bond, and isopeptide (Cys - ASP)
				"3B2M", // isopeptide (Lys - Asn)
				"1CAD", // FE and 4 Cys, cross-link4
				"1FP4", // CFM, HCA, CYS, HIS
				"1M1N", // CFN, HCA, CYS, HIS
				//"1G21", // CFM, HCA, CYS, HIS, (tolerance 0.5)
				//"1M34", // CFM, HCA, CYS, HIS, (tolerance 1.0)
				"1G7K", // CRQ, cross-link1
				"1EMA", // CRO, cross-link1
				"1GGE", // HIS-TYR, cross-link2, (bond length error 0.6)
				"2JE3", // HEC, CYS, CYS, LYS
				//"1MHL", "1MYP" // not work for HEM
				//"3HML", // PQQ, GLU, TYR, (bond length error 2)
				"1FWX","1QNI","2IWF","2IWK", // CU4
				//"1G20", // CLF (bond length error 20)
				"1SU6","1SU7","1SU8","1SUF",
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
//		parser.setbondLengthTolerance(5);
		
		Set<ProteinModification> mods = ProteinModification.getProteinModifications();
//		Set<ProteinModification> mods = ProteinModification.getByCategory(ModificationCategory.ATTACHMENT);
//		Set<ProteinModification> mods = ProteinModification.getByCategory(ModificationCategory.CHEMICAL_MODIFICATION);
//		Set<ProteinModification> mods = ProteinModification.getByResidId("AA0310");
//		Set<ProteinModification> mods = java.util.Collections.singleton(ProteinModification.getById("0102"));
		
		assertFalse(mods.isEmpty());

		int nrmodel = struc.nrModels();
		for (int modelnr=0; modelnr<nrmodel; modelnr++) {
			System.out.println("Model "+(modelnr+1));

			List<ModifiedCompound> mcs = parser.parse(struc, mods, modelnr);

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
		
		if (cat == ModificationCategory.CHEMICAL_MODIFICATION
				|| cat == ModificationCategory.CROSS_LINK_1) {
			Group g = mc.getGroups().get(0);
			Chain chain = g.getParent();
			System.out.println("\t"+g.getPDBName()+"\t"+chain.getName()+"\t"+g.getPDBCode());
		} else {
			List<Atom[]> atomLinkages = mc.getAtomLinkages();
			for (Atom[] atoms : atomLinkages) {
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
