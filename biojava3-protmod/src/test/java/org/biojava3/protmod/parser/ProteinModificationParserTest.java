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
		String[][] names = new String[][] {
				// Attachments
				{"3HN3", "AA0151"}, // NAG
				{"1CPO", "AA0406"}, // XYS
				{"1AL2", "AA0059"}, // MYR
				{"1L9H", "AA0106"}, // PLM
				{"1BDO", "AA0117"}, // BTN
				//{"2IQD", "AA0118"}, // no successful test case for LPA
				{"1AT9", "AA0120"}, // RET
				//{"1DJP", "AA0121"}, // DO2, (bond length error 3.0)
				{"1ALL", "AA0131"}, // CYC
				{"1B8D", "AA0132"}, // PEB
				{"1OGP", "AA0142"}, // MTQ
				{"1EL5", "AA0143"}, // FAD on CYS
				{"1W1O", "AA0144"}, // FAD on HIS
				{"1DII", "AA0145"}, // FAD on TYR
				{"2KJS", "AA0150"}, // PNS
				{"1D7E", "AA0207"}, // HC4
				{"2TMD", "AA0220"}, // FMN
				{"1VAO", "AA0221"}, // FAD on HIS
				{"1PDA", "AA0252"}, // DPM
				{"2J96", "AA0258"}, // PVN
				{"2HIL", "AA0264"}, // OPE
				//{"1RTX", "AA0329"}, // HEM, (bond length error 3.0, much closer to FE)
				{"1FEH", "AA0334"}, // HC1
				//{"2Z6D", "AA0351"}, // FMN, (bond length error 2.0)
				{"1N63", "AA0355"}, // CYS-CUN-MCN
				{"1HXQ", "AA0372"}, // U5P on HIS
				{"1QI9", "AA0395"}, // VO4 on HIS
				{"1XG0", "AA0428"}, // DBV on CYS
				//{"1E9W", "AA0447"}, // TSI on ILE, error when reading
				{"2HIL", "AA0497"}, // OPE on SER
				{"3I3L", "AA0522"}, // FAD on ASP
				
				// Modified resdiues
				{"3MVJ", "AA0037"}, // SEP
				{"3MVJ", "AA0038"}, // TPO
				{"1KZU", "AA0021"}, // FME
				{"1AA6", "AA0022"}, // CSE
				{"1NT0", "AA0026"}, // AHB
				{"1ERM", "AA0027"}, // BHD
				{"1QGW", "AA0028"}, // LYZ
				{"2G66", "AA0029"}, // HY3
				{"2G66", "AA0030"}, // HYP
				{"1A39", "AA0031"}, // PCA
				{"1AG7", "AA0032"}, // CGU
				{"1D5W", "AA0033"}, // PHD
				{"1H9C", "AA0034"}, // CSP
				{"1EUD", "AA0035"}, // NEP
				{"1NSQ", "AA0036"}, // HIP
				{"3LXN", "AA0039"}, // PTR
				{"1ZM2", "AA0040"}, // DDE
				{"1E0Z", "AA0055"}, // ALY
				{"1DM3", "AA0056"}, // SCY
				{"2NPP", "AA0061"}, // MAA
				{"1GK8", "AA0064"}, // MME
				{"1DOJ", "AA0065"}, // MEA
				{"1DOJ", "AA0172"}, // TYS
				{"1G42", "AA0067"}, // 2MR
				{"2B2U", "AA0068"}, // DA2
				{"2B2U", "AA0074"}, // M3L
				{"1ALL", "AA0070"}, // MEN
				{"3FMY", "AA0071"}, // MEQ
				{"1E6Y", "AA0073"}, // MHS
				{"1E6Y", "AA0272"}, // AGM
				{"1IV8", "AA0075"}, // MLY
				{"1IV8", "AA0076"}, // MLZ
				{"1ZTO", "AA0082"}, // AAR
				{"1D7T", "AA0085"}, // CY3
				{"1D5M", "AA0091"}, // CLE
				// {"1XAE", "AA0094"}, // NFA, C-terminal modification, but occurs in non-terminal residue
				{"2H9E", "AA0095"}, // LPD
				//{"2BF9", "AA0099"}, // TYC, error reading PDB file
				{"1YYL", "AA0100"}, // VLM
				{"1AEX", "AA0101"}, // SCH
				{"1OMW", "AA0105"}, // CMT
				{"2C0J", "AA0106"}, // P1L
				{"1AA1", "AA0114"}, // KCX
				{"1O5K", "AA0115"}, // MCL
				{"1A8I", "AA0119"}, // LLP
				{"2J4Y", "AA0120"}, // LYR
				//PVL not exist in PDB
				{"1A2V", "AA0147"}, // TPQ
				{"1JJU", "AA0148"}, // TRQ
				{"1WCT", "AA0155"}, // GTH
				{"1A2C", "AA0172"}, // TYS
				{"1WCT", "AA0179"}, // BTR
				{"1AUK", "AA0185"}, // FGL
				{"148L", "AA0191"}, // DAL
				{"1C4B", "AA0192"}, // DIL
				{"1T5M", "AA0196"}, // DSG
				{"1CZQ", "AA0198"}, // DTR
				{"2JUE", "AA0199"}, // DTH
				{"1A7Y", "AA0200"}, // DVA
				{"1CXP", "AA0205"}, // CSO
				{"1F8W", "AA0205"}, // CSX
				{"1FFV", "AA0215"}, // ARO
				{"1CKN", "AA0228"}, // GPL
				{"1BUW", "AA0230"}, // SNC
				{"1CZI", "AA0234"}, // SMC
				{"1E93", "AA0251"}, // OMT
				{"1ACD", "AA0262"}, // CSD
				{"1C0T", "AA0262"}, // CSW
				{"1E6Y", "AA0265"}, // GL3
				{"1BI0", "AA0269"}, // CSS
				{"1E6Y", "AA0272"}, // AGM
				{"1HBM", "AA0273"}, // MGN
				{"1FFU", "AA0277"}, // CSZ
				{"3H5R", "AA0302"}, // SNN, note: SNN is not at C-terminal in some structures, e.g. 3I4W
				{"1JQ7", "AA0311"}, // DMH
				{"1J6Z", "AA0317"}, // HIC
				{"1B80", "AA0322"}, // HTR
				{"1CWM", "AA0336"}, // IML
				{"1BCK", "AA0337"}, // MLE
				{"1EA7", "AA0361"}, // OSE
				{"1TYS", "AA0363"}, // CXM
				{"1EBV", "AA0364"}, // OAS
				{"2VZK", "AA0423"}, // TH5
				{"2IU4", "AA0431"}, // HIQ
				{"1Y9A", "AA0432"}, // OHS
				{"2IUW", "AA0444"}, // LED
				{"1K83", "AA0449"}, // ILX
				{"2VH3", "AA0458"}, // FGL
				{"2AOC", "AA0464"}, // OLT
				{"1DSR", "AA0478"}, // AHB
				{"1AIQ", "AA0493"}, // CXM
				{"1CF0", "AA0509"}, // IYR
				{"1CTP", "AA0510"}, // TYI
				{"3L4M", "AA0520"}, // 0AF
				{"4ECA", "AA0525"}, // AEI

				// Cross link
				{"3M6S", "AA0025"}, // Disulfide bond
				{"1A6L", "AA0139"}, // F3S
				{"1A70", "AA0137"}, // FES
				{"1RPB", "AA0216"}, // Isopeptide (Cys - ASP)
				{"3B2M", "AA0294"}, // isopeptide (Lys - Asn)
				{"1CAD", "AA0136"}, // FE and 4 Cys, cross-link4
				{"1FP4", "AA0141"}, // CFM, HCA, CYS, HIS
				{"1M1N", "AA0141"}, // CFN, HCA, CYS, HIS
				//{"1G21", "AA0141"}, // CFM, HCA, CYS, HIS, (bond length error 0.5)
				//{"1M34", "AA0141"}, // CFM, HCA, CYS, HIS, (bond length error 1.0)
				{"1G7K", "AA0183"}, // CRQ, cross-link1
				{"1EMA", "AA0183"}, // CRO, cross-link1
				//{"1GGE", "AA0250"}, // HIS-TYR, cross-link2, (bond length error 0.6)
				{"2JE3", "AA0271"}, // HEC, CYS, CYS, LYS
				//{"1MHL", "AA0280"}, // not work for HEM
				//{"1MYP", "AA0280"}, // not work for HEM
				//{"3HML", "AA0283"}, // PQQ, GLU, TYR, (bond length error 2)
				{"1FWX", "AA0298"}, // CU4
				{"1QNI", "AA0298"}, // CU4
				{"2IWF", "AA0298"}, // CU4
				{"2IWK", "AA0298"}, // CU4
				//{"1G20", "AA0300"}, // CLF (bond length error 20)
				{"1SU6", "AA0310"}, // NFS, 5 CYS, HIS
				{"1SU7", "AA0310"}, // NFS, 5 CYS, HIS (looks like 6 CYS are linked)
				//{"1JJU", "AA0313"}, // CYS-TRP, (bond length error 3)
				{"1JJU", "AA0314"}, // CYS-ASP
				{"1JJU", "AA0315"}, // CYS-GLU
				//{"1AJ1", "AA0330"}, // CYS-THR, could not find.
				{"1PXQ", "AA0340"}, // CYS-PHE
				{"1PXQ", "AA0342"}, // CYS-THR
				{"1ITK", "AA0348"}, // MET-TYR-TRP
				//{"1R30", "AA0356"}, // 3 CYS-SF4-SAM (bond length error 0.6)
				{"1R30", "AA0357"}, // 3 CYS-FES-ARG
				// {"1S5L", "AA0366"}, // 2 ASP-3 GLU-HIT-OEC (bond length error 6)
				{"1NGK", "AA0368"}, //TYR-TYR
				{"1YZW", "AA0378"}, // CRU
				{"1XQM", "AA0379"}, // CH6
				{"1UIS", "AA0379"}, // NRQ
				{"2OJK", "AA0380"}, // NYG
				{"2A46", "AA0381"}, // CR7
				{"1YZW", "AA0183"}, // CRU
				{"1XQM", "AA0183"}, // CH6
				{"1UIS", "AA0183"}, // NRQ
				{"2OJK", "AA0183"}, // NYG
				{"2A46", "AA0183"}, // CR7
				{"2AXR", "AA0436"}, // CYS-FAD-HIS
				{"2QH7", "AA0438"}, // 3 CYS-FES-HIS
				//{"2VUM", "AA0451"}, // CYS-TRP (bond length error 2)
				{"3EE4", "AA0490"}, // VAL-TYR
				{"3H8L", "AA0513"}, // CYS-S3H-CYS
		};
		for ( String[] name : names){
			System.out.println("===\n"+name[0]);
			try {
				parserTest(name[0], null);
//				parserTest(name[0], name[1]);
			} catch (Exception e){
				e.printStackTrace();
				fail(e.getMessage());
			}
		}
	}	

	private void parserTest(String pdbId, String residId) throws IOException, StructureException {		
		Structure struc = TmpAtomCache.cache.getStructure(pdbId);

		DefaultProteinModificationParser parser = new DefaultProteinModificationParser();
		parser.setRecordUnidentifiableAtomLinkages(true);
//		parser.setbondLengthTolerance(2);
		
		Set<ProteinModification> mods;
		if (residId==null) {
			mods = ProteinModification.getProteinModifications();
		} else {
			mods = ProteinModification.getByResidId(residId);
		}
		
		assertFalse(mods.isEmpty());

		int nrmodel = struc.nrModels();
		for (int modelnr=0; modelnr<nrmodel; modelnr++) {
			System.out.println("Model "+(modelnr+1));

			parser.parse(struc, mods, modelnr);
			List<ModifiedCompound> mcs = parser.getIdentifiedModifiedCompound();
			List<Atom[]> unidentifiedLinkages = parser.getUnidentifiableAtomLinkages();

			assertFalse(mcs.isEmpty());
			
			int i=0;
			for (ModifiedCompound mc : mcs) {
				System.out.println("Modification #"+(++i)+":");
				printModification(mc);
			}
			
			i = 0;
			for (Atom[] atoms : unidentifiedLinkages) {
				System.out.println("Unidenfied linkage #"+(++i)+":");
				printLinkage(atoms);
			}
		}
	}
	
	private void printModification(ModifiedCompound mc){
		ProteinModification mod = mc.getModification();
		ModificationCategory cat = mod.getCategory();
		System.out.println(cat.label()+": "+mod.getId());
		
		List<Atom[]> atomLinkages = mc.getAtomLinkages();
		if (atomLinkages.isEmpty()) {
			Group g = mc.getGroups().iterator().next();
			Chain chain = g.getParent();
			System.out.println("\t"+g.getPDBName()+"\t"+chain.getName()+"\t"+g.getPDBCode());
		} else {
			
			for (Atom[] atoms : atomLinkages) {
				printLinkage(atoms);
			}
		}
	}
	
	private void printLinkage(Atom[] atoms) {
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
