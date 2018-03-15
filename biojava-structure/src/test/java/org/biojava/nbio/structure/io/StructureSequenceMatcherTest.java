/**
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
 * Created by Spencer Bliven
 *
 */
package org.biojava.nbio.structure.io;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Spencer Bliven
 *
 */
public class StructureSequenceMatcherTest {

	private Structure struct1;
	private String[] pdbNum1;
	private String seq1;

	@Before
	public void setUp() throws IOException, StructureException {
		String name1 = "2PTC";

		AtomCache cache = new AtomCache();

		struct1 = cache.getStructure(name1);
		pdbNum1 = new String[] {
				"16", "17", "18", "19", "20", "21", "22", "23", "24", "25",
				"26", "27", "28", "29", "30", "31", "32", "33", "34", "37",
				"38", "39", "40", "41", "42", "43", "44", "45", "46", "47",
				"48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
				"58", "59", "60", "61", "62", "63", "64", "65", "66", "67",
				"69", "70", "71", "72", "73", "74", "75", "76", "77", "78",
				"79", "80", "81", "82", "83", "84", "85", "86", "87", "88",
				"89", "90", "91", "92", "93", "94", "95", "96", "97", "98",
				"99", "100", "101", "102", "103", "104", "105", "106", "107",
				"108", "109", "110", "111", "112", "113", "114", "115", "116",
				"117", "118", "119", "120", "121", "122", "123", "124", "125",
				"127", "128", "129", "130", "132", "133", "134", "135", "136",
				"137", "138", "139", "140", "141", "142", "143", "144", "145",
				"146", "147", "148", "149", "150", "151", "152", "153", "154",
				"155", "156", "157", "158", "159", "160", "161", "162", "163",
				"164", "165", "166", "167", "168", "169", "170", "171", "172",
				"173", "174", "175", "176", "177", "178", "179", "180", "181",
				"182", "183", "184A", "184", "185", "186", "187", "188A", "188",
				"189", "190", "191", "192", "193", "194", "195", "196", "197",
				"198", "199", "200", "201", "202", "203", "204", "209", "210",
				"211", "212", "213", "214", "215", "216", "217", "219", "220",
				"221A", "221", "222", "223", "224", "225", "226", "227", "228",
				"229", "230", "231", "232", "233", "234", "235", "236", "237",
				"238", "239", "240", "241", "242", "243", "244", "245",
				"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
				"13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
				"23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
				"33", "34", "35", "36", "37", "38", "39", "40", "41", "42",
				"43", "44", "45", "46", "47", "48", "49", "50", "51", "52",
				"53", "54", "55", "56", "57", "58"
		};
		seq1 =
				//>2PTC:E|PDBID|CHAIN|SEQUENCE
				"IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNT"+
				"LNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNM"+
				"FCAGYLEGGKDSCQGDSGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN"+
				//>2PTC:I|PDBID|CHAIN|SEQUENCE
				"RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA";

		Assert.assertTrue(seq1.length() == pdbNum1.length);

		/*// report some stats
		System.out.println("The SEQRES and ATOM information is available via the chains:");
		int modelnr = 0 ; // also is 0 if structure is an XRAY structure.
		List<Chain> chains = struct1.getChains(modelnr);
		for (Chain cha:chains){
			List<Group> agr = cha.getAtomGroups(GroupType.AMINOACID);
			List<Group> hgr = cha.getAtomGroups(GroupType.HETATM);
			List<Group> ngr = cha.getAtomGroups(GroupType.NUCLEOTIDE);

			System.out.print("chain: >"+cha.getChainID()+"<");
			System.out.print(" length SEQRES: " +cha.getSeqResLength());
			System.out.print(" length ATOM: " +cha.getAtomLength());
			System.out.print(" aminos: " +agr.size());
			System.out.print(" hetatms: "+hgr.size());
			System.out.println(" nucleotides: "+ngr.size());
		}
		System.out.println(prot.toString());
		*/
	}

	@Test
	public void testSubstructureMatchingProteinSequence() throws CompoundNotFoundException {
		ProteinSequence seq = new ProteinSequence(seq1.substring(30, 40));
		Structure result = StructureSequenceMatcher.getSubstructureMatchingProteinSequence(seq, struct1);

		Assert.assertEquals("Wrong number of groups", 10, StructureTools.getNrGroups(result));
		Assert.assertEquals("Wrong number of chains", 1, result.getChains().size());
		int i = 0;
		for (Group group : result.getChainByIndex(0).getAtomGroups()) {
			Assert.assertTrue("Contains non-amino acid group", group instanceof AminoAcid);
			AminoAcid aa = (AminoAcid) group;
			char c = StructureTools.get1LetterCodeAmino(aa.getPDBName());
			Assert.assertEquals("Wrong amino acid", seq.getSequenceAsString().charAt(i), c);
			i++;
		}
	}

	@Test
	public void testGetProteinSequenceForStructure() {
		Map<Integer,Group> groupIndexPos = new HashMap<Integer,Group>();
		ProteinSequence prot = StructureSequenceMatcher.getProteinSequenceForStructure(struct1, groupIndexPos);


		// Test returned sequence
		Assert.assertEquals("Unreported residues", seq1.length(), prot.getLength());
		Assert.assertEquals("Modified residues", seq1, prot.toString());

		// Test mapping
		Assert.assertEquals("Missing residues in mapping", seq1.length(), groupIndexPos.size());

		for(int res=0;res<seq1.length();res++) {
			Assert.assertTrue("no mapping for group " + res, groupIndexPos.containsKey(res));
			Group g = groupIndexPos.get(res);

			ResidueNumber resnum = g.getResidueNumber();
			Character aa = StructureTools.get1LetterCodeAmino(g.getPDBName());
			Assert.assertEquals("Wrong PDB number at pos " + res, pdbNum1[res], resnum.toString());
			Assert.assertEquals("Wrong Amino acid at pos " + res,
					Character.valueOf(seq1.charAt(res)), aa);
			//System.out.format("%4d %.5s %s\n", res,resnum.toString(),aa.toString());
		}
	}

	@Test
	public void testMatchSequenceToStructure() throws StructureException, CompoundNotFoundException {
		// create modified sequence by removing 10 residues and adding 3
		String sequenceStr = //>2PTC:E|PDBID|CHAIN|SEQUENCE
			"IVGGYTCGAN" +
			"XXX"+ //added
			"TVPYQVSLNS" +
			//"GYHFCGGSLI" +
			"NSQWVVSAAH" +
			"-CYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNT"+
			"LNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNM"+
			"FCAGYLEGGKDSCQGDSGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN";
		String[] correctResidues = new String[] {
				"16", "17", "18", "19", "20", "21", "22", "23", "24", "25",
				null, null, null,
				"26", "27", "28", "29", "30", "31", "32", "33", "34", "37",
				//"38", "39", "40", "41", "42", "43", "44", "45", "46", "47",
				"48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
				null,"58", "59", "60", "61", "62", "63", "64", "65", "66", "67",
				"69", "70", "71", "72", "73", "74", "75", "76", "77", "78",
				"79", "80", "81", "82", "83", "84", "85", "86", "87", "88",
				"89", "90", "91", "92", "93", "94", "95", "96", "97", "98",

				"99", "100", "101", "102", "103", "104", "105", "106", "107",
				"108", "109", "110", "111", "112", "113", "114", "115", "116",
				"117", "118", "119", "120", "121", "122", "123", "124", "125",
				"127", "128", "129", "130", "132", "133", "134", "135", "136",
				"137", "138", "139", "140", "141", "142", "143", "144", "145",
				"146", "147", "148", "149", "150", "151", "152", "153", "154",
				"155", "156", "157", "158", "159", "160", "161", "162", "163",
				"164", "165", "166", "167", "168", "169", "170", "171", "172",
				"173", "174", "175", "176", "177", "178", "179", "180", "181",
				"182", "183", "184A", "184", "185", "186", "187", "188A", "188",
				"189", "190", "191", "192", "193", "194", "195", "196", "197",
				"198", "199", "200", "201", "202", "203", "204", "209", "210",
				"211", "212", "213", "214", "215", "216", "217", "219", "220",
				"221A", "221", "222", "223", "224", "225", "226", "227", "228",
				"229", "230", "231", "232", "233", "234", "235", "236", "237",
				"238", "239", "240", "241", "242", "243", "244", "245"
		};

		System.err.println("Note: the following 10 warnings about missing residues are expected.");
		ProteinSequence seq = new ProteinSequence(sequenceStr);
		ResidueNumber[] match = StructureSequenceMatcher.matchSequenceToStructure(seq, struct1);

		Assert.assertEquals("Wrong length!", sequenceStr.length(), match.length);
		for(int i=0;i<sequenceStr.length();i++) {
			ResidueNumber res = match[i];
			if( res == null) {
				if(!(sequenceStr.charAt(i) == '-' || sequenceStr.charAt(i) == 'X' )) {
					Assert.fail("Incorrectly marked as missing residue at pos " + i + " aa " + sequenceStr.charAt(i));
				}
			} else {
				Group g = struct1.findGroup(res.getChainName(), res.toString());
				Assert.assertNotNull(g);
				String aa3 = g.getPDBName();
				Assert.assertNotNull(aa3);
				Character aa = StructureTools.get1LetterCodeAmino(aa3);
				Assert.assertEquals("Wrong PDB number at position " + i,
						correctResidues[i], g.getResidueNumber().toString());
				Assert.assertEquals("Wrong amino acid at position " + i,
						Character.valueOf(sequenceStr.charAt(i)), aa);
			}
		}
	}

	@Test
	public void testRemoveGaps1() throws CompoundNotFoundException {
		String ungapped = "ACDEFGHIKLMNPQRSTVWY";
		String gapped = "--ACDE-F-GHI..KLM-NPQRSTVWY--";

		ProteinSequence gappedProt = new ProteinSequence(gapped);
		ProteinSequence ungappedProt = StructureSequenceMatcher.removeGaps(gappedProt);

		Assert.assertEquals(ungapped, ungappedProt.getSequenceAsString());
	}

}
