/*
 *                  BioJava development code
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
 * Created on Jun 10, 2007
 * Author: Andreas Prlic
 * Author: Jules Jacobsen
 *
 */
package org.biojava.nbio.structure.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.JournalArticle;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.test.util.StringManipulationTestsHelper;
import org.junit.Before;
import org.junit.Test;

/**
 * Test the {@link PDBFileParser}.
 * 
 * @author Aleix Lafita
 *
 */
public class PDBFileParserTest {

	private static PDBFileParser parser;

	public static final String newline = System.getProperty("line.separator");

	@Before
	public void setUp(){
		parser = new PDBFileParser();
	}

	@Test
	public void test2LetterResidueName() throws IOException {
		// from 1a4w:
		String t =
				"HETATM 2242 NA    NA L 541       5.845 -14.122  30.560  0.88 23.48          NA"+newline+
				"HETATM 2243 NA    NA L 542      18.411 -16.475  38.464  0.88 24.77          NA"+newline+
				"TER                                                                             "+newline;
		BufferedReader br = new BufferedReader(new StringReader(t));
		Structure s = parser.parsePDBFile(br);
		String pdb = s.toPDB();

		assertEquals("two letter residue names are not dealt with correctly! ",t,pdb);


	}

	@Test
	public void testCorrectFloatingPointDisplay() throws IOException {

		// from 1a4w:


		String t =
				"ATOM      1  N   ASP L   1A     11.095  19.341  20.188  1.00 30.14           N"+newline+
				"ATOM      2  CA  ASP L   1A     10.070  18.634  19.379  1.00 28.34           C"+newline+
				"ATOM      3  C   ASP L   1A      9.846  17.102  19.503  1.00 26.08           C"+newline+
				"ATOM      4  O   ASP L   1A      8.744  16.584  19.162  1.00 23.47           O"+newline+
				"ATOM      5  CB  ASP L   1A     10.255  18.858  17.853  1.00 37.55           C"+newline+
				"ATOM      6  CG  ASP L   1A      8.836  19.264  17.401  1.00 42.76           C"+newline+
				"ATOM      7  OD1 ASP L   1A      8.058  19.292  18.400  1.00 44.03           O"+newline+
				"ATOM      8  OD2 ASP L   1A      8.616  19.668  16.244  1.00 46.88           O"+newline+
				"ATOM      9  N   CYS L   1      10.835  16.440  20.113  1.00 23.72           N"+newline+
				"ATOM     10  CA  CYS L   1      10.769  14.970  20.210  1.00 20.89           C"+newline+
				"ATOM     11  C   CYS L   1       9.580  14.524  21.006  1.00 18.64           C"+newline+
				"ATOM     12  O   CYS L   1       9.110  15.220  21.912  1.00 19.03           O"+newline+
				"ATOM     13  CB  CYS L   1      12.117  14.468  20.771  1.00 21.77           C"+newline+
				"ATOM     14  SG  CYS L   1      12.247  14.885  22.538  1.00 20.55           S"+newline+
				"TER                                                                             "+newline;

		BufferedReader br = new BufferedReader(new StringReader(t));

		Structure s = parser.parsePDBFile(br);
		String pdb = s.toPDB();

		assertTrue("the created PDB file does not match the input file", pdb.equals(t));


	}

	@Test
	public void testPDBHeader() throws IOException {

		String t =
				"HEADER    COMPLEX (SERINE PROTEASE/INHIBITORS)    06-FEB-98   1A4W "+newline+
				"TITLE     CRYSTAL STRUCTURES OF THROMBIN WITH THIAZOLE-CONTAINING  "+newline+
				"TITLE    2 INHIBITORS: PROBES OF THE S1' BINDING SITE              "+newline+
				"EXPDTA    X-RAY DIFFRACTION                                        "+newline+
				"AUTHOR    J.H.MATTHEWS,R.KRISHNAN,M.J.COSTANZO,B.E.MARYANOFF,      "+newline+
				"AUTHOR   2 A.TULINSKY                                              "+newline;
				//"REMARK   2 RESOLUTION. 2.00 ANGSTROMS.                             "+newline;

		BufferedReader br = new BufferedReader(new StringReader(t));

		Structure s = parser.parsePDBFile(br);
		String pdb = s.toPDB();

		if (! pdb.equalsIgnoreCase(t)){
			StringManipulationTestsHelper.compareString(t, pdb);
			System.out.println(t);
			System.out.println(pdb);
		}

		// we ignore the case here, since the month FEB is written as Feb, which should be ok...
		assertTrue("the created header does not match the PDB file" ,pdb.equalsIgnoreCase(t));



	}

	@Test
	public void testREMARK200() throws IOException {

		// test that the resolution is only read from REMARK 3 lines
		String w1 =
				"REMARK 200  RESOLUTION RANGE HIGH      (A) : 1.20"+newline+
				"REMARK 200  RESOLUTION RANGE LOW       (A) : 20.00"+
				"REMARK   200 RESOLUTION9.9  ANGSTROMS."; // this line could give wrong resolution info, but it should not be parsed;

		BufferedReader br = new BufferedReader(new StringReader(w1));

		Structure s = parser.parsePDBFile(br);

		float resolution =s.getPDBHeader().getResolution();
		assertEquals(resolution,PDBHeader.DEFAULT_RESOLUTION,0.01);
	}

	@Test
	public void testREMARK3() throws IOException {

		// note that we used to parse resolution from REMARK 2, now from REMARK 3 since it is more complete
		// and more consistent with info in mmCIF
		// taken from 4tnl
		String w2 =
				"REMARK   3  DATA USED IN REFINEMENT.                  "+newline+
				"REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.80 "+newline+
				"REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 34.27"+newline+
				"REMARK   3   MIN(FOBS/SIGMA_FOBS)              : 1.421"+newline+
				"REMARK   3   COMPLETENESS FOR RANGE        (%) : 99.9 "+newline+
				"REMARK   3   NUMBER OF REFLECTIONS             : 58500"+newline+
				"REMARK   3                                            "+newline+
				"REMARK   3  FIT TO DATA USED IN REFINEMENT.           "+newline+
				"REMARK   3   R VALUE     (WORKING + TEST SET) : 0.213 "+newline+
				"REMARK   3   R VALUE            (WORKING SET) : 0.212 "+newline+
				"REMARK   3   FREE R VALUE                     : 0.233 "+newline+
				"REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.236 "+newline+
				"REMARK   3   FREE R VALUE TEST SET COUNT      : 3063  ";
		BufferedReader br = new BufferedReader(new StringReader(w2));

		Structure s = parser.parsePDBFile(br);

		float resolution = s.getPDBHeader().getResolution();
		float rfree = s.getPDBHeader().getRfree();


		assertEquals(resolution,1.8, 0.00001);
		assertEquals(rfree,   0.233, 0.00001);
	}

	//        @Test
	//        public void testREMARK800() {
	//                        // from 1a4w:
	//			String t =
	//                                "REMARK 800                                                                      " + newline +
	//                                "REMARK 800 SITE                                                                 " + newline +
	//                                "REMARK 800 SITE_IDENTIFIER: CAT                                                 " + newline +
	//                                "REMARK 800 EVIDENCE_CODE: UNKNOWN                                               " + newline +
	//                                "REMARK 800 SITE_DESCRIPTION: ACTIVE SITE                                        " + newline +
	//                                "REMARK 800 SITE_IDENTIFIER: AC1                                                 " + newline +
	//                                "REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
	//                                "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE NA H 541                  " + newline +
	//                                "REMARK 800 SITE_IDENTIFIER: AC2                                                 " + newline +
	//                                "REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
	//                                "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE NA H 542                  " + newline +
	//                                "REMARK 800 SITE_IDENTIFIER: AC3                                                 " + newline +
	//                                "REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
	//                                "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE ANS H 373                 " + newline +
	//                                "REMARK 800 SITE_IDENTIFIER: AC4                                                 " + newline +
	//                                "REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
	//                                "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE DAR H 350                 " + newline +
	//                                "REMARK 800 SITE_IDENTIFIER: AC5                                                 " + newline +
	//                                "REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
	//                                "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE 2EP H 375                 " + newline +
	//                                "REMARK 800 SITE_IDENTIFIER: AC6                                                 " + newline +
	//                                "REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
	//                                "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE KTH H 377                 " + newline;
	//
	//			BufferedReader br = new BufferedReader(new StringReader(t));
	//		try {
	//			Structure s = parser.parsePDBFile(br);
	//			String pdb = s.toPDB();
	//                        System.out.println("testREMARK800: " + newline  + pdb);
	//			assertTrue("the created PDB file does not match the input file", pdb.equals(t));
	//		} catch (Exception e){
	//			fail(e.getMessage());
	//		}
	//        }

	@Test
	public void testSITE() throws IOException {
		// from 1a4w:
		String remark800Test =
				//don't add these here - they are present in the PDB file but are added in the structure.toPDB()
				//                                "REMARK 800                                                                      " + newline +
				//                                "REMARK 800 SITE                                                                 " + newline +
				"REMARK 800 SITE_IDENTIFIER: CAT                                                 " + newline +
				"REMARK 800 EVIDENCE_CODE: UNKNOWN                                               " + newline +
				"REMARK 800 SITE_DESCRIPTION: ACTIVE SITE                                        " + newline +
				"REMARK 800 SITE_IDENTIFIER: AC1                                                 " + newline +
				"REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
				"REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE NA H 541                  " + newline +
				"REMARK 800 SITE_IDENTIFIER: AC2                                                 " + newline +
				"REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
				"REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE NA H 542                  " + newline +
				"REMARK 800 SITE_IDENTIFIER: AC3                                                 " + newline +
				"REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
				"REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE ANS H 373                 " + newline +
				"REMARK 800 SITE_IDENTIFIER: AC4                                                 " + newline +
				"REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
				"REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE DAR H 350                 " + newline +
				"REMARK 800 SITE_IDENTIFIER: AC5                                                 " + newline +
				"REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
				"REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE 2EP H 375                 " + newline +
				"REMARK 800 SITE_IDENTIFIER: AC6                                                 " + newline +
				"REMARK 800 EVIDENCE_CODE: SOFTWARE                                              " + newline +
				"REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE KTH H 377                 " + newline;
		String sitesTest =
				"SITE     1 CAT  3 HIS H  57  ASP H 102  SER H 195                               " + newline +
				"SITE     1 AC1  6 ARG H 221A LYS H 224  HOH H 403  HOH H 460                    " + newline +
				"SITE     2 AC1  6 HOH H 464  HOH H 497                                          " + newline +
				"SITE     1 AC2  6 LYS H 169  THR H 172  PHE H 204A HOH H 467                    " + newline +
				"SITE     2 AC2  6 HOH H 515  HOH H 518                                          " + newline +
				"SITE     1 AC3  8 TYR H  60A GLU H  97A LEU H  99  TRP H 215                    " + newline +
				"SITE     2 AC3  8 GLY H 216  GLY H 219  DAR H 350  HOH H 448                    " + newline +
				"SITE     1 AC4  8 ASP H 189  ALA H 190  TRP H 215  GLY H 216                    " + newline +
				"SITE     2 AC4  8 GLY H 219  ANS H 373  2EP H 375  HOH H 471                    " + newline +
				"SITE     1 AC5  4 TRP H  60D LEU H  99  DAR H 350  KTH H 377                    " + newline +
				"SITE     1 AC6  5 HIS H  57  TRP H  60D LYS H  60F 2EP H 375                    " + newline +
				"SITE     2 AC6  5 HOH H 532                                                     " + newline;
		InputStream inStream = this.getClass().getResourceAsStream("/1a4w.pdb");
		Structure s = parser.parsePDBFile(inStream);
		//                        System.out.print(s.getSites());
		Chain chain = new ChainImpl();
		chain.setName("H");
		for (Site site : s.getSites()) {
			//System.out.println("Site: " + site.getSiteID());
			for (Group group : site.getGroups()) {
				//manually add the chain as this is added later once the ATOM recrds are parsed usually
				group.setChain(chain);
				//					System.out.println("    PDBName: " + group.getPDBName());
				//					System.out.println("    PDBCode: " + group.getPDBCode());
				//					System.out.println("    Type: " + group.getType());
				//					System.out.println("    Parent: " + group.getChainName());
			}

		}
		StringBuilder remark800 = new StringBuilder();
		StringBuilder sites = new StringBuilder();

		for (Site site : s.getSites()) {
			remark800.append(site.remark800toPDB());
			sites.append(site.toPDB());
		}

		//                        System.out.println("testSITE: " + newline  + pdb);
		if (!remark800.toString().equals(remark800Test)) {
			//				System.out.println("Expected:");
			//				System.out.println(remark800Test);
			//				System.out.println("Got:");
			//				System.out.println(remark800.toString());
		}
		assertTrue("the created PDB REMARK800 section does not match the input file", remark800.toString().equals(remark800Test));

		if (!sites.toString().equals(sitesTest)) {
			System.out.println("Expected:");
			System.out.println(sitesTest);
			System.out.println("Got:");
			System.out.println(sites.toString());
		}
		assertEquals("the created PDB SITE section does not match the input file", sites.toString(),sitesTest);



	}

	@Test
	public void testMultiLineJRNL() throws IOException {
		//            System.out.println("Testing JRNL record parsing from 3pfk");
		String jrnlString =
				"JRNL        AUTH   P.R.EVANS,G.W.FARRANTS,P.J.HUDSON                            " + newline +
				"JRNL        TITL   PHOSPHOFRUCTOKINASE: STRUCTURE AND CONTROL.                  " + newline +
				"JRNL        REF    PHILOS.TRANS.R.SOC.LONDON,    V. 293    53 1981              " + newline +
				"JRNL        REF  2 SER.B                                                        " + newline +
				"JRNL        REFN                   ISSN 0080-4622                               " + newline +
				"JRNL        PMID   6115424                                                      ";


		BufferedReader br = new BufferedReader(new StringReader(jrnlString));
		Structure s = null;

		s = parser.parsePDBFile(br);

		// String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("293", journalArticle.getVolume());
		assertEquals("53", journalArticle.getStartPage());
		assertEquals(1981, journalArticle.getPublicationDate());
		assertEquals("PHILOS.TRANS.R.SOC.LONDON, SER.B", journalArticle.getJournalName());
	}

	@Test
	public void testIncorrectDateFormatMultiLineJRNL() throws IOException{
		//            System.out.println("Testing JRNL record parsing from 3pfk");
		String jrnlString =
				"JRNL        AUTH   P.R.EVANS,G.W.FARRANTS,P.J.HUDSON                            " + newline +
				"JRNL        TITL   PHOSPHOFRUCTOKINASE: STRUCTURE AND CONTROL.                  " + newline +
				"JRNL        REF    PHILOS.TRANS.R.SOC.LONDON,    V. 293    53 19SE              " + newline +
				"JRNL        REF  2 SER.B                                                        " + newline +
				"JRNL        REFN                   ISSN 0080-4622                               " + newline +
				"JRNL        PMID   6115424                                                      ";


		BufferedReader br = new BufferedReader(new StringReader(jrnlString));
		Structure s = null;
		s = parser.parsePDBFile(br);
		//    String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("293", journalArticle.getVolume());
		assertEquals("53", journalArticle.getStartPage());
		assertEquals(0, journalArticle.getPublicationDate());
		assertEquals("PHILOS.TRANS.R.SOC.LONDON, SER.B", journalArticle.getJournalName());
	}

	@Test
	public void testInvalidFormatREFsectionJRNL() throws IOException{
		//            System.out.println("Testing JRNL record parsing from 3pfk");
		String jrnlString =
				"JRNL        AUTH   P.R.EVANS,G.W.FARRANTS,P.J.HUDSON                            " + newline +
				//            "JRNL        TITL   PHOSPHOFRUCTOKINASE: STRUCTURE AND CONTROL.                  " + newline +
				"JRNL        REF    INTERESTING TIMES                                            " + newline +
				"JRNL        REFN                   ISSN 0080-4622                               " + newline +
				"JRNL        PMID   6115424                                                      ";


		BufferedReader br = new BufferedReader(new StringReader(jrnlString));
		Structure s = null;
		s = parser.parsePDBFile(br);
		// String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("", journalArticle.getVolume());
		assertEquals("", journalArticle.getStartPage());
		assertEquals(0, journalArticle.getPublicationDate());
		assertEquals("", journalArticle.getJournalName());
	}

	@Test
	public void testSecondMultiLineJRNL() throws IOException{
		//            System.out.println("Testing JRNL record parsing from 1gpb");
		String jrnlString =
				"JRNL        AUTH   K.R.ACHARYA,D.I.STUART,K.M.VARVILL,L.N.JOHNSON               " + newline +
				"JRNL        TITL   GLYCOGEN PHOSPHORYLASE B: DESCRIPTION OF THE                 " + newline +
				"JRNL        TITL 2 PROTEIN STRUCTURE                                            " + newline +
				"JRNL        REF    GLYCOGEN PHOSPHORYLASE B:                1 1991              " + newline +
				"JRNL        REF  2 DESCRIPTION OF THE PROTEIN                                   " + newline +
				"JRNL        REF  3 STRUCTURE                                                    " + newline +
				"JRNL        PUBL   WORLD SCIENTIFIC PUBLISHING CO.,SINGAPORE                    " + newline +
				"JRNL        REFN                                                                ";


		BufferedReader br = new BufferedReader(new StringReader(jrnlString));
		Structure s = null;
		s = parser.parsePDBFile(br);
		// String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("", journalArticle.getVolume());
		assertEquals("1", journalArticle.getStartPage());
		assertEquals(1991, journalArticle.getPublicationDate());
		assertEquals("GLYCOGEN PHOSPHORYLASE B: DESCRIPTION OF THE PROTEIN STRUCTURE", journalArticle.getJournalName());
	}

	@Test
	public void testSingleLineJRNL() throws IOException{
		//            System.out.println("Testing JRNL record parsing from 2bln");
		String jrnlString =
				"JRNL        AUTH   G.J.WILLIAMS,S.D.BREAZEALE,C.R.H.RAETZ,J.H.NAISMITH          " + newline +
				"JRNL        TITL   STRUCTURE AND FUNCTION OF BOTH DOMAINS OF ARNA, A            " + newline +
				"JRNL        TITL 2 DUAL FUNCTION DECARBOXYLASE AND A                            " + newline +
				"JRNL        TITL 3 FORMYLTRANSFERASE, INVOLVED IN 4-AMINO-4-DEOXY-L-            " + newline +
				"JRNL        TITL 4 ARABINOSE BIOSYNTHESIS.                                      " + newline +
				"JRNL        REF    J.BIOL.CHEM.                  V. 280 23000 2005              " + newline +
				"JRNL        REFN                   ISSN 0021-9258                               " + newline +
				"JRNL        PMID   15809294                                                     " + newline +
				"JRNL        DOI    10.1074/JBC.M501534200                                       ";


		BufferedReader br = new BufferedReader(new StringReader(jrnlString));
		Structure s = null;
		s = parser.parsePDBFile(br);
		//   String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		//            System.out.println(journalArticle.getRef());
		assertEquals("280", journalArticle.getVolume());
		assertEquals("23000", journalArticle.getStartPage());
		assertEquals(2005, journalArticle.getPublicationDate());
		assertEquals("J.BIOL.CHEM.", journalArticle.getJournalName());
	}

	@Test
	public void testToBePublishedJRNL() throws IOException{
		//            System.out.println("Testing JRNL record parsing from 1i2c");
		String jrnlString =
				"JRNL        AUTH   M.J.THEISEN,S.L.SANDA,S.L.GINELL,C.BENNING,                  " + newline +
				"JRNL        AUTH 2 R.M.GARAVITO                                                 " + newline +
				"JRNL        TITL   CHARACTERIZATION OF THE ACTIVE SITE OF                       " + newline +
				"JRNL        TITL 2 UDP-SULFOQUINOVOSE SYNTHASE: FORMATION OF THE                " + newline +
				"JRNL        TITL 3 SULFONIC ACID PRODUCT IN THE CRYSTALLINE STATE.              " + newline +
				"JRNL        REF    TO BE PUBLISHED                                              " + newline +
				"JRNL        REFN                                                                ";


		BufferedReader br = new BufferedReader(new StringReader(jrnlString));
		Structure s = null;
		s = parser.parsePDBFile(br);

		JournalArticle journalArticle = s.getJournalArticle();
		assertNull(journalArticle.getVolume());
		assertNull(journalArticle.getStartPage());
		assertEquals(0,journalArticle.getPublicationDate());
		assertEquals("TO BE PUBLISHED", journalArticle.getJournalName());
	}

	@Test
	public void test4hhbAcceptedAtomNames() throws IOException, StructureException{

		FileParsingParameters params = new FileParsingParameters();

		String[] acceptedAtoms = {StructureTools.CA_ATOM_NAME, StructureTools.CB_ATOM_NAME};
		params.setAcceptedAtomNames(acceptedAtoms);
		AtomCache cache = new AtomCache();

		cache.setFileParsingParams(params);
		Structure s = cache.getStructure("4HHB");

		//System.out.println(s.toPDB());
		Atom[] ca = StructureTools.getRepresentativeAtomArray(s);
		Atom[] cb = StructureTools.getAtomArray(s, new String[]{StructureTools.CB_ATOM_NAME,});

		// gly does not have cb...
		assertTrue(ca.length > cb.length);
		assertTrue(cb.length > 500);

		Atom[] allwithoutGLY = StructureTools.getAtomArray(s, acceptedAtoms);

		assertEquals(allwithoutGLY.length, (ca.length + cb.length - (ca.length-cb.length)));

	}

	@Test
	public void testCorrectAtomNamePadding() throws IOException {

		// from 1a4w:
		String atomLines =
				"HETATM 2242 NA    NA H 541       5.845 -14.122  30.560  0.88 23.48          NA"+newline+
				"HETATM 2243 NA    NA H 542      18.411 -16.475  38.464  0.88 24.77          NA"+newline+
				"HETATM 2244  C1  QWE H 373      17.735 -17.178  22.612  1.00 26.29           C"+newline+
				"HETATM 2245  C2  QWE H 373      18.543 -17.350  21.462  1.00 26.22           C"+newline+
				"HETATM 2246  C3  QWE H 373      19.877 -17.046  21.558  1.00 26.06           C"+newline+
				"HETATM 2247  C4  QWE H 373      20.401 -16.566  22.766  1.00 26.60           C"+newline+
				"HETATM 2248  C4A QWE H 373      19.613 -16.358  23.893  1.00 27.46           C"+newline+
				"HETATM 2249  C5  QWE H 373      20.084 -15.791  25.046  1.00 28.32           C"+newline+
				"HETATM 2250  C6  QWE H 373      19.241 -15.568  26.114  1.00 27.47           C"+newline+
				"HETATM 2251  C7  QWE H 373      17.893 -15.971  26.028  1.00 27.79           C"+newline+
				"HETATM 2252  C8  QWE H 373      17.370 -16.510  24.879  1.00 26.51           C"+newline+
				"HETATM 2253  C8A QWE H 373      18.205 -16.689  23.809  1.00 26.70           C"+newline+
				"HETATM 2254  N   QWE H 373      21.383 -15.285  25.104  1.00 29.72           N"+newline+
				"HETATM 2255  CM1 QWE H 373      21.898 -14.217  24.214  1.00 30.75           C"+newline+
				"HETATM 2256  CM2 QWE H 373      22.635 -16.213  25.496  1.00 31.28           C"+newline+
				"HETATM 2257  S   QWE H 373      15.909 -17.581  22.417  1.00 26.82           S"+newline+
				"HETATM 2258  O1S QWE H 373      15.988 -18.013  21.010  1.00 27.08           O"+newline+
				"HETATM 2259  O2S QWE H 373      14.998 -18.191  23.402  1.00 25.52           O"+newline+
				"HETATM 2260  N1  QWE H 373      14.893 -15.970  22.485  1.00 26.45           N"+newline+
				"HETATM 2261  CA  QWE H 373      15.193 -15.114  21.296  1.00 26.19           C"+newline+
				"HETATM 2262  C   QWE H 373      16.307 -14.167  21.752  1.00 26.27           C"+newline+
				"HETATM 2263  O   QWE H 373      16.212 -13.781  22.914  1.00 24.70           O"+newline+
				"HETATM 2264  CB  QWE H 373      13.942 -14.354  20.879  1.00 23.76           C"+newline+
				"HETATM 2265  CG  QWE H 373      13.232 -13.686  22.058  1.00 25.01           C"+newline+
				"HETATM 2266  CD  QWE H 373      11.892 -13.133  21.595  1.00 23.01           C"+newline+
				"HETATM 2267  NE  QWE H 373      11.218 -12.647  22.780  1.00 21.45           N"+newline+
				"HETATM 2268  CZ  QWE H 373      11.457 -11.554  23.483  1.00 20.16           C"+newline+
				"HETATM 2269  NH1 QWE H 373      12.359 -10.674  23.045  1.00 20.04           N"+newline+
				"HETATM 2270  NH2 QWE H 373      10.723 -11.336  24.566  1.00 20.85           N"+newline+
				"HETATM 2271  N11 QWE H 373      17.280 -13.692  20.830  1.00 27.91           N"+newline+
				"HETATM 2272  C21 QWE H 373      17.043 -14.033  19.327  1.00 30.31           C"+newline+
				"HETATM 2273  C31 QWE H 373      18.421 -14.472  18.809  1.00 29.73           C"+newline+
				"HETATM 2274  C41 QWE H 373      19.610 -13.417  19.194  1.00 29.33           C"+newline+
				"HETATM 2275  C51 QWE H 373      19.696 -13.083  20.765  1.00 28.57           C"+newline+
				"HETATM 2276  C61 QWE H 373      18.351 -12.680  21.310  1.00 28.50           C"+newline+
				"HETATM 2277  C1' QWE H 373      16.593 -13.014  18.550  1.00 33.08           C"+newline+
				"HETATM 2278  C2' QWE H 373      15.941 -13.425  17.179  1.00 37.01           C"+newline+
				"HETATM 2279  S1  QWE H 373      15.783 -14.904  14.458  1.00 43.62           S"+newline+
				"HETATM 2280  O2  QWE H 373      17.488 -11.922  16.287  1.00 41.38           O"+newline+
				"HETATM 2281  C52 QWE H 373      17.059 -15.583  13.508  1.00 43.74           C"+newline+
				"HETATM 2282  C22 QWE H 373      16.864 -13.556  14.739  1.00 42.63           C"+newline+
				"HETATM 2283 C2'1 QWE H 373      16.825 -12.903  16.107  1.00 40.59           C"+newline+
				"HETATM 2284  C42 QWE H 373      18.146 -14.734  13.451  1.00 43.96           C"+newline+
				"HETATM 2285  N3  QWE H 373      18.049 -13.554  14.106  1.00 43.46           N"+newline+
				"TER                                                                             "+newline;

		BufferedReader br = new BufferedReader(new StringReader(atomLines));

		Structure s = parser.parsePDBFile(br);
		String pdb = s.toPDB();


		assertTrue("the created PDB file does not match the input file", pdb.equals(atomLines));

	}
	
	/**
	 * Test handling of missing Element column. Issue 537 in github.
	 */
	@Test
	public void testMissingElements() throws IOException {
		
		// A two residue structure without Element column
		String missingElement =
				"ATOM      1  N   ASP L   1A     11.095  19.341  20.188  1.00 30.14"+newline+
				"ATOM      2  CA  ASP L   1A     10.070  18.634  19.379  1.00 28.34"+newline+
				"ATOM      3  C   ASP L   1A      9.846  17.102  19.503  1.00 26.08"+newline+
				"ATOM      4  O   ASP L   1A      8.744  16.584  19.162  1.00 23.47"+newline+
				"ATOM      5  CB  ASP L   1A     10.255  18.858  17.853  1.00 37.55"+newline+
				"ATOM      6  CG  ASP L   1A      8.836  19.264  17.401  1.00 42.76"+newline+
				"ATOM      7  OD1 ASP L   1A      8.058  19.292  18.400  1.00 44.03"+newline+
				"ATOM      8  OD2 ASP L   1A      8.616  19.668  16.244  1.00 46.88"+newline+
				"ATOM      9  N   CYS L   1      10.835  16.440  20.113  1.00 23.72"+newline+
				"ATOM     10 CA   CYS L   1      10.769  14.970  20.210  1.00 20.89"+newline+
				"ATOM     11  C   CYS L   1       9.580  14.524  21.006  1.00 18.64"+newline+
				"ATOM     12  O   CYS L   1       9.110  15.220  21.912  1.00 19.03"+newline+
				"ATOM     13  CB  CYS L   1      12.117  14.468  20.771  1.00 21.77"+newline+
				"ATOM     14  SG  CYS L   1      12.247  14.885  22.538  1.00 20.55"+newline+
				"TER                                                               "+newline;
		
		// A two residue structure with empty Element column
		String emptyElement =
				"ATOM      1  N   ASP L   1A     11.095  19.341  20.188  1.00 30.14            "+newline+
				"ATOM      2  CA  ASP L   1A     10.070  18.634  19.379  1.00 28.34            "+newline+
				"ATOM      3  C   ASP L   1A      9.846  17.102  19.503  1.00 26.08            "+newline+
				"ATOM      4  O   ASP L   1A      8.744  16.584  19.162  1.00 23.47            "+newline+
				"ATOM      5  CB  ASP L   1A     10.255  18.858  17.853  1.00 37.55            "+newline+
				"ATOM      6  CG  ASP L   1A      8.836  19.264  17.401  1.00 42.76            "+newline+
				"ATOM      7  OD1 ASP L   1A      8.058  19.292  18.400  1.00 44.03            "+newline+
				"ATOM      8  OD2 ASP L   1A      8.616  19.668  16.244  1.00 46.88            "+newline+
				"ATOM      9  N   CYS L   1      10.835  16.440  20.113  1.00 23.72            "+newline+
				"ATOM     10  CA  CYS L   1      10.769  14.970  20.210  1.00 20.89            "+newline+
				"ATOM     11  C   CYS L   1       9.580  14.524  21.006  1.00 18.64            "+newline+
				"ATOM     12  O   CYS L   1       9.110  15.220  21.912  1.00 19.03            "+newline+
				"ATOM     13  CB  CYS L   1      12.117  14.468  20.771  1.00 21.77            "+newline+
				"ATOM     14  SG  CYS L   1      12.247  14.885  22.538  1.00 20.55            "+newline+
				"TER                                                                             "+newline;
		
		String original =
				"ATOM      1  N   ASP L   1A     11.095  19.341  20.188  1.00 30.14           N"+newline+
				"ATOM      2  CA  ASP L   1A     10.070  18.634  19.379  1.00 28.34           C"+newline+
				"ATOM      3  C   ASP L   1A      9.846  17.102  19.503  1.00 26.08           C"+newline+
				"ATOM      4  O   ASP L   1A      8.744  16.584  19.162  1.00 23.47           O"+newline+
				"ATOM      5  CB  ASP L   1A     10.255  18.858  17.853  1.00 37.55           C"+newline+
				"ATOM      6  CG  ASP L   1A      8.836  19.264  17.401  1.00 42.76           C"+newline+
				"ATOM      7  OD1 ASP L   1A      8.058  19.292  18.400  1.00 44.03           O"+newline+
				"ATOM      8  OD2 ASP L   1A      8.616  19.668  16.244  1.00 46.88           O"+newline+
				"ATOM      9  N   CYS L   1      10.835  16.440  20.113  1.00 23.72           N"+newline+
				"ATOM     10  CA  CYS L   1      10.769  14.970  20.210  1.00 20.89           C"+newline+
				"ATOM     11  C   CYS L   1       9.580  14.524  21.006  1.00 18.64           C"+newline+
				"ATOM     12  O   CYS L   1       9.110  15.220  21.912  1.00 19.03           O"+newline+
				"ATOM     13  CB  CYS L   1      12.117  14.468  20.771  1.00 21.77           C"+newline+
				"ATOM     14  SG  CYS L   1      12.247  14.885  22.538  1.00 20.55           S"+newline+
				"TER                                                                             "+newline;

		
		BufferedReader br = new BufferedReader(new StringReader(missingElement));
		Structure s = parser.parsePDBFile(br);
		String pdb = s.toPDB();
		assertTrue("the Element column has not been filled correctly", pdb.equals(original));
		
		
		br = new BufferedReader(new StringReader(emptyElement));
		s = parser.parsePDBFile(br);
		pdb = s.toPDB();
		assertTrue("the Element column has not been filled correctly", pdb.equals(original));
		
	}
	
	/**
	 * Test the parsing of release and last modified dates.
	 */
	@Test
	public void testDates() throws IOException {
		
		String revisionDates = 
		"REVDAT   5   13-JUL-11 1STP    1       VERSN                                    "+newline+
		"REVDAT   4   24-FEB-09 1STP    1       VERSN                                    " + newline+
		"REVDAT   3   01-APR-03 1STP    1       JRNL                                     " + newline+
		"REVDAT   2   15-OCT-94 1STP    1       AUTHOR                                   " + newline+
		"REVDAT   1   15-OCT-92 1STP    0                                                " + newline;
	
		BufferedReader br = new BufferedReader(new StringReader(revisionDates));
		Structure s = parser.parsePDBFile(br);
		
		// The latest modified date should be 2011
		assertEquals(s.getPDBHeader().getModDate().getYear() + 1900, 2011);
		
		// The release date should be 1992
		assertEquals(s.getPDBHeader().getRelDate().getYear() + 1900, 1992);
	
	}
}
