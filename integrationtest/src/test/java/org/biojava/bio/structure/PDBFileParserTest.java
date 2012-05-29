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
package org.biojava.bio.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileParser;
import org.biojava.bio.structure.util.StringManipulationTestsHelper;


import junit.framework.TestCase;

public class PDBFileParserTest extends TestCase {

	PDBFileParser parser;

	public static final String newline = System.getProperty("line.separator");

	protected void setUp(){
		parser = new PDBFileParser();
	}




	/** parse the remark lines and return the resolution
	 *
	 * @param fakeFile
	 * @return the resolution as a Float or null if no resolution found
	 * @throws Exception
	 */
	@SuppressWarnings("deprecation")
	private Object testREMARK2Parsing(String fakeFile) throws Exception{
		BufferedReader br = new BufferedReader(new StringReader(fakeFile));

		Object resolution = null;

		Structure s = parser.parsePDBFile(br);
		Map<String, Object> m = s.getHeader();
		resolution =  m.get("resolution");

		return resolution;
	}

	public void test2LetterResidueName() {
		try {

			// from 1a4w:
			String t =
					"HETATM 2242 NA    NA   541       5.845 -14.122  30.560  0.88 23.48          NA"+newline+
					"HETATM 2243 NA    NA   542      18.411 -16.475  38.464  0.88 24.77          NA"+newline;
			BufferedReader br = new BufferedReader(new StringReader(t));
			Structure s = parser.parsePDBFile(br);
			String pdb = s.toPDB();

			assertEquals("two letter residue names are not dealt with correctly! ",t,pdb);


		} catch (Exception e){
			fail(e.getMessage());
		}

	}

	public void testCorrectFloatingPointDisplay() {

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
				"ATOM     14  SG  CYS L   1      12.247  14.885  22.538  1.00 20.55           S"+newline;

		BufferedReader br = new BufferedReader(new StringReader(t));
		try {
			Structure s = parser.parsePDBFile(br);
			String pdb = s.toPDB();

			assertTrue("the created PDB file does not match the input file", pdb.equals(t));
		} catch (Exception e){
			fail(e.getMessage());
		}

	}

	public void testPDBHeader(){

		String t =
				"HEADER    COMPLEX (SERINE PROTEASE/INHIBITORS)    06-FEB-98   1A4W "+newline+
				"TITLE     CRYSTAL STRUCTURES OF THROMBIN WITH THIAZOLE-CONTAINING  "+newline+
				"TITLE    2 INHIBITORS: PROBES OF THE S1' BINDING SITE              "+newline+
				"EXPDTA    X-RAY DIFFRACTION                                        "+newline+
				"AUTHOR    J.H.MATTHEWS,R.KRISHNAN,M.J.COSTANZO,B.E.MARYANOFF,      "+newline+
				"AUTHOR   2 A.TULINSKY                                              "+newline+
				"REMARK   2 RESOLUTION. 2.00 ANGSTROMS.                             "+newline;

		BufferedReader br = new BufferedReader(new StringReader(t));
		try {
			Structure s = parser.parsePDBFile(br);
			String pdb = s.toPDB();

			if (! pdb.equalsIgnoreCase(t)){
				StringManipulationTestsHelper.compareString(t, pdb);
				System.out.println(t);
				System.out.println(pdb);
			}

			// we ignore the case here, since the month FEB is written as Feb, which should be ok...
			assertTrue("the created header does not match the PDB file" ,pdb.equalsIgnoreCase(t));

		} catch (Exception e){
			fail(e.getMessage());
		}

	}

	public void testREMARK200() {

		// test that the resolution is only read from REMARK 2 lines
		String w1 = "REMARK 200  RESOLUTION RANGE HIGH      (A) : 1.20"+newline+
				"REMARK 200  RESOLUTION RANGE LOW       (A) : 20.00"+
				"REMARK   200 RESOLUTION9.9  ANGSTROMS."; // this line could give wrong resolution info, but it should not be parsed;
		boolean parsingOK = true;
		String errorMsg   = "";

		try {
			Object resolution = testREMARK2Parsing(w1);
			assertEquals(resolution,null);
		} catch (Exception e){
			parsingOK = false;
			//e.printStackTrace();
			errorMsg = e.getMessage();
		}


		assertEquals("parsing failed with error " + errorMsg, parsingOK, true);
	}

	public void testREMARK2(){

		String w2 =
				"REMARK   2 "+newline+
				"REMARK   2 RESOLUTION. 1.2  ANGSTROMS."; // the correct line

		boolean parsingOK = true;
		String errorMsg   = "";
		try {
			Object resolution = testREMARK2Parsing(w2);
			assertEquals(resolution,new Float(1.2));
		} catch (Exception e){
			parsingOK = false;
			//e.printStackTrace();
			errorMsg = e.getMessage();
		}


		assertEquals("parsing failed with error " + errorMsg, parsingOK, true);
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

	public void testSITE() {
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
		try {
			Structure s = parser.parsePDBFile(inStream);
			//                        System.out.print(s.getSites());
			Chain chain = new ChainImpl();
			chain.setChainID("H");
			for (Site site : s.getSites()) {
				//System.out.println("Site: " + site.getSiteID());
				for (Group group : site.getGroups()) {
					//manually add the chain as this is added later once the ATOM recrds are parsed usually
					group.setChain(chain);
					//					System.out.println("    PDBName: " + group.getPDBName());
					//					System.out.println("    PDBCode: " + group.getPDBCode());
					//					System.out.println("    Type: " + group.getType());
					//					System.out.println("    Parent: " + group.getChainId());
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


		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

	public void testMultiLineJRNL(){
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
		try {
			s = parser.parsePDBFile(br);
		} catch (IOException ex) {
			Logger.getLogger(PDBFileParserTest.class.getName()).log(Level.SEVERE, null, ex);
		}
		// String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("293", journalArticle.getVolume());
		assertEquals("53", journalArticle.getStartPage());
		assertEquals(1981, journalArticle.getPublicationDate());
		assertEquals("PHILOS.TRANS.R.SOC.LONDON, SER.B", journalArticle.getJournalName());
	}

	public void testIncorrectDateFormatMultiLineJRNL(){
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
		try {
			s = parser.parsePDBFile(br);
		} catch (IOException ex) {
			Logger.getLogger(PDBFileParserTest.class.getName()).log(Level.SEVERE, null, ex);
		}
		//    String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("293", journalArticle.getVolume());
		assertEquals("53", journalArticle.getStartPage());
		assertEquals(0, journalArticle.getPublicationDate());
		assertEquals("PHILOS.TRANS.R.SOC.LONDON, SER.B", journalArticle.getJournalName());
	}

	public void testInvalidFormatREFsectionJRNL(){
		//            System.out.println("Testing JRNL record parsing from 3pfk");
		String jrnlString =
				"JRNL        AUTH   P.R.EVANS,G.W.FARRANTS,P.J.HUDSON                            " + newline +
				//            "JRNL        TITL   PHOSPHOFRUCTOKINASE: STRUCTURE AND CONTROL.                  " + newline +
				"JRNL        REF    INTERESTING TIMES                                            " + newline +
				"JRNL        REFN                   ISSN 0080-4622                               " + newline +
				"JRNL        PMID   6115424                                                      ";


		BufferedReader br = new BufferedReader(new StringReader(jrnlString));
		Structure s = null;
		try {
			s = parser.parsePDBFile(br);
		} catch (IOException ex) {
			Logger.getLogger(PDBFileParserTest.class.getName()).log(Level.SEVERE, null, ex);
		}
		// String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("", journalArticle.getVolume());
		assertEquals("", journalArticle.getStartPage());
		assertEquals(0, journalArticle.getPublicationDate());
		assertEquals("", journalArticle.getJournalName());
	}


	public void testSecondMultiLineJRNL(){
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
		try {
			s = parser.parsePDBFile(br);
		} catch (IOException ex) {
			Logger.getLogger(PDBFileParserTest.class.getName()).log(Level.SEVERE, null, ex);
		}
		// String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		assertEquals("", journalArticle.getVolume());
		assertEquals("1", journalArticle.getStartPage());
		assertEquals(1991, journalArticle.getPublicationDate());
		assertEquals("GLYCOGEN PHOSPHORYLASE B: DESCRIPTION OF THE PROTEIN STRUCTURE", journalArticle.getJournalName());
	}

	public void testSingleLineJRNL(){
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
		try {
			s = parser.parsePDBFile(br);
		} catch (IOException ex) {
			Logger.getLogger(PDBFileParserTest.class.getName()).log(Level.SEVERE, null, ex);
		}
		//   String jrnl = s.getJournalArticle().toString();
		//            System.out.println(jrnl);
		JournalArticle journalArticle = s.getJournalArticle();
		//            System.out.println(journalArticle.getRef());
		assertEquals("280", journalArticle.getVolume());
		assertEquals("23000", journalArticle.getStartPage());
		assertEquals(2005, journalArticle.getPublicationDate());
		assertEquals("J.BIOL.CHEM.", journalArticle.getJournalName());
	}

	public void testToBePublishedJRNL(){
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
		try {
			s = parser.parsePDBFile(br);
		} catch (IOException ex) {
			Logger.getLogger(PDBFileParserTest.class.getName()).log(Level.SEVERE, null, ex);
		}

		JournalArticle journalArticle = s.getJournalArticle();
		assertNull(journalArticle.getVolume());
		assertNull(journalArticle.getStartPage());
		assertEquals(0,journalArticle.getPublicationDate());
		assertEquals("TO BE PUBLISHED", journalArticle.getJournalName());
	}

	public void test4hhbAcceptedAtomNames(){

		FileParsingParameters params = new FileParsingParameters();

		String[] acceptedAtoms = {StructureTools.caAtomName, " CB "};
		params.setAcceptedAtomNames(acceptedAtoms);
		AtomCache cache = new AtomCache();

		cache.setFileParsingParams(params);
		try {
			Structure s = cache.getStructure("4HHB");

			//System.out.println(s.toPDB());
			Atom[] ca = StructureTools.getAtomCAArray(s);
			Atom[] cb = StructureTools.getAtomArray(s, new String[]{" CB ",});

			// gly does not have cb...
			assertTrue(ca.length > cb.length);
			assertTrue(cb.length > 500);

			Atom[] allwithoutGLY = StructureTools.getAtomArray(s, acceptedAtoms);

			assertEquals(allwithoutGLY.length, (ca.length + cb.length - (ca.length-cb.length)));

		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
}
