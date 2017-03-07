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
 */
package org.biojava.nbio.structure.io;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TestHeaderOnly {

	private static final Logger logger = LoggerFactory.getLogger(TestHeaderOnly.class);


	private final String pdbID = "1REP";

	/**
	 * All groups are expected to be empty.
	 *
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testHeaderOnly() throws StructureException, IOException {
		// Get either PDB or mmCIF with a headerOnly = true.

		// Test 1: with PDB
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);

		FileParsingParameters params = new FileParsingParameters();
		params.setHeaderOnly(true);
		// params.setAlignSeqRes(true);  // Now this is default.
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		Structure sPDB = StructureIO.getStructure(pdbID);

		Assert.assertEquals(false, doSeqResHaveAtoms(sPDB));
		
		// Test 2: with mmCIF
		cache.setUseMmCif(true);

		Structure sCIF = StructureIO.getStructure(pdbID);
		Assert.assertEquals(false, doSeqResHaveAtoms(sCIF));

		
	}

	/**
	 * Test that with alignSeqRes, expected Group(s) have Atoms, while others
	 * are present with correct sequence but empty.
	 *
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testAlignSeqres() throws StructureException, IOException {
		// Get either PDB or mmCIF with a headerOnly = false.

		// Test 1: with PDB
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);

		FileParsingParameters params = new FileParsingParameters();
		params.setHeaderOnly(false);
		// params.setAlignSeqRes(true);  // Now this is default.
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		Structure sPDB = StructureIO.getStructure(pdbID);
		Assert.assertEquals(true, doSeqResHaveAtoms(sPDB));
		check1REPChainC(sPDB); // Check particular residues to be aligned.

		// Test 2: with mmCIF
		cache.setUseMmCif(true);

		Structure sCIF = StructureIO.getStructure(pdbID);
		Assert.assertEquals(true, doSeqResHaveAtoms(sCIF));
		check1REPChainC(sCIF); // Check particular residues to be aligned.
	}

	// A better test follows that uses local files.
	// @Test
	public void testSpeed() {
		// Force using a file reader.
		MMCIFFileReader fr = new MMCIFFileReader();
		FileParsingParameters par = new FileParsingParameters();
		//par.setAlignSeqRes(true);
		// par.setHeaderOnly(true);
		par.setHeaderOnly(false);
		fr.setFileParsingParameters(par);
		fr.setFetchBehavior(FetchBehavior.FETCH_FILES);

		Structure s = null;
		long start = System.nanoTime();
		try {
			// Medium sized structure parsed in 0.549s (no header) vs .676s (header) ~ 20% faster
			s = fr.getStructureById("4WZ6");
			// A larger structure could be parsed ~ 4.991s (no header) vs 5.867s (header) ~ 16% faster
			// s = fr.getStructureById("4V60");
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		long stop = System.nanoTime();
		double diff = (stop - start) / 1000000000.0;
		logger.info(String.format("[%s] Elapsed time: %.3f s", s.getIdentifier(), diff));
	}

	// Test using local files.
	@Test
	public void testSpeed2() throws StructureException, IOException {
		// Test the file parsing speed when the files are already downloaded.

		InputStream cifStream = new GZIPInputStream(this.getClass().getResourceAsStream("/4hhb.cif.gz"));
		InputStream pdbStream = new GZIPInputStream(this.getClass().getResourceAsStream("/4hhb.pdb.gz"));

		assertNotNull(cifStream);

		FileParsingParameters params = new FileParsingParameters();
		params.setHeaderOnly(true);  // Flip this true/false to compare parsing speed.

		logger.info("Testing PDB parsing speed");
		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setFileParsingParameters(params);
		//pdbpars.setLoadChemCompInfo(true);
		long start = System.nanoTime();
		Structure s1 = pdbpars.parsePDBFile(pdbStream) ;
		long stop = System.nanoTime();
		double diff = (stop - start) / 1000000000.0;
		logger.info(String.format("[%s] Elapsed time: %.3f s", s1.getIdentifier(), diff));

		MMcifParser mmcifpars = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		consumer.setFileParsingParameters(params);
		mmcifpars.addMMcifConsumer(consumer);

		logger.info("Testing mmCIF parsing speed");
		start = System.nanoTime();
		mmcifpars.parse(cifStream) ;
		Structure s2 = consumer.getStructure();
		stop = System.nanoTime();
		diff = (stop - start) / 1000000000.0;
		logger.info(String.format("[%s] Elapsed time: %.3f s", s2.getIdentifier(), diff));

		/* Running from an SSD..
		 * PDB .165s (all atom) -> 0.009s (only header)  95% faster.
		 * mmCIF 0.323s (no header) -> 0.175s (only header) 45% faster.
		 */
	}

	/**
	 * Scan through SeqResGroups, returns true if any have Atoms.
	 * @param s
	 * @return
	 */
	public boolean doSeqResHaveAtoms(Structure s) {
		for (int i = 0; i < s.nrModels(); i++) {
			for (Chain c : s.getChains(i)) {
				for (Group g : c.getSeqResGroups()) {
					if (hasAtoms(g)) return true; // Found some Atoms in a Seqres group.
				}
			}
		}

		return false;
	}

	/**
	 * Does a group have any Atom(s)?
	 *
	 * @param g : a group
	 * @return true if has any Atom(s)
	 */
	public boolean hasAtoms(Group g) {
		if (g.getAtoms().size() > 0) return true;
		return false;
	}

	/**
	 * Check that the gapped residues have no atoms, but that ungapped residues
	 * have atoms.
	 *
	 * @param s: Structure to test.
	 */
	public void check1REPChainC(Structure s) throws StructureException {
		String sequence = "MAETAVINHKKRKNSPRIVQSNDLTEAAYSLSRDQKRMLYLFVDQIRK" +
				"SDGTLQEHDGICEIHVAKYAEIFGLTSAEASKDIRQALKSFAGKEVVFYRPEEDAGDE" +
				"KGYESFPWFIKPAHSPSRGLYSVHINPYLIPFFIGLQNRFTQFRLSETKEITNPYAMR" +
				"LYESLCQYRKPDGSGIVSLKIDWIIERYQLPQSYQRMPDFRRRFLQVCVNEINSRTPM" +
				"RLSYIEKKKGRQTTHIVFSFRDITSMTTG";

		boolean [] shouldMatch = new boolean[sequence.length()];
		for (int i = 0; i < sequence.length(); i++) shouldMatch[i] = true;

		// 1-14 is gap
		for (int i = 0; i < 14; i++) shouldMatch[i] = false;

		// 50-55 is gap
		for (int i = 49; i < 55; i++) shouldMatch[i] = false;

		// 98-109 is gap
		for (int i = 97; i < 109; i++) shouldMatch[i] = false;

		// 247-251 is gap
		for (int i = 246; i < 251; i++) shouldMatch[i] = false;

		Chain c = s.getPolyChainByPDB("C");

		List<Group> seqres = c.getSeqResGroups();

		// Check lengths
		Assert.assertEquals(sequence.length(), seqres.size());

		// Check sequences.
		Assert.assertEquals(sequence, c.getSeqResSequence());

		for (int i = 0; i < sequence.length(); i++) {
			Assert.assertEquals(shouldMatch[i], hasAtoms(seqres.get(i)));
		}
	}

	/**
	 *
	 * @param seqres : a list of Group(s)
	 * @return a String representing these Groups
	 */
	public String getSequenceString(List<Group> seqres) {
		StringBuilder sb = new StringBuilder();

		for (Group g : seqres) {
			ChemComp c = g.getChemComp();
			sb.append(c.getOne_letter_code());
		}

		return sb.toString();
	}
}
