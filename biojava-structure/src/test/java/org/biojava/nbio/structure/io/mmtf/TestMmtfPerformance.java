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
package org.biojava.nbio.structure.io.mmtf;

import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.assertTrue;

/**
 * Test the performance of MMTF format in BioJava.
 *
 * @author Andreas Prlic
 * on 1/9/17.
 *
 */
public class TestMmtfPerformance {

	private static final Logger logger = LoggerFactory.getLogger(TestMmtfPerformance.class);

	private static final int NUMBER_OF_REPEATS = 10;

	private static String convertStreamToString(java.io.InputStream is) {
		try (
		java.util.Scanner s = new java.util.Scanner(is)){
			return s.useDelimiter("\\A").hasNext() ? s.next() : "";
		}
	}

	private byte[] getByteArrayFromInputStream(InputStream is) throws IOException {
		ByteArrayOutputStream buffer = new ByteArrayOutputStream();

		int nRead;
		byte[] data = new byte[16384];

		while ((nRead = is.read(data, 0, data.length)) != -1) {
			buffer.write(data, 0, nRead);
		}

		buffer.flush();

		return buffer.toByteArray();
	}

	@Test
	public void test3HBX() throws Exception{
		String pdbId = "3hbx";

		pdbId = pdbId.toUpperCase();

		URL url = new URL("https://files.rcsb.org/download/"+pdbId+".pdb.gz");

		String pdbFile = convertStreamToString(new GZIPInputStream(url.openStream()));

		long totalTimePDB = 0;
		long totalTimeMMTF = 0;

		byte[] pdbBytes = pdbFile.getBytes();

		PDBFileParser parser = new PDBFileParser();

		URL mmtfURL = new URL("https://mmtf.rcsb.org/v1.0/full/" + pdbId + ".mmtf.gz");

		byte[] mmtfdata = getByteArrayFromInputStream(new GZIPInputStream((mmtfURL.openStream())));

		// first make sure chemcomp cache is warmed up (chemcomp files are parsed). Like that we count the parsing time without the influence of chemcomp parsing
		MmtfActions.readFromInputStream(new ByteArrayInputStream(mmtfdata));
		parser.parsePDBFile(new ByteArrayInputStream(pdbBytes));

		for ( int i =0 ; i< NUMBER_OF_REPEATS ; i++) {

			long mmtfStart = System.nanoTime();
			MmtfActions.readFromInputStream(new ByteArrayInputStream(mmtfdata));
			long mmtfEnd = System.nanoTime();

			long pdbStart = System.nanoTime();
			parser.parsePDBFile(new ByteArrayInputStream(pdbBytes));
			long pdbEnd = System.nanoTime();

			totalTimePDB += (pdbEnd - pdbStart);

			totalTimeMMTF += (mmtfEnd-mmtfStart);
		}

		long timePDB = (totalTimePDB/NUMBER_OF_REPEATS);
		long timeMMTF = (totalTimeMMTF/NUMBER_OF_REPEATS);

		logger.info("average time to parse mmtf: " + timeMMTF/(1000*1000) + " ms.");
		logger.info("average time to parse PDB : " + timePDB/(1000*1000) + " ms. ");

		assertTrue( "It should not be the case, but it is faster to parse a PDB file ("+timePDB+" ns.) than MMTF ("+( timeMMTF)+" ns.)!",( timePDB) > ( timeMMTF));

	}
}
