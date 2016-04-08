package org.biojava.nbio.structure.mmtf;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmtf.MmtfStructureReader;
import org.junit.Test;
import org.rcsb.mmtf.examples.HandleIO;

/**
 * Tests for the MMTF parser in Biojava
 * @author Anthony Bradley
 *
 */
public class TestMmtfParser {

	/**
	 * Can we parse an MMTF file and get the right number of chains.
	 * @throws IOException 
	 */
	@Test
	public void basicTest() throws IOException {
	    HandleIO gbjs = new HandleIO();
	    byte[] inputByteArr = gbjs.getFromUrl("1qmz");
	    // Now get the biojava structure.
	    Structure biojavaStruct = MmtfStructureReader.getBiojavaStruct(inputByteArr);
	    assertEquals(biojavaStruct.getChains().size(), 6);
	    assertEquals(biojavaStruct.getPDBCode().toLowerCase(), "1qmz");
	}
}
