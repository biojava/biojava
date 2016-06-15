package org.biojava.nbio.structure.io.mmtf;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;
import org.rcsb.mmtf.decoder.StructureDataToAdapter;
import org.rcsb.mmtf.encoder.AdapterToStructureData;

/**
 * Tests to see if roundtripping of MMTF can be done.
 * @author Anthony Bradley
 *
 */
public class TestMmtfRoundTrip {

	/**
	 * Test that we can round trip a simple structure.
	 * @throws IOException an error reading the file
	 * @throws StructureException an error parsing the structure
	 */
	@Test
	public void testRoundTrip() throws IOException, StructureException {
		Structure structure = StructureIO.getStructure("4CUP");
		AdapterToStructureData writerToEncoder = new AdapterToStructureData();
		new MmtfStructureWriter(structure, writerToEncoder);
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		new StructureDataToAdapter(writerToEncoder, mmtfStructureReader);
		mmtfStructureReader.getStructure();
	}
}
