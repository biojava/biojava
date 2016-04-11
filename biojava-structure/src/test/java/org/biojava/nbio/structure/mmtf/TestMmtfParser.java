package org.biojava.nbio.structure.mmtf;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmtf.MmtfStructureReader;
import org.junit.Test;
import org.rcsb.mmtf.dataholders.MmtfBean;
import org.rcsb.mmtf.decoder.BeanToDataApi;
import org.rcsb.mmtf.decoder.DataApiToReader;
import org.rcsb.mmtf.deserializers.MessagePackDeserializer;
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
	    byte[] inputbyteArr = gbjs.getFromUrl("1qmz");
		// Get the reader - this is the bit that people need to implement.
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		// Set up the deserializer
		MessagePackDeserializer messagePackDeserializer = new MessagePackDeserializer();
		// Get the data
		MmtfBean mmtfBean = messagePackDeserializer.deserialize(inputbyteArr);
		// Set up the data API
		BeanToDataApi beanToGet = new BeanToDataApi(mmtfBean);
		// Set up the inflator
		DataApiToReader getToInflator = new DataApiToReader();
		// Do the inflation
		getToInflator.read(beanToGet, mmtfStructureReader);
		// Get the structue
		Structure biojavaStruct = mmtfStructureReader.getStructure();
	    assertEquals(biojavaStruct.getChains().size(), 6);
	    assertEquals(biojavaStruct.getPDBCode().toLowerCase(), "1qmz");
	}
}
