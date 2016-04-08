package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmtf.MmtfStructureReader;
import org.rcsb.mmtf.dataholders.MmtfBean;
import org.rcsb.mmtf.decoder.BeanToGet;
import org.rcsb.mmtf.decoder.GetToInflator;
import org.rcsb.mmtf.deserializers.MessagePackDeserializer;

public class DemoMmtfReader {

	/**
	 * Utility function to get a Biojava structure from a byte array.
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @param parsingParams
	 * @return
	 * @throws IOException 
	 */
	public static Structure getBiojavaStruct(byte[] inputByteArray) throws IOException {
		// Get the reader - this is the bit that people need to implement.
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		// Set up the deserializer
		MessagePackDeserializer messagePackDeserializer = new MessagePackDeserializer();
		// Get the data
		MmtfBean mmtfBean = messagePackDeserializer.deserialize(inputByteArray);
		// Set up the data API
		BeanToGet beanToGet = new BeanToGet(mmtfBean);
		// Set up the inflator
		GetToInflator getToInflator = new GetToInflator();
		// Do the inflation
		getToInflator.read(beanToGet, mmtfStructureReader);
		// Get the structue
		return mmtfStructureReader.getStructure();
	}
}
