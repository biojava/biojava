package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmtf.MmtfStructureWriter;
import org.rcsb.mmtf.encoder.GetToBean;
import org.rcsb.mmtf.encoder.InflatorToGet;
import org.rcsb.mmtf.serializers.MessagePackSerializer;

public class DemoMmtfWriter {

	/**
	 * Utility function to get a byte array from a Biojava structure
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @param parsingParams
	 * @return
	 * @throws IOException 
	 */
	public static byte[] getByteArray(Structure structure) throws IOException {
		// Set up the transform from the inflator to the get api
		InflatorToGet inflatorToGet = new InflatorToGet();
		// Get the writer - this is what people implement
		MmtfStructureWriter mmtfStructureWriter = new MmtfStructureWriter(structure);
		// Now deflate
		mmtfStructureWriter.write(inflatorToGet);
		// Get to bean
		GetToBean getToBean = new GetToBean(inflatorToGet);
		MessagePackSerializer messagePackSerializer = new MessagePackSerializer();
		return messagePackSerializer.serialize(getToBean.getMmtfBean());
	}
}
