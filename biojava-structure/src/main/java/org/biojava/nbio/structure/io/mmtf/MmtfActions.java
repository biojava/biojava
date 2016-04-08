package org.biojava.nbio.structure.io.mmtf;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.rcsb.mmtf.dataholders.MmtfBean;
import org.rcsb.mmtf.decoder.BeanToGet;
import org.rcsb.mmtf.decoder.GetToInflator;
import org.rcsb.mmtf.deserializers.MessagePackDeserializer;
import org.rcsb.mmtf.encoder.GetToBean;
import org.rcsb.mmtf.encoder.InflatorToGet;
import org.rcsb.mmtf.serializers.MessagePackSerializer;

public class MmtfActions {

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
	
	/**
	 * 
	 * @param structure
	 * @return
	 * @throws IOException
	 */
	public static Structure roundTrip(Structure structure) throws IOException {
		return getBiojavaStruct(getByteArray(structure));
	}
}
