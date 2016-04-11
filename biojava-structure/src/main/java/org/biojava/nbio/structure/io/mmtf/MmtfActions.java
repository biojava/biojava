package org.biojava.nbio.structure.io.mmtf;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.rcsb.mmtf.api.MmtfDecodedDataInterface;
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
	 * 
	 * @param pdbId
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 */
	public static byte[] getByteArray(String pdbId) throws IOException, StructureException {
		MmtfUtils.setUpBioJava();
		return getByteArray(StructureIO.getStructure(pdbId));
	}

	/**
	 * Utility function to get a byte array from a Biojava structure
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @param parsingParams
	 * @return
	 * @throws IOException 
	 */
	public static byte[] getByteArray(Structure structure) throws IOException {
		MessagePackSerializer messagePackSerializer = new MessagePackSerializer();
		return messagePackSerializer.serialize(getBean(structure));
	}

	/**
	 * Utility function to get an mmtf bean from a Biojava structure
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @param parsingParams
	 * @return
	 * @throws IOException 
	 */
	public static MmtfBean getBean(Structure structure) throws IOException {
		// Get to bean
		GetToBean getToBean = new GetToBean(getApi(structure));
		return getToBean.getMmtfBean();
	}
	
	/**
	 * Utility function to get an mmtf bean from a Biojava structure
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @param parsingParams
	 * @return
	 * @throws IOException 
	 * @throws StructureException 
	 */
	public static MmtfBean getBean(String pdbId) throws IOException, StructureException {
		// Get to bean
		GetToBean getToBean = new GetToBean(getApi(pdbId));
		return getToBean.getMmtfBean();
	}
	
	/**
	 * Utility function to get an mmtf bean from a Biojava structure
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @param parsingParams
	 * @return
	 * @throws IOException 
	 */
	public static MmtfDecodedDataInterface getApi(Structure structure) throws IOException {
		// Set up the transform from the inflator to the get api
		InflatorToGet inflatorToGet = new InflatorToGet();
		// Get the writer - this is what people implement
		MmtfStructureWriter mmtfStructureWriter = new MmtfStructureWriter(structure);
		// Now deflate
		mmtfStructureWriter.write(inflatorToGet);
		// Return the API
		return inflatorToGet;
	}

	/**
	 * Utility function to get an mmtf bean from a Biojava structure
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @param parsingParams
	 * @return
	 * @throws IOException 
	 * @throws StructureException 
	 */
	public static MmtfDecodedDataInterface getApi(String pdbId) throws IOException, StructureException {
		// Set up the transform from the inflator to the get api
		InflatorToGet inflatorToGet = new InflatorToGet();
		// Get the writer - this is what people implement
		MmtfStructureWriter mmtfStructureWriter = new MmtfStructureWriter(
				StructureIO.getStructure(pdbId));
		// Now deflate
		mmtfStructureWriter.write(inflatorToGet);
		// Get to bean
		return inflatorToGet;
	}

	/**
	 * Round trip a given structure. Mainly for testing purposes
	 * @param structure the input structure
	 * @return the round tripped structure (conversion to messge pack mmtf and back).
	 * @throws IOException
	 */
	public static Structure roundTrip(Structure structure) throws IOException {
		return getBiojavaStruct(getByteArray(structure));
	}
}
