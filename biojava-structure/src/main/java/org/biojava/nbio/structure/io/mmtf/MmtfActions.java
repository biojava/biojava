package org.biojava.nbio.structure.io.mmtf;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.rcsb.mmtf.api.MmtfDecodedDataInterface;
import org.rcsb.mmtf.dataholders.MmtfBean;
import org.rcsb.mmtf.decoder.BeanToDataApi;
import org.rcsb.mmtf.decoder.DataApiToReader;
import org.rcsb.mmtf.deserializers.MessagePackDeserializer;
import org.rcsb.mmtf.encoder.DataApiToBean;
import org.rcsb.mmtf.encoder.WriterToDataApi;
import org.rcsb.mmtf.serializers.MessagePackSerializer;
import org.rcsb.mmtf.utils.DownloadUtils;

public class MmtfActions {

	/**
	 * Utility function to get a Biojava structure from a PDB code.
	 * @param pdbCode the pdb code of the structure you desire
	 * @return a Biojava structure object relating to the input PDB code
	 * @throws IOException 
	 */
	public static Structure getBiojavaStruct(String pdbCode) throws IOException {
		return getBiojavaStruct(DownloadUtils.getDataFromUrl(pdbCode));
	}
	
	/**
	 * Utility function to get a Biojava structure from a byte array.
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @return a Biojava structure object relating to the input byte array.
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
		BeanToDataApi beanToGet = new BeanToDataApi(mmtfBean);
		// Set up the inflator
		DataApiToReader getToInflator = new DataApiToReader();
		// Do the inflation
		getToInflator.read(beanToGet, mmtfStructureReader);
		// Get the structue
		return mmtfStructureReader.getStructure();
	}

	/**
	 * Get the byte array (messagepack encoded) of the PDB code you desire.
	 * @param pdbId the input PDB id for the structure you want
	 * @return the byte array of the structure compressed and message pack encoded
	 * @throws IOException
	 * @throws StructureException
	 */
	public static byte[] getByteArray(String pdbId) throws IOException, StructureException {
		MmtfUtils.setUpBioJava();
		return getByteArray(StructureIO.getStructure(pdbId));
	}

	/**
	 * Utility function to get a byte array from a Biojava structure.
	 * @param structure the input Biojava structure object
	 * @return the byte array of the structure compressed and message pack encoded
	 * @throws IOException 
	 */
	public static byte[] getByteArray(Structure structure) throws IOException {
		MessagePackSerializer messagePackSerializer = new MessagePackSerializer();
		return messagePackSerializer.serialize(getBean(structure));
	}

	/**
	 * Utility function to get an mmtf bean from a Biojava structure.
	 * @param structure the input Biojava structure object
	 * @return the raw (compressed) data as an MmtfBean object
	 * @throws IOException 
	 */
	public static MmtfBean getBean(Structure structure) throws IOException {
		// Get to bean
		DataApiToBean getToBean = new DataApiToBean(getApi(structure));
		return getToBean.getMmtfBean();
	}
	
	/**
	 * Utility function to get an mmtf bean from a PDB id.
	 * @param pdbId the input PDB id for the structure you want
	 * @return the byte array of the structure compressed and message pack encoded
	 * @throws IOException 
	 * @throws StructureException 
	 */
	public static MmtfBean getBean(String pdbId) throws IOException, StructureException {
		// Get to bean
		DataApiToBean getToBean = new DataApiToBean(getApi(pdbId));
		return getToBean.getMmtfBean();
	}
	
	/**
	 * Utility function to get an API to the data  from a Biojava structure
	 * @param structure the input Biojava structure object
	 * @return the API to the data in the form of an MmtfDecodedDataInterface
	 * @throws IOException 
	 */
	public static MmtfDecodedDataInterface getApi(Structure structure) throws IOException {
		// Set up the transform from the inflator to the get api
		WriterToDataApi inflatorToGet = new WriterToDataApi();
		// Get the writer - this is what people implement
		MmtfStructureWriter mmtfStructureWriter = new MmtfStructureWriter(structure);
		// Now deflate
		mmtfStructureWriter.write(inflatorToGet);
		// Return the API
		return inflatorToGet;
	}

	/**
	 * Utility function to get an API to the data from a PDB code.
	 * @param pdbId the input PDB id for the structure you want
	 * @return the API to the data in the form of an MmtfDecodedDataInterface
	 * @throws IOException 
	 * @throws StructureException 
	 */
	public static MmtfDecodedDataInterface getApi(String pdbId) throws IOException, StructureException {
		// Set up the transform from the inflator to the get api
		WriterToDataApi inflatorToGet = new WriterToDataApi();
		// Get the writer - this is what people implement
		MmtfStructureWriter mmtfStructureWriter = new MmtfStructureWriter(
				StructureIO.getStructure(pdbId));
		// Now deflate
		mmtfStructureWriter.write(inflatorToGet);
		// Get to bean
		return inflatorToGet;
	}

	/**
	 * Round trip a given structure. Mainly for testing purposes.
	 * @param structure the input Biojava structure object
	 * @return the round tripped structure (conversion to messge pack mmtf and back).
	 * @throws IOException
	 */
	public static Structure roundTrip(Structure structure) throws IOException {
		return getBiojavaStruct(getByteArray(structure));
	}
}
