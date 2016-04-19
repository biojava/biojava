package org.biojava.nbio.structure.io.mmtf;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.rcsb.mmtf.decoder.DefaultDecoder;
import org.rcsb.mmtf.decoder.DecoderToReader;
import org.rcsb.mmtf.decoder.ReaderUtils;
import org.rcsb.mmtf.encoder.WriterToEncoder;
import org.rcsb.mmtf.encoder.WriterUtils;

/**
 * A class of functions for reading and writing Biojava structures using MMTF
 * @author Anthony Bradley
 *
 */
public class MmtfActions {
	
	/**
	 * Utility function to get a Biojava structure from a file.
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @return a Biojava structure object relating to the input byte array.
	 * @throws IOException 
	 */
	public static Structure readFromFile(String filePath) throws IOException {
		// Get the reader - this is the bit that people need to implement.
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		// Do the inflation
		new DecoderToReader(new DefaultDecoder(ReaderUtils.getDataFromFile(filePath)), mmtfStructureReader);
		// Get the structue
		return mmtfStructureReader.getStructure();
	}
	
	/**
	 * Utility function to write a Biojava structure object to a path.
	 * @param structure the structure to write
	 * @param path the full path of the file to write
	 * @throws IOException
	 */
	public static void writeToFile(Structure structure, String path) throws IOException {
		// Set up this writer
		WriterToEncoder writerToEncoder = new WriterToEncoder();
		// Get the writer - this is what people implement
		new MmtfStructureWriter(structure, writerToEncoder);
		// Now write this data to file
		WriterUtils.writeDataToFile(writerToEncoder, path);
	}

	
	/**
	 * Utility function to get a Biojava structure from the REST service.
	 * @param pdbId the PDB code of the required structure
	 * @return a Biojava structure object relating to the input byte array
	 * @throws IOException 
	 */
	public static Structure readFromWeb(String pdbId) throws IOException {
		// Get the reader - this is the bit that people need to implement.
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		// Do the inflation
		new DecoderToReader(new DefaultDecoder(ReaderUtils.getDataFromUrl(pdbId)), mmtfStructureReader);
		// Get the structue
		return mmtfStructureReader.getStructure();
	}
}