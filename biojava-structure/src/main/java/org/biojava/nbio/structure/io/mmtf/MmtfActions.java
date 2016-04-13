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
	 * Utility function to get a Biojava structure from a byte array.
	 * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
	 * @return a Biojava structure object relating to the input byte array.
	 * @throws IOException 
	 */
	public static Structure readBiojavaStruct(String filePath) throws IOException {
		// Get the reader - this is the bit that people need to implement.
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		// Set up the inflator
		DecoderToReader getToInflator = new DecoderToReader();
		// Do the inflation
		getToInflator.read(new DefaultDecoder(ReaderUtils.getDataFromFile(filePath)), mmtfStructureReader);
		// Get the structue
		return mmtfStructureReader.getStructure();
	}
	
	/**
	 * Utility function to write a Biojava structure object to a path.
	 * @param structure the structure to write
	 * @param path the full path of the file to write
	 * @throws IOException
	 */
	public static void writeBiojavaStruct(Structure structure, String path) throws IOException {
		// Set up this writer
		WriterToEncoder inflatorToGet = new WriterToEncoder();
		// Get the writer - this is what people implement
		MmtfStructureWriter mmtfStructureWriter = new MmtfStructureWriter(structure);
		// Now pass to the get API
		mmtfStructureWriter.write(inflatorToGet);
		// Now write this dat to file
		WriterUtils.writeDataToFile(inflatorToGet, path);
	}

}