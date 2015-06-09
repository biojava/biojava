package org.biojava.nbio.structure.io;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.junit.Test;

public class TestMMCIFWriting {

	@Test
	public void test1SMT() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(true);
		
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		
		Structure originalStruct = StructureIO.getStructure("1SMT");
				
		File outputFile = File.createTempFile("biojava_testing_", ".cif");
		
		
		FileWriter fw = new FileWriter(outputFile);
		fw.write(originalStruct.toMMCIF());
		fw.close();
		
		
		MMcifParser parser = new SimpleMMcifParser();

		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();

		FileParsingParameters fileParsingParams = new FileParsingParameters();
		fileParsingParams.setAlignSeqRes(true);

		consumer.setFileParsingParameters(fileParsingParams);

		parser.addMMcifConsumer(consumer);

		//parser.parse(new BufferedReader(new FileReader(new File("/home/duarte_j/test.cif")))); 
		parser.parse(new BufferedReader(new FileReader(outputFile)));
		
		Structure readStruct = consumer.getStructure();
		
		assertNotNull(readStruct); 
		
		assertEquals(originalStruct.getChains().size(), readStruct.getChains().size());
		
		for (int i=0;i<originalStruct.getChains().size();i++) {
			assertEquals(originalStruct.getChains().get(i).getAtomGroups().size(),
							readStruct.getChains().get(i).getAtomGroups().size());
			
			Chain origChain = originalStruct.getChains().get(i);
			Chain readChain = readStruct.getChains().get(i);
			
			assertEquals(origChain.getAtomGroups().size(), readChain.getAtomGroups().size());
			//assertEquals(origChain.getSeqResGroups().size(), readChain.getSeqResGroups().size());
		}
		
	}

}
