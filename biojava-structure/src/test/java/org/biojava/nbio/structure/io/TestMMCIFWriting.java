package org.biojava.nbio.structure.io;

import static org.junit.Assert.*;

import java.io.FileWriter;
import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
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
		
		Structure sCif = StructureIO.getStructure("1SMT");
		
		assertNotNull(sCif);
		

		FileConvert fc = new FileConvert(sCif);
		
		FileWriter fw = new FileWriter("/home/duarte_j/test.cif");
		//System.out.println(fc.toMmCif());
		fw.write(fc.toMmCif());
		fw.close();
	}

}
