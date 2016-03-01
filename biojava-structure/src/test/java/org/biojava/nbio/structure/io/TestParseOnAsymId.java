package org.biojava.nbio.structure.io;

import java.io.IOException;
import java.util.List;
import static org.junit.Assert.*;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestParseOnAsymId {



	@Test
	public void test4cup() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setUseInternalChainId(true);
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		// Set the test lists
		String[] asymChainList = {"A","B","C","D","E","F"};
		String[] authChainList = {"A","A","A","A","A","A"};
		String[] asymChainListTest = new String[6];
		String[] authChainListTest = new String[6];
		// Get the structure
		Structure bioJavaStruct = StructureIO.getStructure("4cup");
		List<Chain> chainList = bioJavaStruct.getChains();
		assertEquals(6,chainList.size());
		for(int i=0; i<chainList.size();i++){
			Chain c = chainList.get(i);
			authChainListTest[i] = c.getInternalChainID();
			asymChainListTest[i] = c.getChainID();
		}
		// Now check both lists are the same
		assertArrayEquals(authChainListTest, authChainList);
		assertArrayEquals(asymChainListTest, asymChainList);

		
		params.setUseInternalChainId(false);
		Structure bioJavaStructDiff = StructureIO.getStructure("4cup");
		List<Chain> chainListDiff = bioJavaStructDiff.getChains();
		assertEquals(1,chainListDiff.size());
		
		String[] authChainListTestDiff = new String[1];
		String[] asymChainListTestDiff = new String[1];
		for(int i=0; i<chainListDiff.size();i++){
			Chain c = chainListDiff.get(i);
			authChainListTestDiff[i] = c.getInternalChainID();
			asymChainListTestDiff[i] = c.getChainID();
		}
		assertArrayEquals(authChainListTestDiff, asymChainListTestDiff);
		assertArrayEquals(new String[] {"A"}, asymChainListTestDiff);
		
	}

}