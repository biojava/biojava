package org.biojava.bio.structure;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileParser;
import org.junit.Test;

public class TestCompoundHeuristics {
	
	private static final String PATH_TO_TEST_FILES = "/org/biojava/bio/structure/io/";
	
	//private static final String[] FILES_TO_TEST = 
	//	{"1b8g_raw.pdb.gz", "2m7y_raw.pdb.gz", "3c5f_raw.pdb.gz", "4b19_raw.pdb.gz"};


	
	//public void testCompoundHeuristics() throws IOException {
	//	for (String fileName:FILES_TO_TEST) {
	//		testFile(fileName);
	//	}
	//}
	
	@Test
	public void test1b8gRaw() throws IOException, StructureException { 
		
		Structure s = getStructure("1b8g_raw.pdb.gz", true);
		
		assertEquals(1,s.getCompounds().size());
		
		Chain chainA = s.getChainByPDB("A");
		
		assertEquals(2, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());

		checkCompoundsNumbered(s.getCompounds());
		
		
		s = getStructure("1b8g_raw.pdb.gz", false);
		
		assertEquals(1,s.getCompounds().size());
		
		chainA = s.getChainByPDB("A");
		
		assertEquals(2, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());

		checkCompoundsNumbered(s.getCompounds());
	}

	
	@Test
	public void test2m7yRaw() throws IOException, StructureException { 
		
		Structure s = getStructure("2m7y_raw.pdb.gz", true);
		
		assertEquals(1,s.getCompounds().size());
		
		Chain chainA = s.getChainByPDB("A");
		
		assertEquals(1, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());

		checkCompoundsNumbered(s.getCompounds());
		
		
		s = getStructure("2m7y_raw.pdb.gz", false);
		
		assertEquals(1,s.getCompounds().size());
		
		chainA = s.getChainByPDB("A");
		
		assertEquals(1, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());

		checkCompoundsNumbered(s.getCompounds());
	}
	
	@Test
	public void test3c5fRaw() throws IOException, StructureException { 
		
		Structure s = getStructure("3c5f_raw.pdb.gz", true);
		
		assertEquals(4,s.getCompounds().size());
		
		Chain chainA = s.getChainByPDB("A");
		
		assertEquals(2, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());

		checkCompoundsNumbered(s.getCompounds());
		
		
		s = getStructure("3c5f_raw.pdb.gz", false);
		
		assertEquals(4,s.getCompounds().size());
		
		chainA = s.getChainByPDB("A");
		
		assertEquals(2, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());

		checkCompoundsNumbered(s.getCompounds());
	}
	
	@Test
	public void test4b19Raw() throws IOException, StructureException { 
		
		Structure s = getStructure("4b19_raw.pdb.gz", true);
		
		assertEquals(1,s.getCompounds().size());
		
		Chain chainA = s.getChainByPDB("A");
		
		assertEquals(1, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());
		
		checkCompoundsNumbered(s.getCompounds());
		
		
		s = getStructure("4b19_raw.pdb.gz", false);
		
		assertEquals(1,s.getCompounds().size());
		
		chainA = s.getChainByPDB("A");
		
		assertEquals(1, chainA.getCompound().getChains().size());
		
		assertEquals(chainA,chainA.getCompound().getRepresentative());
		
		checkCompoundsNumbered(s.getCompounds());

	}

	private Structure getStructure(String fileName, boolean setAlignSeqRes) throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream(PATH_TO_TEST_FILES+fileName));
		
		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(setAlignSeqRes);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		System.out.println("Entities for file: "+fileName);
		for (Compound ent:s.getCompounds()) {
			System.out.print(ent.getRepresentative().getChainID()+":");
			for (Chain c:ent.getChains()) {
				System.out.print(" "+c.getChainID());
			}
			System.out.println();
		}
		
		return s;
	}

	private void checkCompoundsNumbered(List<Compound> compounds) {
		
		Collections.sort(compounds, new Comparator<Compound>() {

			@Override
			public int compare(Compound o1, Compound o2) {
				return new Integer(o1.getMolId()).compareTo(o2.getMolId());
			}
		});
		
		int id = 1;
		for (Compound compound:compounds) {
			assertEquals(id,compound.getMolId());
			id++;
		}

		
	}
}
