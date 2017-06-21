package org.biojava.nbio.structure.io;

import static org.junit.Assert.fail;

import java.io.IOException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestWriteLargeCoordinatePDB {
	
	// This test checks that 'grouping' characters such as commas are not
	// incorrectly introduced into formatted PDB coordinate fields.
	// See FileConvert.d3 formatter.
	@Test
	public void TestWrite5D9Q() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);

		FileParsingParameters params = new FileParsingParameters();
		params.setHeaderOnly(false);
		cache.setFileParsingParams(params);

		StructureIO.setAtomCache(cache);

		// Example structure with large coordinates in PDB file.
		Structure sPDB = StructureIO.getStructure("5D9Q");
		
		// If 48 column for a ATOM/HETATM has a comma, fail.
		for (Group g : sPDB.getChain("K").getAtomGroups()) {
			for (Atom a : g.getAtoms()) {
				if (a.toPDB().contains(",")) {
					fail("Comma present in ATOM/HETATM record.");
				}
			}
		}
		
		//try (PrintWriter p = new PrintWriter(new FileWriter(new File("/tmp/test.pdb")))) {
		//	p.print(sPDB.toPDB());
		//}
	}
	
}
