package org.biojava.nbio.structure;

import static org.junit.Assert.*;

import java.io.IOException;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.Arrays;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TestURLIdentifier {
	private static final Logger logger = LoggerFactory.getLogger(TestURLIdentifier.class);
	@Test
	public void testFileIdentifiers() throws StructureException, IOException {
		AtomCache cache = new AtomCache();
		
		URL url;
		Structure full, reduced;
		URLIdentifier id;
		
		url = getClass().getResource("/2pos.pdb");
		id = new URLIdentifier(url);
		
		full = id.loadStructure(cache);
		assertNotNull("PDB file didn't work",full);
		
		reduced = id.reduce(full);
		assertTrue(Arrays.deepEquals(StructureTools.getAllAtomArray(full),
				StructureTools.getAllAtomArray(reduced)));
		
		url = getClass().getResource("/4hhb.cif.gz");
		id = new URLIdentifier(url);
		
		full = id.loadStructure(cache);
		assertNotNull("CIF file didn't work",full);
		
		reduced = id.reduce(full);
		assertTrue(Arrays.deepEquals(StructureTools.getAllAtomArray(full),
				StructureTools.getAllAtomArray(reduced)));
		
	}
	
	@Test
	public void testURLParameters() throws StructureException, IOException {
		AtomCache cache = new AtomCache();
		
		URL url;
		Structure full, reduced;
		URLIdentifier id;
		
		String base = getClass().getResource("/2pos.pdb").getPath();
		
		url = new URL("file://" + base + "?format=PDB");
		id = new URLIdentifier(url);

		full = id.loadStructure(cache);
		assertNotNull(full);
		assertEquals("2POS",id.toCanonical().getPdbId());
//		assertEquals("2POS",full.getName()); // What should this get set to with identifiers?

		url = new URL("file://" + base + "?residues=A:1-5");
		id = new URLIdentifier(url);
		assertEquals("wrong canonical for residues=A:1-5","2POS.A_1-5",id.toCanonical().toString());
		
		full = id.loadStructure(cache);
		assertNotNull(full);
		reduced = id.reduce(full);
		assertEquals("wrong length for residues=A:1-5", 5, StructureTools.getRepresentativeAtomArray(reduced).length);

		url = new URL("file://" + base + "?chainId=A");
		id = new URLIdentifier(url);
		assertEquals("wrong canonical for chainId=A","2POS.A",id.toCanonical().toString());

		full = id.loadStructure(cache);
		assertNotNull(full);
		reduced = id.reduce(full);
		assertEquals("wrong length for chainId=A", 94, StructureTools.getRepresentativeAtomArray(reduced).length);

		try {
			url = new URL("http://www.rcsb.org/pdb/files/1B8G.pdb.gz");
			id = new URLIdentifier(url);

			full = id.loadStructure(cache);
			assertNotNull(full);
			assertEquals("1B8G",id.toCanonical().getPdbId());
		} catch(UnknownHostException e) {
			logger.error("Unable to connect to rcsb.org");
			// still pass
		}
	}
}
