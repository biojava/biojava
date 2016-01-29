package org.biojava.nbio.structure;

import static org.junit.Assert.*;
import static org.junit.Assume.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;

import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestURLIdentifier {
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
		assertEquals("2pos",id.toCanonical().getPdbId());
//		assertEquals("2pos",full.getName()); // What should this get set to with identifiers?

		url = new URL("file://" + base + "?residues=A:1-5");
		id = new URLIdentifier(url);
		assertEquals("wrong canonical for residues=A:1-5","2pos.A_1-5",id.toCanonical().toString());
		
		full = id.loadStructure(cache);
		assertNotNull(full);
		reduced = id.reduce(full);
		assertEquals("wrong length for residues=A:1-5", 5, StructureTools.getRepresentativeAtomArray(reduced).length);

		url = new URL("file://" + base + "?chainId=A");
		id = new URLIdentifier(url);
		assertEquals("wrong canonical for chainId=A","2pos.A",id.toCanonical().toString());

		full = id.loadStructure(cache);
		assertNotNull(full);
		reduced = id.reduce(full);
		assertEquals("wrong length for chainId=A", 94, StructureTools.getRepresentativeAtomArray(reduced).length);

	}
	
	@Test
	public void testHomeDir() throws URISyntaxException, IOException, StructureException {
		// Sanity check. Maybe some weird systems like Travis don't have homes?
		assumeNotNull(System.getProperty("user.home"));
		assumeTrue(new File(System.getProperty("user.home")).exists());
		
		// Test that ~ gets expanded to the home directory properly
		StructureName name = new StructureName("~/");
		assertTrue(name.isFile());
		
		name = new StructureName("~/1abc_this_file_should_not_exits_98nb2.pdb");
		assertTrue(name.isFile());
		
		// Should guess the PDB code, but of course nothing to load
		SubstructureIdentifier canon = name.toCanonical();
		assertEquals("1abc", canon.getPdbId());
		
		try {
			name.loadStructure(new AtomCache());
			fail("Our oddly-named test file exists!?!");
		} catch(FileNotFoundException e) {}
	}
}
