package org.biojava.nbio.structure.test.io;

import static org.junit.Assume.assumeNotNull;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.cif.CifStructureConverter;
import org.junit.jupiter.api.Test;

/**
 * Test that various identifiers in Structure are parsed consistently across
 * file types.
 * 
 * These tests are not intended to be prescriptive, but merely to document the
 * current state and force changes to be committed to git.
 */
public class TestIdentifiersInStructure {
	@Test
	public void testCifIdentifiers() throws IOException, StructureException {
		String noid = "/5pti_noid.cif";
		String yesid = "/5pti.cif";

		// No ID, from file
		URL url = this.getClass().getResource(noid);
		String file = url.getPath();
		Structure s = CifStructureConverter.fromURL(url);
		assertNull(s.getPdbId());
		assertEquals("", s.getIdentifier());
		assertEquals("", s.getName());
		assertNull(s.getStructureIdentifier());

		// No ID, from StructureIO
		s = StructureIO.getStructure(file);
		assertNull(s.getPdbId());
		assertEquals(file, s.getIdentifier());
		assertEquals("5PTI", s.getName());
		StructureName sid = (StructureName) s.getStructureIdentifier();
		assertEquals(file, sid.getIdentifier());

		// No id, from stream
		InputStream stream = this.getClass().getResourceAsStream(noid);
		s = CifStructureConverter.fromInputStream(stream);
		assertNull(s.getPdbId());
		assertEquals("", s.getIdentifier());
		assertEquals("", s.getName());
		assertNull(s.getStructureIdentifier());

		// With ID, from file
		url = this.getClass().getResource(yesid);
		file = url.getPath();
		s = CifStructureConverter.fromURL(url);
		assertEquals("5PTI", s.getPdbId().getId());
		assertEquals("", s.getIdentifier());
		assertEquals("", s.getName());
		assertNull(s.getStructureIdentifier());

		// With ID, from StructureIO
		s = StructureIO.getStructure(file);
		assertEquals("5PTI", s.getPdbId().getId());
		assertEquals(file, s.getIdentifier());
		assertEquals("5PTI", s.getName());
		sid = (StructureName) s.getStructureIdentifier();
		assertEquals(file, sid.getIdentifier());

		// With id, from stream
		stream = this.getClass().getResourceAsStream(yesid);
		s = CifStructureConverter.fromInputStream(stream);
		assertEquals("5PTI", s.getPdbId().getId());
		assertEquals("", s.getIdentifier());
		assertEquals("", s.getName());
		assertNull(s.getStructureIdentifier());
	}

	@Test
	public void testPDBIdentifiers() throws IOException, StructureException {
		String yesid = "/5pti.pdb";
		String noid = "/noid.pdb";

		assumeNotNull(this.getClass().getResource(yesid));
		assumeNotNull(this.getClass().getResource(noid));

		PDBFileParser parser = new PDBFileParser();

		// No id, from StructureIO
		String file = this.getClass().getResource(noid).getPath();
		// Structure s = parser.parsePDBFile(file);
		Structure s = StructureIO.getStructure(file);
		assertNull(s.getPdbId());
		assertEquals(file, s.getIdentifier());
		assertEquals("", s.getName()); // differs from CIF behavior
		StructureName sid = (StructureName) s.getStructureIdentifier();
		assertEquals(file, sid.getIdentifier());

		// No id, from stream
		InputStream stream = this.getClass().getResourceAsStream(noid);
		s = parser.parsePDBFile(new BufferedReader(new InputStreamReader(stream)));
		assertNull(s.getPdbId());
		assertEquals("", s.getIdentifier());
		assertEquals("", s.getName());
		assertNull(s.getStructureIdentifier());
		
		// With id, from StructureIO
		file = this.getClass().getResource(yesid).getPath();
		s = StructureIO.getStructure(file);
		assertEquals("5PTI", s.getPdbId().toString());
		assertEquals(file, s.getIdentifier());
		assertEquals("5PTI", s.getName());
		sid = (StructureName) s.getStructureIdentifier();
		assertEquals(file, sid.getIdentifier());


		// With id, from stream
		stream = this.getClass().getResourceAsStream(yesid);
		s = parser.parsePDBFile(new BufferedReader(new InputStreamReader(stream)));
		assertEquals("5PTI", s.getPdbId().toString());
		assertEquals("", s.getIdentifier());
		assertEquals("", s.getName());
		assertNull(s.getStructureIdentifier());
	}
}