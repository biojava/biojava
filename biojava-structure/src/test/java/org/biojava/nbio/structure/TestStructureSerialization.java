/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import static org.junit.Assert.*;

/**
 * Test the serialization and deserialization of BioJava structure objects.
 * 
 * @author Aleix Lafita
 *
 */
public class TestStructureSerialization {

	@Test
	public void testSerializeStructure() throws IOException, StructureException, ClassNotFoundException {

		PDBFileReader reader = new PDBFileReader();
		reader.getFileParsingParameters().setParseSecStruc(true);
		Structure sin = reader.getStructure("src/test/resources/2gox.pdb");

		// Serialize the structure object and keep it in memory
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		ObjectOutputStream objectOut = new ObjectOutputStream(baos);
		objectOut.writeObject(sin);
		objectOut.close();
		byte[] bytes = baos.toByteArray();
		
		// Deserialize the bytes back into a structure object
		ByteArrayInputStream bais = new ByteArrayInputStream(bytes);
		ObjectInputStream objectIn = new ObjectInputStream(bais);
		Structure sout = (Structure) objectIn.readObject();
		objectIn.close();
		
		// Test properties of the structures before and after serialization
		assertEquals(sin.nrModels(), sout.nrModels());
		assertEquals(sin.getChains().size(), sout.getChains().size());
		assertEquals(sin.getName(), sout.getName());
		
		// Test equal string representations
		assertEquals(sin.toString(), sout.toString());
		
	}
}
