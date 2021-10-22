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

import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertSame;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;

import org.junit.jupiter.api.Test;


/**
 * Tests for {@link PdbId} parsing and its usability to convert between PDB ID formats
 * @author Amr ALHOSSARY
 * @since 6.0.0
 */
public class TestPdbId {


	@Test
	public void testGetIdInDefaultFormat() {
		PdbId pdbId;
		String id;

		pdbId = new PdbId("1abc");
		id = pdbId.getId();
		assertEquals(id, "1ABC");

		pdbId = new PdbId("PDB_55551abc");
		id = pdbId.getId();
		assertEquals(id, "PDB_55551ABC");
	}

	@Test
	public void testGetIdPrefereShortFormat() {
		PdbId pdbId;
		String id;

		pdbId = new PdbId("1abc");
		id = pdbId.getId(PdbId.Behavior.PREFER_SHORT);
		assertEquals(id, "1ABC");

		pdbId = new PdbId("PDB_55551abc");
		id = pdbId.getId(PdbId.Behavior.PREFER_SHORT);
		assertEquals(id, "PDB_55551ABC");
	}

	@Test
	public void testGetIdPrefereExtendedFormat() {
		PdbId pdbId;
		String id;

		pdbId = new PdbId("1abc");
		id = pdbId.getId(PdbId.Behavior.PREFER_EXTENDED);
		assertEquals(id, "PDB_00001ABC");

		pdbId = new PdbId("PDB_55551abc");
		id = pdbId.getId(PdbId.Behavior.PREFER_EXTENDED);
		assertEquals(id, "PDB_55551ABC");
	}
	
	@Test
	public void testGetIdInShortFormat() {
		assertDoesNotThrow(() -> {
			PdbId pdbId = new PdbId("1abc");
			String id = pdbId.getShortId();
			assertEquals(id, "1ABC");
		}, "Unexpected Exception thrown");

		assertThrows(StructureException.class, () -> {
			PdbId pdbId = new PdbId("PDB_55551abc");
			pdbId.getShortId();
		}, "wrongly shortened a non-shortable ID");
	}
	
	
	@Test
	public void testIsShortPDBID() {
		assertTrue(PdbId.isValidShortPdbId("1abc"), "Didn't accept lower case");
		assertTrue(PdbId.isValidShortPdbId("4HHB"), "Didn't accept upper case");
		assertFalse(PdbId.isValidShortPdbId("HHHB"), "Accepted wrong format");
		assertFalse(PdbId.isValidShortPdbId("PDB_00001ABC"), "Accepted extended format");
	}
	
	@Test
	public void testIsExtendedPDBID() {
		assertTrue(PdbId.isValidExtendedPdbId("PDB_00001abc"), "Didn't accept lower case");
		assertTrue(PdbId.isValidExtendedPdbId("PDB_00004HHB"), "Didn't accept upper case");
		assertTrue(PdbId.isValidExtendedPdbId("PDB_22224HHB"), "Didn't accept upper case");
		assertFalse(PdbId.isValidExtendedPdbId("PDB_AAAA4HHB"), "Accepted wrong format");
		assertFalse(PdbId.isValidExtendedPdbId("1ABC"), "Accepted short format");
	}

	@Test
	public void testIsShortCompatible() {
		assertTrue(PdbId.isShortCompatible("PDB_00001abc"), "Didn't accept lower case");
		assertTrue(PdbId.isShortCompatible("PDB_00004HHB"), "Didn't accept upper case");
		assertFalse(PdbId.isShortCompatible("1ABC"), "Accepted short format");
		assertFalse(PdbId.isShortCompatible("PDB_AAAA4HHB"), "Accepted wrong format");

		//Although this is wrong, returning true is the expected behavior of 
		// this method; because it does NOT validate the passed in string.
		assertTrue(PdbId.isShortCompatible("PDB_0000XXXXXXXXXXXXX"), "Accepted wrong format");
	}
	
	@Test
	public void testToExtendedFormat() {
		assertDoesNotThrow(() -> {
			assertEquals(PdbId.toExtendedId("1abc"), "PDB_00001ABC");
		}, "Couldn't extend Id");

		assertDoesNotThrow(() -> {
			assertEquals(PdbId.toExtendedId("PDB_00001abc"), "PDB_00001ABC");
		}, "Didn't recognize extended format");
		
		assertThrows(StructureException.class, () -> {
			PdbId.toExtendedId("PDB_aaaa1abc");
		}, "Accepted wrong format");
	}
	
	@Test
	public void testToShortFormat() {
		assertDoesNotThrow(() -> {
			assertEquals(PdbId.toShortId("PDB_00001ABC"), "1ABC");
		}, "Couldn't shorten Id");
		
		assertDoesNotThrow(() -> {
			assertEquals(PdbId.toShortId("1abc"), "1ABC");
		}, "Didn't recognize short format");
		
		assertThrows(StructureException.class, () -> {
			PdbId.toShortId("PDB_aaaa1abc");
		}, "Accepted wrong format");
		
		assertThrows(StructureException.class, () -> {
			PdbId.toShortId("aabc");
		}, "Accepted wrong format");
	}
	
	@Test
	public void testHashCodeAndEquals() {
		PdbId id1, id2, id3/* , id4 */;
		PdbId other;
		id1 = new PdbId("1abc");
		id2 = new PdbId("PDB_00001ABC");
		id3 = new PdbId("1ABC");
//		id4 = new PdbId("pdb_00001abc");
		other = new PdbId("2ABC");
		
		assertEquals(id1.hashCode(), id2.hashCode());
		assertEquals(id1.hashCode(), id3.hashCode());
//		assertEquals(id1.hashCode(), id4.hashCode());
		assertNotEquals(id1.hashCode(), other.hashCode());
		
		assertTrue(id1.equals(id2));
		assertTrue(id1.equals(id3));
//		assertTrue(id1.equals(id4));
		assertFalse(id1.equals(other));
	}

	@Test
	public void testXXXX() {
		PdbId x1 = new PdbId(PdbId.XXXX_STRING);
		PdbId x2 = new PdbId(PdbId.XXXX_STRING);
		assertFalse(x1 == x2);
		assertEquals(x1.getId(), x2.getId());
		assertEquals(x1.hashCode(), x2.hashCode());
		assertEquals(x1, x1);
		assertEquals(x2, x2);
		assertNotEquals(x1, x2);
	}

	@Test
	public void testClone() {
		assertDoesNotThrow(() -> {
			PdbId id1 = new PdbId("1abc");
			PdbId clone = (PdbId) id1.clone();
			
			assertNotSame(id1, clone);
			assertEquals(id1, clone);
			assertEquals(id1.hashCode(), clone.hashCode());
		}, "unexpected exception thrown while cloning");
	}
	

	@Test
	public void testCompareTo() {
		PdbId id1, id2, id3, id4, id5 ;
		PdbId[] array, expected;
		id1 = new PdbId("1abc");
		id2 = new PdbId("PDB_00011ABC");
		id3 = new PdbId("2ABC");
		id4 = new PdbId("PDB_00001ABA");
		id5 = new PdbId("1100");
		
		array = new PdbId[] {id1, id2, id3, id4, id5};
		System.out.println(Arrays.deepToString(array));
		Arrays.sort(array);
		System.out.println(Arrays.deepToString(array));
		expected = new PdbId[] {id5, id4, id1, id3, id2};
		System.out.println(Arrays.deepToString(expected));
		assertArrayEquals(expected, array);
		
		
		//let's try to have some "distinct but equal" objects.
		id1 = new PdbId("1abc");
		id2 = new PdbId("PDB_00011ABC");
		id3 = new PdbId("2ABC");
		id4 = new PdbId("PDB_00001ABA");
		id5 = new PdbId("1ABA");
		
		array = new PdbId[] {id1, id2, id3, id4, id5};
//		System.out.println(Arrays.deepToString(array));
		Arrays.sort(array);
//		System.out.println(Arrays.deepToString(array));
		expected = new PdbId[] {id5, id4, id1, id3, id2};
//		System.out.println(Arrays.deepToString(expected));
		assertArrayEquals(expected, array); // They should be.
		//Now let the real test begins
		for (int i = 0; i < 2; i++) {
			assertNotSame("Couldn't detect 2 objects that are equal but not the same", expected[i], array[i]);
		}
		for (int i = 2; i < expected.length; i++) {
			assertSame(expected[i], array[i]);
		}
	}
}
