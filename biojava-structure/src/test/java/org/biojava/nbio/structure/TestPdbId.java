package org.biojava.nbio.structure;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

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
		try {
			PdbId pdbId = new PdbId("1abc");
			String id = pdbId.getShortId();
			assertEquals(id, "1ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Unexpected Exception thrown");
		}
		
		try {
			PdbId pdbId = new PdbId("PDB_55551abc");
			pdbId.getShortId();
			fail("wrongly shortened a non-shortable ID");
		} catch (StructureException e) {}
	}
	
	
	@Test
	public void testIsShortPDBID() {
		assertTrue("Didn't accept lower case", PdbId.isValidShortPdbId("1abc"));
		assertTrue("Didn't accept upper case", PdbId.isValidShortPdbId("4HHB"));
		assertFalse("Accepted wrong format", PdbId.isValidShortPdbId("HHHB"));
		assertFalse("Accepted extended format", PdbId.isValidShortPdbId("PDB_00001ABC"));
	}
	
	@Test
	public void testIsExtendedPDBID() {
		assertTrue("Didn't accept lower case", PdbId.isValidExtendedPdbId("PDB_00001abc"));
		assertTrue("Didn't accept upper case", PdbId.isValidExtendedPdbId("PDB_00004HHB"));
		assertTrue("Didn't accept upper case", PdbId.isValidExtendedPdbId("PDB_22224HHB"));
		assertFalse("Accepted wrong format", PdbId.isValidExtendedPdbId("PDB_AAAA4HHB"));
		assertFalse("Accepted short format", PdbId.isValidExtendedPdbId("1ABC"));
	}

	@Test
	public void testIsShortCompatible() {
		assertTrue("Didn't accept lower case", PdbId.isShortCompatible("PDB_00001abc"));
		assertTrue("Didn't accept upper case", PdbId.isShortCompatible("PDB_00004HHB"));
		assertFalse("Accepted short format", PdbId.isShortCompatible("1ABC"));
		assertFalse("Accepted wrong format", PdbId.isShortCompatible("PDB_AAAA4HHB"));

		//Although this is wrong, returning true is the expected behavior of 
		// this method; because it does NOT validate the passed in string.
		assertTrue("Accepted wrong format", PdbId.isShortCompatible("PDB_0000XXXXXXXXXXXXX"));
	}
	
	@Test
	public void testToExtendedFormat() {
		try {
			assertEquals(PdbId.toExtendedId("1abc"), "PDB_00001ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Couldn't extend Id");
		}

		try {
			assertEquals(PdbId.toExtendedId("PDB_00001abc"), "PDB_00001ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Didn't recognize extended format");
		}
		
		try {
			PdbId.toExtendedId("PDB_aaaa1abc");
			fail("Accepted wrong format");
		} catch (StructureException e) {}
	}
	
	@Test
	public void testToShortFormat() {
		try {
			assertEquals(PdbId.toShortId("PDB_00001ABC"), "1ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Couldn't shorten Id");
		}
		
		try {
			assertEquals(PdbId.toShortId("1abc"), "1ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Didn't recognize short format");
		}
		
		try {
			PdbId.toShortId("PDB_aaaa1abc");
			fail("Accepted wrong format");
		} catch (StructureException e) {}
		
		try {
			PdbId.toShortId("aabc");
			fail("Accepted wrong format");
		} catch (StructureException e) {}
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
	public void testClone() {
		PdbId id1, clone;
		try {
			id1 = new PdbId("1abc");
			clone = (PdbId) id1.clone();
			
			assertFalse(id1 == clone);
			assertTrue(id1.equals(clone));
			assertEquals(id1.hashCode(), clone.hashCode());
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
			fail("unexpected exception thrown while cloning");
		}
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
		boolean foundEqualButNotTheSame = false;
		for (int i = 0; i < expected.length; i++) {
			if(expected[i] != array[i]) {
				foundEqualButNotTheSame = true;
				break;
			}
		}
		assertTrue("Couldn't detect 2 objects that are equal but not the same", foundEqualButNotTheSame);
	}
}
