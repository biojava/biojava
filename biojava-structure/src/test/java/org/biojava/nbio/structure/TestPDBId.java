package org.biojava.nbio.structure;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

public class TestPDBId {


	@Test
	public void testGetIdInDefaultFormat() {
		PDBId pdbId;
		String id;

		pdbId = new PDBId("1abc");
		id = pdbId.getId();
		assertEquals(id, "1ABC");

		pdbId = new PDBId("PDB_55551abc");
		id = pdbId.getId();
		assertEquals(id, "PDB_55551ABC");
	}

	@Test
	public void testGetIdPrefereShortFormat() {
		PDBId pdbId;
		String id;

		pdbId = new PDBId("1abc");
		id = pdbId.getId(PDBId.Behavior.PREFER_SHORT);
		assertEquals(id, "1ABC");

		pdbId = new PDBId("PDB_55551abc");
		id = pdbId.getId(PDBId.Behavior.PREFER_SHORT);
		assertEquals(id, "PDB_55551ABC");
	}

	@Test
	public void testGetIdPrefereExtendedFormat() {
		PDBId pdbId;
		String id;

		pdbId = new PDBId("1abc");
		id = pdbId.getId(PDBId.Behavior.PREFER_EXTENDED);
		assertEquals(id, "PDB_00001ABC");

		pdbId = new PDBId("PDB_55551abc");
		id = pdbId.getId(PDBId.Behavior.PREFER_EXTENDED);
		assertEquals(id, "PDB_55551ABC");
	}
	
	@Test
	public void testGetIdInShortFormat() {
		try {
			PDBId pdbId = new PDBId("1abc");
			String id = pdbId.getShortId();
			assertEquals(id, "1ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Unexpected Exception thrown");
		}
		
		try {
			PDBId pdbId = new PDBId("PDB_55551abc");
			pdbId.getShortId();
			fail("wrongly shortened a non-shortable ID");
		} catch (StructureException e) {}
	}
	
	
	@Test
	public void testIsShortPDBID() {
		assertTrue("Didn't accept lower case", PDBId.isValidShortPDBID("1abc"));
		assertTrue("Didn't accept upper case", PDBId.isValidShortPDBID("4HHB"));
		assertFalse("Accepted wrong format", PDBId.isValidShortPDBID("HHHB"));
		assertFalse("Accepted extended format", PDBId.isValidShortPDBID("PDB_00001ABC"));
	}
	
	@Test
	public void testIsExtendedPDBID() {
		assertTrue("Didn't accept lower case", PDBId.isValidExtendedPDBID("PDB_00001abc"));
		assertTrue("Didn't accept upper case", PDBId.isValidExtendedPDBID("PDB_00004HHB"));
		assertTrue("Didn't accept upper case", PDBId.isValidExtendedPDBID("PDB_22224HHB"));
		assertFalse("Accepted wrong format", PDBId.isValidExtendedPDBID("PDB_AAAA4HHB"));
		assertFalse("Accepted short format", PDBId.isValidExtendedPDBID("1ABC"));
	}

	@Test
	public void testIsShortCompatible() {
		assertTrue("Didn't accept lower case", PDBId.isShortCompatible("PDB_00001abc"));
		assertTrue("Didn't accept upper case", PDBId.isShortCompatible("PDB_00004HHB"));
		assertFalse("Accepted short format", PDBId.isShortCompatible("1ABC"));
		assertFalse("Accepted wrong format", PDBId.isShortCompatible("PDB_AAAA4HHB"));

		//Although this is wrong, returning true is the expected behavior of 
		// this method; because it does NOT validate the passed in string.
		assertTrue("Accepted wrong format", PDBId.isShortCompatible("PDB_0000XXXXXXXXXXXXX"));
	}
	
	@Test
	public void testToExtendedFormat() {
		try {
			assertEquals(PDBId.toExtendedId("1abc"), "PDB_00001ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Couldn't extend Id");
		}

		try {
			assertEquals(PDBId.toExtendedId("PDB_00001abc"), "PDB_00001ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Didn't recognize extended format");
		}
		
		try {
			PDBId.toExtendedId("PDB_aaaa1abc");
			fail("Accepted wrong format");
		} catch (StructureException e) {}
	}
	
	@Test
	public void testToShortFormat() {
		try {
			assertEquals(PDBId.toShortId("PDB_00001ABC"), "1ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Couldn't shorten Id");
		}
		
		try {
			assertEquals(PDBId.toShortId("1abc"), "1ABC");
		} catch (StructureException e) {
			e.printStackTrace();
			fail("Didn't recognize short format");
		}
		
		try {
			PDBId.toShortId("PDB_aaaa1abc");
			fail("Accepted wrong format");
		} catch (StructureException e) {}
		
		try {
			PDBId.toShortId("aabc");
			fail("Accepted wrong format");
		} catch (StructureException e) {}
	}
	
	@Test
	public void testHashCodeAndEquals() {
		PDBId id1, id2, id3/* , id4 */;
		PDBId other;
		id1 = new PDBId("1abc");
		id2 = new PDBId("PDB_00001ABC");
		id3 = new PDBId("1ABC");
//		id4 = new PDBId("pdb_00001abc");
		other = new PDBId("2ABC");
		
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
		PDBId id1, clone;
		try {
			id1 = new PDBId("1abc");
			clone = (PDBId) id1.clone();
			
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
		PDBId id1, id2, id3, id4, id5 ;
		PDBId[] array, expected;
		id1 = new PDBId("1abc");
		id2 = new PDBId("PDB_00011ABC");
		id3 = new PDBId("2ABC");
		id4 = new PDBId("PDB_00001ABA");
		id5 = new PDBId("1100");
		
		array = new PDBId[] {id1, id2, id3, id4, id5};
		System.out.println(Arrays.deepToString(array));
		Arrays.sort(array);
		System.out.println(Arrays.deepToString(array));
		expected = new PDBId[] {id5, id4, id1, id3, id2};
		System.out.println(Arrays.deepToString(expected));
		assertArrayEquals(expected, array);
		
		
		//let's try to have some "distinct but equal" objects.
		id1 = new PDBId("1abc");
		id2 = new PDBId("PDB_00011ABC");
		id3 = new PDBId("2ABC");
		id4 = new PDBId("PDB_00001ABA");
		id5 = new PDBId("1ABA");
		
		array = new PDBId[] {id1, id2, id3, id4, id5};
//		System.out.println(Arrays.deepToString(array));
		Arrays.sort(array);
//		System.out.println(Arrays.deepToString(array));
		expected = new PDBId[] {id5, id4, id1, id3, id2};
//		System.out.println(Arrays.deepToString(expected));
		assertArrayEquals(expected, array); // They should be.
		//Now let the real test begins
		boolean foundEqualButNotTheSame = false;
		for (int i = 0; i < expected.length; i++) {
			if(expected[i] == array[i]) {
				foundEqualButNotTheSame = true;
				break;
			}
		}
		assertTrue("Couldn't detect 2 objects that are equal but not the same", foundEqualButNotTheSame);
	}
}
