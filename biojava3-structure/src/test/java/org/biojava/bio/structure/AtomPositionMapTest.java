/**
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
 * Created on 2012-10-27
 * Created by dmyerstu
 *
 */
package org.biojava.bio.structure;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.NavigableMap;

import org.biojava.bio.structure.AtomPositionMap.GroupMatcher;
import org.biojava.bio.structure.align.util.AtomCache;
import org.junit.Before;
import org.junit.Test;

/**
 * A unit test for {@link AtomPositionMap}. Make sure to change the AtomCache directory in {@link #setUp()}.
 * 
 * @author dmyerstu
 * @since 3.0.6
 */
public class AtomPositionMapTest {

	private AtomCache cache;

	@Before
	public void setUp() throws Exception {
		cache = new AtomCache("/tmp/PDB", false);
	}

	/**
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testBasic() throws IOException, StructureException { // no insertion codes
		String pdbId = "1w0p";
		int length = 92;
		ResidueNumber start = new ResidueNumber("A", 25, null);
		ResidueNumber end = new ResidueNumber("A", 117, null);
		AtomPositionMap map = AtomPositionMap.ofAminoAcids(cache.getAtoms(pdbId));
		NavigableMap<ResidueNumber, Integer> navMap = map.getNavMap();
		// assertEquals("The maps have different sizes", navMap.size(), map.size());
		for (ResidueNumber n : navMap.keySet()) {
			assertEquals("An element is missing", map.getPosition(n), navMap.get(n).intValue());
		}
		int realLength = map.calcLength(start, end);
		assertEquals("Real atom length is wrong", length, realLength);
	}

	/**
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testWithCustomMatcher() throws IOException, StructureException {
		String pdbId = "1w0p";
		GroupMatcher matcher = new GroupMatcher() {
			int i = 0;

			@Override
			public boolean matches(Group group) {
				i++;
				if (group.getType().equals(GroupType.AMINOACID)) {
					AminoAcid aa = (AminoAcid) group;
					if (aa.getAminoType().equals('K')) {
						// System.out.println(i);
						return true;
					}
				}
				return false;
			}
		};
		int length = 7;
		ResidueNumber start = new ResidueNumber("A", 42, null); // Calpha = 130
		ResidueNumber end = new ResidueNumber("A", 323, null); // Calpha = 2274
		Structure structure = cache.getStructure(pdbId);
		Atom[] atoms = StructureTools.getAllAtomArray(structure);
		AtomPositionMap map = new AtomPositionMap(atoms, matcher);
		int realLength = map.calcLength(start, end);
		assertEquals("Real atom length is wrong", length, realLength);
		int lengthWithAtoms = map.calcLength(130, 2274, "A");
		assertEquals("Real atom length is wrong", length, lengthWithAtoms);
	}

	/**
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testWithInsertionCodes() throws IOException, StructureException {

		String pdbId = "1qdm";
		// has 2 insertion code regions:
		// P at the beginning starting at 6P and ending at 27P, where 2 starts
		// S between 247 (before:1S) and 248 (after:104S)

		AtomPositionMap map = AtomPositionMap.ofAminoAcids(cache.getAtoms(pdbId));
		NavigableMap<ResidueNumber, Integer> navMap = map.getNavMap();
		// assertEquals("The maps have different sizes", navMap.size(), map.size());
		for (ResidueNumber n : navMap.keySet()) {
			assertEquals("An element is missing", map.getPosition(n), navMap.get(n).intValue());
		}

		int length1 = 59;
		int length2 = 131;
		ResidueNumber start = new ResidueNumber("A", 246, null);
		ResidueNumber mid = new ResidueNumber("A", 85, 'S');
		ResidueNumber end = new ResidueNumber("A", 300, null);
		int realLength1 = map.calcLength(start, mid);
		assertEquals("Real atom length is wrong", length1, realLength1);
		int realLength2 = map.calcLength(start, end);
		assertEquals("Real atom length is wrong", length2, realLength2);
	}

}
