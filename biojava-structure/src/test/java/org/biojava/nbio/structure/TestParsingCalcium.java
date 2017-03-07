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
 * Created on Nov 1, 2013
 * Author: andreas
 *
 */

package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class TestParsingCalcium {


	@Test
	public void testCalciumParsing() throws StructureException, IOException {

		String pdbID = "1SU4";

		// Calcium is at position 995
		// HETATM 7673 CA    CA A 995      64.194  12.588   7.315  1.00 41.55          CA

		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure(pdbID);
		cache.setUseMmCif(true);

		Structure m = cache.getStructure(pdbID);

		Group g1 = s.getNonPolyChainsByPDB("A").get(0).getGroupByPDB(new ResidueNumber("A",995,null));
		Group g2 = m.getNonPolyChainsByPDB("A").get(0).getGroupByPDB(new ResidueNumber("A",995,null));

		// can't do that! the atom index is not the same!
		//assertEquals(g1.getAtom(0).toPDB(),g2.getAtom(0).toPDB());

		assertEquals(g1.getAtom(0).getName(),"CA");
		assertEquals(g1.getAtom(0).getElement(),Element.Ca);
		assertEquals(g1.getAtom(0).getName(), g2.getAtom(0).getName());


	}

	@Test
	public void testCAreturnsCalpha() throws IOException, StructureException {

		// there's an ambiguity in PDB names between the C alpha of an aminoacid, named "CA"
		// and the Calcium of a Calcium ion, named "CA" (or potentially any calcium present in a group)
		// in PDB files both are distinguished by different paddings: " CA " (Calpha) vs "CA  " (Calcium)


		String pdbID = "1SU4";

		// Calcium is at position 995
		// HETATM 7673 CA    CA A 995      64.194  12.588   7.315  1.00 41.55          CA

		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure(pdbID);

		Atom[] atoms = StructureTools.getRepresentativeAtomArray(s);
		for (Atom atom:atoms) {
			assertTrue("atom "+atom.getPDBserial()+" of residue "+atom.getGroup().getResidueNumber()+"-"+atom.getGroup().getPDBName()+
					" is not a Carbon alpha", atom.getElement()==Element.C);
		}

		for (Chain c:s.getChains()) {
			atoms = StructureTools.getRepresentativeAtomArray(c);
			for (Atom atom:atoms) {
				assertTrue("atom "+atom.getPDBserial()+" of residue "+atom.getGroup().getResidueNumber()+"-"+atom.getGroup().getPDBName()+
						" is not a Carbon alpha", atom.getElement()==Element.C);
			}
		}

		atoms = StructureTools.getBackboneAtomArray(s);

		boolean hasGlycine = false;

		for (Atom atom:atoms) {
			assertTrue("atom "+atom.getPDBserial()+" of residue "+atom.getGroup().getResidueNumber()+"-"+atom.getGroup().getPDBName()+
					" is not a backbone atom",
					atom.getElement()==Element.C ||
					atom.getElement()==Element.O ||
					atom.getElement()==Element.N    );

			assertTrue("backbone atoms should not contain CB atoms",!atom.getName().equals("CB"));

			if (atom.getGroup().getPDBName().equals("GLY")) {
				hasGlycine = true;
			}

		}
		assertTrue("the backbone atoms should contain atoms from at least 1 glycine",hasGlycine);

	}
}
