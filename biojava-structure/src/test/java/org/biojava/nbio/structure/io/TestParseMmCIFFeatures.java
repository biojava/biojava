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
package org.biojava.nbio.structure.io;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

/**
 * Created by larsonmattr on 10/31/2015.
 */
public class TestParseMmCIFFeatures {
	@Test
	public void testSSBond()throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		FileParsingParameters params = new FileParsingParameters();
		params.setCreateAtomBonds(true);
		cache.setFileParsingParams(params);
		Structure sCif = StructureIO.getStructure("2OYA");

		assertNotNull(sCif);

		// After it has read the file, it should check that expected SSBONDs are present.
		List<Bond> bonds = sCif.getSSBonds();

		// 2OYA has 6 ssbonds, 3 on molecule A and 3 on molecule B
		assertEquals(6, bonds.size());



		// Check the bonds
		assertDisulfideBond("A", "A", 446, 507, bonds.get(0));
		assertDisulfideBond("A", "A", 459, 517, bonds.get(1));
		assertDisulfideBond("A", "A", 487, 497, bonds.get(2));
		assertDisulfideBond("B", "B", 446, 507, bonds.get(3));
		assertDisulfideBond("B", "B", 459, 517, bonds.get(4));
		assertDisulfideBond("B", "B", 487, 497, bonds.get(5));

	}

	@Test
	public void testSSBondAltLocs() throws IOException, StructureException {
		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		FileParsingParameters params = new FileParsingParameters();
		params.setCreateAtomBonds(true);
		cache.setFileParsingParams(params);
		Structure sCif = StructureIO.getStructure("3DVF");

		assertNotNull(sCif);

		// After it has read the file, it should check that expected SSBONDs are present.
		List<Bond> bonds = sCif.getSSBonds();

		// 3DVF has 2 ssbonds including 1 in an alt loc
		assertEquals(2, bonds.size());


		// Check the bonds
		assertDisulfideBond("A", "A", 23, 88, bonds.get(0));
		assertDisulfideBond("A", "A", 23, 88, bonds.get(1));

		// check that we have the bonds correctly assigned to both atom and its altgroup atom
		Group g = sCif.getPolyChainByPDB("A").getGroupByPDB(new ResidueNumber("A", 23, ' '));
		Group altG = g.getAltLocGroup('B');
		Group g2 = sCif.getPolyChainByPDB("A").getGroupByPDB(new ResidueNumber("A", 88, ' '));
		List<Bond> bs = g.getAtom("SG").getBonds();
		Bond b = null;
		for (Bond bn:bs){
			if (bn.getAtomA().getName().equals("SG") && bn.getAtomB().getName().equals("SG")) {
				b = bn;
			}
		}
		List<Bond> altbs = altG.getAtom("SG").getBonds();
		Bond bAltG = null;
		for (Bond bn:altbs){
			if (bn.getAtomA().getName().equals("SG") && bn.getAtomB().getName().equals("SG")) {
				bAltG = bn;
			}
		}

		assertNotNull(b);
		assertNotNull(bAltG);

		//System.out.println(g);
		//System.out.println(altG);
		assertSame(g.getAtom("SG") , b.getAtomA());
		assertSame(altG.getAtom("SG"), bAltG.getAtomA());
		assertSame(g2.getAtom("SG"), b.getAtomB());
		assertSame(g2.getAtom("SG"), bAltG.getAtomB());


	}

	private void assertDisulfideBond(String expectedChainId1, String expectedChainId2, int expectedResSerial1, int expectedResSerial2, Bond bond) {
		String chainId1 = bond.getAtomA().getGroup().getChainId();
		String chainId2 = bond.getAtomB().getGroup().getChainId();
		ResidueNumber resNum1 = bond.getAtomA().getGroup().getResidueNumber();
		ResidueNumber resNum2 = bond.getAtomB().getGroup().getResidueNumber();
		assertEquals("disulfide bond first chain id failed ", expectedChainId1, chainId1);
		assertEquals("disulfide bond second chain id failed ", expectedChainId2, chainId2);
		assertEquals("disulfide bond failed first residue number failed ", new ResidueNumber(expectedChainId1, expectedResSerial1, null), resNum1);
		assertEquals("disulfide bond failed second residue number failed ", new ResidueNumber(expectedChainId2, expectedResSerial2, null), resNum2);
	}


	@Test
	public void testSites()throws IOException, StructureException {

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("4HHB");

		assertNotNull(sCif);

		// After it has read the file, it should check that expected SITES are present.
		List<Site> sites = sCif.getSites();

		// 4HHB has 6 sites from ligands.
		assertEquals(6, sites.size());

		// Check for each site that it has parsed all residues.
		assertEquals(1, getGroupsInSite(sCif, "AC1"));
		assertEquals(1, getGroupsInSite(sCif, "AC2"));
		assertEquals(16, getGroupsInSite(sCif, "AC3"));
		assertEquals(13, getGroupsInSite(sCif, "AC4"));
		assertEquals(15, getGroupsInSite(sCif, "AC5"));
		assertEquals(7, getGroupsInSite(sCif, "AC6"));

		// Check that struct_site parsing worked, and they have descriptions.
		assertEquals(getDescription(sCif, "AC1"), "BINDING SITE FOR RESIDUE PO4 D 147");
		assertEquals(getDescription(sCif, "AC2"), "BINDING SITE FOR RESIDUE PO4 B 147");
		assertEquals(getDescription(sCif, "AC3"), "BINDING SITE FOR RESIDUE HEM A 142");
		assertEquals(getDescription(sCif, "AC4"), "BINDING SITE FOR RESIDUE HEM B 148");
		assertEquals(getDescription(sCif, "AC5"), "BINDING SITE FOR RESIDUE HEM C 142");
		assertEquals(getDescription(sCif, "AC6"), "BINDING SITE FOR RESIDUE HEM D 148");
	}

	@Test
	public void testSites1a4w()throws IOException, StructureException {
		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("1A4W");

		assertNotNull(sCif);

		// After it has read the file, it should check that expected SITES are present.
		List<Site> sites = sCif.getSites();

		// 1a4w has 5 sites from ligands.
		assertEquals(5, sites.size());

		// Check for each site that it has parsed all residues.
		assertEquals(3, getGroupsInSite(sCif, "CAT"));
		assertEquals(6, getGroupsInSite(sCif, "AC1")); // Site has residue with insertion code.
		assertEquals(6, getGroupsInSite(sCif, "AC2"));
		assertEquals(14, getGroupsInSite(sCif, "AC3"));
		assertEquals(14, getGroupsInSite(sCif, "AC4"));

	}

	private int getGroupsInSite(Structure structure, String site) {
		for (Site a_site : structure.getSites()) {
			if (a_site.getSiteID().equals(site)) {
				return a_site.getGroups().size();
			}
		}
		return 0;
	}

	private String getDescription(Structure structure, String site) {
		for (Site a_site : structure.getSites()) {
			if (a_site.getSiteID().equals(site)) {
				return a_site.getDescription();
			}
		}
		return "";
	}
}
