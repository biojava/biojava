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
 * Created on 2013-06-24
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.rcsb;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.List;

import org.junit.Test;



public class RCSBLigandsFactoryTest {


	private static final String TEST_DIR = "src/test/resources/";
	
	/**
	 * Opens the file as a {@link FileInputStream}. Copied from ResourceList, which is not in biojava.
	 */
	private FileInputStream openStream(String filename) {
		File file = new File(TEST_DIR + filename);
		FileInputStream fis;
		try {
			fis = new FileInputStream(file);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
		return fis;
	}

	/**
	 * Tests on the live database. Just makes sure the resource can be found.
	 * If this test fails, it may be because the database changed.
	 */
	@Test
	public void testUrl() {
		RCSBLigands ligands = RCSBLigandsFactory.get("1w0p");
		assertEquals(4, ligands.getLigands().size());
		assertEquals("CA", ligands.getLigands().get(0).getId());
	}

	/**
	 * Covers all the basic features, including EC numbers. Does not cover multiple polymers or multiple chains.
	 */
	@Test
	public void test1() {
		RCSBLigands description = RCSBLigandsFactory.get(openStream("describeMol/4hhb_ligands.xml"));
		
		assertEquals("4HHB", description.getPdbId());
		List<RCSBLigand> ligands = description.getLigands();
		assertEquals(2, ligands.size());
		
		RCSBLigand ligand;
		
		ligand = ligands.get(0);
		assertEquals("HEM", ligand.getId());
		assertEquals("non-polymer", ligand.getType());
		assertEquals(616.487, ligand.getWeight(), 0.0);
		assertEquals("PROTOPORPHYRIN IX CONTAINING FE", ligand.getName());
		assertEquals("C34 H32 FE N4 O4", ligand.getFormula());
		assertEquals("FEDYMSUPMFCVOD-UJJXFSCMSA-N", ligand.getInChIKey());
		assertEquals("InChI=1S/C34H34N4O4/c1-7-21-17(3)25-13-26-19(5)23(9-11-33(39)40)31(37-26)16-32-24(10-12-34(41)42)20(6)28(38-32)15-30-22(8-2)18(4)27(36-30)14-29(21)35-25/h7-8,13-16,36-37H,1-2,9-12H2,3-6H3,(H,39,40)(H,41,42)/b25-13-,26-13-,27-14-,28-15-,29-14-,30-15-,31-16-,32-16-", ligand.getInChI());
		assertEquals("Cc1c2/cc/3\\nc(/cc\\4/c(c(/c(/[nH]4)c/c5n/c(c\\c(c1CCC(=O)O)[nH]2)/C(=C5C)CCC(=O)O)C=C)C)C(=C3C)C=C", ligand.getSmiles());

		ligand = ligands.get(1);
		assertEquals("PO4", ligand.getId());
		assertEquals("non-polymer", ligand.getType());
		assertEquals(94.971, ligand.getWeight(), 0.0);
		assertEquals("PHOSPHATE ION", ligand.getName());
		assertEquals("O4 P -3", ligand.getFormula());
		assertEquals("NBIIXXVUZAFLBC-UHFFFAOYSA-K", ligand.getInChIKey());
		assertEquals("InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-3", ligand.getInChI());
		assertEquals("[O-]P(=O)([O-])[O-]", ligand.getSmiles());
	}

}
