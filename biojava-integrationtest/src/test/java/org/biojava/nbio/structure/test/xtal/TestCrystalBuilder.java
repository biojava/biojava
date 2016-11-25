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
package org.biojava.nbio.structure.test.xtal;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.biojava.nbio.structure.xtal.CrystalBuilder;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

public class TestCrystalBuilder {

	@Test
	public void test1NMR() throws IOException, StructureException {

		// a monomer NMR entry: must have no interfaces

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure s1 = StructureIO.getStructure("1NMR");

		CrystalBuilder cb = new CrystalBuilder(s1);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()==0);

	}

	@Test
	public void test1B8G() throws IOException, StructureException {

		// a crystallographic entry: several interfaces

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure s1 = StructureIO.getStructure("1B8G");
		CrystalBuilder cb = new CrystalBuilder(s1);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()>1);


	}

	@Test
	public void test2MFZ() throws IOException, StructureException {

		// a dimer NMR entry: must have 1 interface

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure s1 = StructureIO.getStructure("2MFZ");
		CrystalBuilder cb = new CrystalBuilder(s1);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()==1);

	}

	@Test
	public void test4MF8() throws IOException, StructureException {

		// a crystallographic structure with protein+DNA: has only 3 prot-prot interfaces the rest are DNA-involving ones

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure s1 = StructureIO.getStructure("4MF8");
		CrystalBuilder cb = new CrystalBuilder(s1);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()>3);

	}

	@Test
	public void test2H2Z() throws IOException, StructureException {

		// a crystallographic structure C 1 2 1.
		// Should have a minimum number of contacts of 3, from the C number given in:
		// Wukowitz & Yeates, Nature Structural Biology, 1995
		// the molecule happens to be placed quite far from the origin, so this tests if we really capture all contacts

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure s1 = StructureIO.getStructure("2H2Z");
		CrystalBuilder cb = new CrystalBuilder(s1);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		assertTrue(interfaces.size()>=3);

	}
	
	@Test
	public void test4HHB() throws IOException, StructureException {

		// 4hhb is a very old entry with a non-standard coordinate frame convention, we should calculate only AU contacts
		
		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure s1 = StructureIO.getStructure("4HHB");
		CrystalBuilder cb = new CrystalBuilder(s1);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		// 5 interfaces in the AU: the 4 of the tetramer + 1 cross-interface 
		assertEquals(5, interfaces.size());

	}

}
