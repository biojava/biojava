/*
 *                  BioJava development code
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
 * Created on Jun 8, 2007
 *
 */
package org.biojava.nbio.structure.test;

import junit.framework.TestCase;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;

import java.io.IOException;
import java.io.InputStream;

import static org.junit.Assume.assumeNoException;

public class StructureToolsTest extends TestCase {

	Structure structure, structure2, structure3, structure4;

	@Override
	protected void setUp() throws IOException
	{
		InputStream inStream = this.getClass().getResourceAsStream("/5pti.pdb");
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();

		structure = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(structure);

		assertEquals("structure does not contain one chain ", 1 ,structure.size());

		// since biojava 5, chains contain either only polymers or only nonpolymers: here we get the first protein chain with 58 residues
		Chain chain = structure.getChainByIndex(0);
		assertEquals("Wrong number of residues.",58,chain.getAtomLength());

		inStream.close();

		// Load structure2
		inStream = this.getClass().getResourceAsStream("/1lnl.pdb");
		assertNotNull(inStream);

		structure2 = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(structure2);

		assertEquals("structure does not contain 3 chains ", 3 ,structure2.size());

		inStream.close();

		// Load structure3
		inStream = this.getClass().getResourceAsStream("/1a4w.pdb");
		assertNotNull(inStream);

		structure3 = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(structure3);

		assertEquals("structure does not contain 3 chains ", 3 ,structure3.size());

		inStream.close();
	}


	public void testGetCAAtoms(){
		Atom[] cas = StructureTools.getRepresentativeAtomArray(structure);
		assertEquals("did not find the expected number of Atoms (58), but got " + cas.length,58,cas.length);
	}

	public void testGetAtomsConsistency() throws IOException, StructureException{

		//Save the existing ChemCompProvider
		ChemCompProvider provider = ChemCompGroupFactory.getChemCompProvider();
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		cache.setFileParsingParams(params);

		Structure hivA = cache.getStructure("1hiv.A");
		Atom[] caSa = StructureTools.getRepresentativeAtomArray(hivA);
		Atom[] caCa = StructureTools.getRepresentativeAtomArray(hivA.getChainByIndex(0));
		assertEquals("did not find the same number of Atoms from structure and from chain..",
				caSa.length,caCa.length);
		Structure hivB = cache.getStructure("1hiv.B");
		Atom[] caSb = StructureTools.getRepresentativeAtomArray(hivB);
		Atom[] caCb = StructureTools.getRepresentativeAtomArray(hivB.getChainByIndex(0));
		assertEquals("did not find the same number of Atoms from structure and from chain..",
				caSb.length,caCb.length);
		//Both chains have to be the same size (A and B)
		assertEquals(99,caSa.length);
		assertEquals("did not find the same number of Atoms in both chains...",
				caSa.length,caCb.length);
		assertEquals(99,caSa.length);

		ChemCompGroupFactory.setChemCompProvider(provider);
	}

	public void testGetNrAtoms(){
		int length = StructureTools.getNrAtoms(structure);
		assertEquals("did not find the expected number of Atoms (1087), but got " + length,1087,length);
	}

	@SuppressWarnings("deprecation")
	public void testGetSubRanges() throws StructureException {
		String range;
		Structure substr;
		Chain chain;

		// normal substructures
		range = "A:3-7";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 5, chain.getAtomLength() );

		// full chains
		range = "A";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		// since biojava 5, chains contain either only polymers or only nonpolymers: here we get the first protein chain with 408 residues
		assertEquals("Did not find the expected number of residues in "+range, 408, chain.getAtomLength() );
		//assertEquals("subrange doesn't equal original chain A.", structure2.getChainByPDB("A"), chain);

		// full chains
		range = "A:";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		// since biojava 5, chains contain either only polymers or only nonpolymers: here we get the first protein chain with 408 residues
		assertEquals("Did not find the expected number of residues in "+range, 408, chain.getAtomLength() );
		//assertEquals("subrange doesn't equal original chain A.", structure2.getChainByPDB("A"), chain);

		// combined ranges
		range = "A:3-7,B:8-12";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 2, substr.size());

		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in first chain of "+range, 5, chain.getAtomLength() );

		chain = substr.getChainByIndex(1);
		assertEquals("Did not find the expected number of residues in second chain of "+range, 5, chain.getAtomLength() );

		// combined ranges
		range = "A,B:8-12";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 2, substr.size());

		// since biojava 5, chains contain either only polymers or only nonpolymers: here we get the first protein chain with 408 residues
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in first chain of "+range, 408, chain.getAtomLength() );

		chain = substr.getChainByIndex(1);
		assertEquals("Did not find the expected number of residues in second chain of "+range, 5, chain.getAtomLength() );

		// parentheses
		range = "(A:3-7)";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 5, chain.getAtomLength() );

		// single residue
		range = "A:3";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 1, chain.getAtomLength() );

		
		// negative residues
		range = "A:-3";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 1, chain.getAtomLength() );
		range = "A:-3--1";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );
		range = "A:-3-1";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 4, chain.getAtomLength() );
		
		// Special '-' case
		range = "-";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Should have gotten whole structure",structure2, substr);

		// Test single-chain syntax
		range = "_:";
		substr = StructureTools.getSubRanges(structure, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		// since biojava 5, chains contain either only polymers or only nonpolymers: here we get the first protein chain with 58 residues
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in first chain of "+range, 58, chain.getAtomLength() );

		// Test single-chain syntax in a multi-chain structure. Should give chain A.
		range = "_:";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		// since biojava 5, chains contain either only polymers or only nonpolymers: here we get the first protein chain with 408 residues
		chain = substr.getChainByIndex(0);
		assertEquals("Chain _ not converted to chain A.","A",chain.getChainID());
		assertEquals("Did not find the expected number of residues in first chain of "+range, 408, chain.getAtomLength() );

		try {
			// Illegal chain name
			range = "X:";
			substr = StructureTools.getSubRanges(structure2, range);
			fail("Illegal chain name in '"+range+"'. Should throw StructureException");
		} catch(StructureException ex) {} //expected

		// some negative tests
		try {
			range = "7-10";
			substr = StructureTools.getSubRanges(structure2, range);
			fail("Illegal range '"+range+"'. Should throw IllegalArgumentException");
		} catch(IllegalArgumentException ex) {} //expected
		try {
			range = "A7-10";
			substr = StructureTools.getSubRanges(structure2, range);
			fail("Illegal range '"+range+"'. Should throw IllegalArgumentException");
		} catch(IllegalArgumentException ex) {} //expected
	}

	public void testRevisedConvention() throws IOException, StructureException{

		AtomCache cache = new AtomCache();


		String name11 = "4hhb.A";
		Structure s = cache.getStructure(name11);
		assertEquals(1,s.getPolyChains().size());
		assertEquals(3,s.getChains().size()); // protein, HEM, water


		String name12 = "4hhb.A:";
		s = cache.getStructure(name12);
		assertEquals(1,s.getPolyChains().size());
		assertEquals(3,s.getChains().size());

		String name13 = "4hhb.A_";
		s = cache.getStructure(name13);
		assertEquals(1,s.getPolyChains().size());
		assertEquals(3,s.getChains().size());

		String name9 = "4hhb.C_1-83";
		String chainId = "C";
		s = cache.getStructure(name9);
		assertEquals(1,s.getPolyChains().size());
		assertEquals(2,s.getChains().size()); // drops waters

		Chain c = s.getPolyChainByPDB(chainId);
		assertEquals(c.getName(),chainId);
		Atom[] ca = StructureTools.getRepresentativeAtomArray(s);
		assertEquals(83,ca.length);

		String name10 = "4hhb.C_1-83,A_1-10";
		s = cache.getStructure(name10);
		assertEquals(2,s.getPolyChains().size());
		assertEquals(3,s.getChains().size()); // Includes C heme
		ca = StructureTools.getRepresentativeAtomArray(s);
		assertEquals(93, ca.length);


	}


	// this will get replaced by #81
	//	public void testStructureToolsRegexp(){
	//
	//
	//		Pattern p =  ResidueRange.RANGE_REGEX;
	//
	//		String t2 = "A_10-20";
	//		Matcher m2 = p.matcher(t2);
	//		assertNotNull(m2);
	//		assertTrue(m2.find());
	//		assertTrue(m2.matches());
	//
	//		//	for (int i=0;i< m2.groupCount();i++){
	//		//		String s = m2.group(i);
	//		//		System.out.println(s);
	//		//	}
	//		assertEquals(3,m2.groupCount());
	//
	//
	//		String t1 = "A:10-20";
	//		Matcher m1  = p.matcher(t1);
	//		assertNotNull(m1);
	//		assertTrue(m1.find());
	//		assertTrue(m1.matches());
	//		assertEquals(3,m1.groupCount());
	//
	//
	//		String t3 = "A";
	//		Matcher m3  = p.matcher(t3);
	//
	//		assertNotNull(m3);
	//		assertTrue(m3.find());
	//		assertTrue(m3.matches());
	//		assertEquals(3,m3.groupCount());
	//
	//
	//	}

	/**
	 * Test some subranges that we used to have problems with
	 * @throws StructureException
	 */
	@SuppressWarnings("deprecation")
	public void testGetSubRangesExtended() throws StructureException {
		String range;
		Structure substr;
		Chain chain;

		// negative indices
		range = "A:-3-7";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		// Note residue 0 is missing from 1lnl
		assertEquals("Did not find the expected number of residues in "+range, 10, chain.getAtomLength() );

		// double negative indices
		range = "A:-3--1";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );

		// mixed indices
		range = "A:-3-+1";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 4, chain.getAtomLength() );

		// positive indices
		range = "A:+1-6";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 6, chain.getAtomLength() );

		// partial ranges
		range = "A:-+1";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 4, chain.getAtomLength() );
		
		range = "A:--1";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );

		range = "A:^-+1";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 4, chain.getAtomLength() );
		
		range = "A:^-$";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.getPolyChains().size());
		chain = substr.getPolyChains().get(0);
		
		range = "A:400-";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 6, chain.getAtomLength() );

		range = "A:400-$";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());
		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in "+range, 6, chain.getAtomLength() );

		// whitespace
		range = "A:3-7, B:8-12";
		substr = StructureTools.getSubRanges(structure2, range);
		assertEquals("Wrong number of chains in "+range, 2, substr.size());

		chain = substr.getChainByIndex(0);
		assertEquals("Did not find the expected number of residues in first chain of "+range, 5, chain.getAtomLength() );

		chain = substr.getChainByIndex(1);
		assertEquals("Did not find the expected number of residues in second chain of "+range, 5, chain.getAtomLength() );

	}

	/**
	 * Test insertion codes
	 * @throws StructureException
	 */
	@SuppressWarnings("deprecation")
	public void testGetSubRangesInsertionCodes() throws StructureException {
		String range;
		Structure substr;
		Chain chain;

		// range including insertion
		range = "H:35-37"; //includes 36A
		substr = StructureTools.getSubRanges(structure3, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 4, chain.getAtomLength() );


		// end with insertion
		range = "H:35-36A";
		substr = StructureTools.getSubRanges(structure3, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );

		// begin with insertion
		range = "H:36A-38"; //includes 36A
		substr = StructureTools.getSubRanges(structure3, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 3, chain.getAtomLength() );

		// within insertion
		range = "L:14-14K";
		substr = StructureTools.getSubRanges(structure3, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 12, chain.getAtomLength() );

		// within insertion
		range = "L:14C-14J";
		substr = StructureTools.getSubRanges(structure3, range);
		assertEquals("Wrong number of chains in "+range, 1, substr.size());

		chain = substr.getChainByIndex(0);

		assertEquals("Did not find the expected number of residues in "+range, 8, chain.getAtomLength() );
	}

	public void testGroupsWithinShell() {
		//TODO
	}

	public void testCAmmCIF() throws StructureException {

		//Save the existing ChemCompProvider
		ChemCompProvider provider = ChemCompGroupFactory.getChemCompProvider();
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		//mmCIF files left justify their atom names (eg "CA  "), so can have different behavior
		AtomCache pdbCache = new AtomCache();
		pdbCache.setUseMmCif(false);
		FileParsingParameters params = new FileParsingParameters();
		pdbCache.setFileParsingParams(params);

		AtomCache mmcifCache = new AtomCache();
		mmcifCache.setUseMmCif(true);
		FileParsingParameters params2 = new FileParsingParameters();
		mmcifCache.setFileParsingParams(params2);


		Structure pdb=null, mmcif=null;

		String name = "3PIU";
		try {
			pdb = pdbCache.getStructure(name);
			mmcif = mmcifCache.getStructure(name);
		} catch (IOException e) {
			assumeNoException(e);
		}

		Atom[] pdbCA = StructureTools.getRepresentativeAtomArray(pdb);
		Atom[] mmcifCA = StructureTools.getRepresentativeAtomArray(mmcif);

		assertEquals("PDB has wrong length",409,pdbCA.length);
		assertEquals("PDB has wrong length",409,mmcifCA.length);

		ChemCompGroupFactory.setChemCompProvider(provider);
	}

}
