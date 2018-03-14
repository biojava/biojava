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

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;

import static org.junit.Assume.assumeNoException;
import static org.junit.Assert.*;

public class StructureToolsTest {

	Structure structure, structure2, structure3, structure4;

	@Before
	public void setUp() throws IOException
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

	@Test
	public void testGetCAAtoms(){
		Atom[] cas = StructureTools.getRepresentativeAtomArray(structure);
		assertEquals("did not find the expected number of Atoms (58), but got " + cas.length,58,cas.length);
	}

	@Test
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

	@Test
	public void testGetNrAtoms(){
		int length = StructureTools.getNrAtoms(structure);
		assertEquals("did not find the expected number of Atoms (1087), but got " + length,1087,length);
	}

	
	@Test
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

	public void testGroupsWithinShell() {
		//TODO
	}

	@Test
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
	
	@Test
	public void testGetRepresentativeAtomsProtein() throws StructureException, IOException {
		Structure s = StructureIO.getStructure("1smt");
		Chain c = s.getChainByIndex(0);
		Atom[] atoms = StructureTools.getRepresentativeAtomArray(c);
		assertEquals(98,atoms.length);
		
		Chain clonedChain = (Chain)c.clone();
		atoms = StructureTools.getRepresentativeAtomArray(clonedChain); 
		assertEquals(98,atoms.length);
	}

	/**
	 * See https://github.com/biojava/biojava/issues/631
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testGetRepresentativeAtomsDna() throws StructureException, IOException {
	
		Structure s = StructureIO.getStructure("2pvi");
		Chain c = s.getPolyChainByPDB("C");
		Atom[] atoms = StructureTools.getRepresentativeAtomArray(c); // chain C (1st nucleotide chain)
		// actually it should be 13, but at the moment one of the nucleotides is not caught correctly because it's non-standard
		assertEquals(12,atoms.length);
		
		Chain clonedChain = (Chain)c.clone();
		atoms = StructureTools.getRepresentativeAtomArray(clonedChain); // chain C (1st nucleotide chain)
		assertEquals(12,atoms.length);
		
	}
}
