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
 * Created on Jul 26, 2007
 *
 */

package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ReducedChemCompProvider;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import org.junit.Test;

import static org.junit.Assert.*;


public class PdbFileFormat30Test {




	/**
	 * There is a file format change in v3.0 of the PDB file format
	 * this test makes sure that the atom name changes are being processed correctly
	 * @throws IOException
	 */
	@Test
	public void testRead30File() throws IOException{
		Structure s = getStructure("/388d_v30.pdb");
		int nrNuc = getNrNucleotides(s);

		// there are 4 nucleotides less in the new version
		// some chemically modified nucleotides residues have been declared to be HETATOMS

		int shouldNr = 20;
		assertEquals("structure does not contain the right number of nucleotides ", shouldNr ,nrNuc);

		List<EntityInfo> compounds= s.getEntityInfos();
		// from Biojava 5.0 on we are creating entities whenever an entity is found to be without an assigned compound 
		// in the file, for polymer entities, nonpolymer entities and water entities.
		// For this file: 1 dna polymeric entity, 1 MG nonpolymeric entity, 1 water
		// see issues https://github.com/biojava/biojava/issues/305 and https://github.com/biojava/biojava/pull/394
		assertEquals(3, compounds.size());
		EntityInfo mol = compounds.get(0);
		assertTrue(mol.getDescription().startsWith("DNA"));


		Structure s2 = getStructure("/104D_v30.pdb");

		int nrNuc2 = getNrNucleotides(s2);
		int shouldNr2 = 24;
		assertEquals("structure does not contain the right number of nucleotides ", shouldNr2 , nrNuc2);


	}

	@Test
	public void testRead23File() throws IOException{

		Structure s = getStructure("/388d_v23.pdb");
		int nrNuc = getNrNucleotides(s);
		int shouldNr = 24;
		assertEquals("structure does not contain the right number of nucleotides ", shouldNr , nrNuc);

		List<EntityInfo> compounds= s.getEntityInfos();
		// from Biojava 5.0 on we are creating entities whenever an entity is found to be without an assigned compound 
		// in the file, for polymer entities, nonpolymer entities and water entities. 
		// For this entry: we have 1 dna polymeric entity, 1 FLO nonpoly entity, 1 MO6 nonpoly entity, 1 water entity
		// see issues https://github.com/biojava/biojava/issues/305 and https://github.com/biojava/biojava/pull/394
		assertEquals(4, compounds.size());
		EntityInfo mol = compounds.get(0);

		assertTrue(mol.getDescription().startsWith("DNA"));


		Structure s2 = getStructure("/104D_v23.pdb");

		int nrNuc2 = getNrNucleotides(s2);
		int shouldNr2 = 24;
		assertEquals("structure does not contain the right number of nucleotides ", shouldNr2 , nrNuc2);

	}

	private Structure getStructure(String fileName) throws IOException{

		InputStream inStream = this.getClass().getResourceAsStream(fileName);
		assertNotNull(inStream);

		ChemCompGroupFactory.setChemCompProvider(new ReducedChemCompProvider());
		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(false);
		pdbpars.setFileParsingParameters(params);
		Structure structure = null;

		structure = pdbpars.parsePDBFile(inStream) ;

		return structure;
	}

	private int getNrNucleotides(Structure s){
		GroupIterator iter = new GroupIterator(s);
		int nr = 0;
		while(iter.hasNext()){
			Group g = iter.next();

			if (g.getType().equals(GroupType.NUCLEOTIDE)){
				nr ++;
			} else {
				//System.out.println(g.getType() + g.getPDBName());
			}

		}
		return nr;
	}

	/**
	 * Checks that the legacy file check is working and that that non-legacy
	 * files have the correct number of chains when the line length is over
	 * 72 characters.
	 * @throws IOException
	 */
	@Test
	public void testIsLegacyFormat_pdb_COMPND_handler() throws IOException{

		Structure s = getStructure("/3mk3.pdb");

		List<EntityInfo> compounds= s.getEntityInfos();
		// we are doing heuristics to find missing compounds not specified in header
		// thus here we have the one specified in header plus a SO4 nonpolymer entity
		assertEquals(2, compounds.size());
		EntityInfo mol = compounds.get(0);
		assertTrue(mol.getDescription().equals("6,7-DIMETHYL-8-RIBITYLLUMAZINE SYNTHASE"));
		assertEquals(60, mol.getChainIds().size());
		assertEquals(60, mol.getChains().size());
		assertTrue(isChainNameInEntity(mol,"S"));
		assertTrue(isChainNameInEntity(mol,"T"));
		assertTrue(isChainNameInEntity(mol,"U"));
		assertTrue(isChainNameInEntity(mol,"g"));
		assertTrue(isChainNameInEntity(mol,"h"));
		assertTrue(isChainNameInEntity(mol,"i"));
	}
	
	private boolean isChainNameInEntity(EntityInfo e, String chainName) {
		for (Chain c:e.getChains()) {
			if (c.getName().equals(chainName)) return true;
		}
		return false;
	}
}
