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

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.junit.Test;

/**
 * Created by edlunde-dnastar
 * @since 10/30/2015.
 */
public class TestParseMmCIFLigands {

	private static final int HEM_COUNT_4HHB = 172;	//Number of atoms in HEM groups of 4HHB (manually determined from CIF file)
	private static final int ATOM_COUNT_3UCB = 114;      //number of atoms in 3UCB, including alternate ligand conformations


	@Test
	public void testLigandConnections()throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		// This needs MMCIF
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		StructureIO.setAtomCache(cache);

		Structure sCif = StructureIO.getStructure("4HHB");

		//Verify that we have all HEM atoms from the CIF file.
		assertEquals( HEM_COUNT_4HHB, countBondedAtomsInLigandGroups(sCif) );

	}

	private int countBondedAtomsInLigandGroups(Structure s){

		int count = 0;

		for (Chain c:s.getChains()) {
			for (Group g:c.getAtomGroups()) {
				if (!g.isWater() && !PolymerType.ALL_POLYMER_TYPES.contains(g.getChemComp().getPolymerType())) {

					for (Atom a:g.getAtoms()) {
						if (a.getBonds()!=null) count++;
					}
					for (Group altg : g.getAltLocs()) {
						for (Atom a:altg.getAtoms()) {
							if (a.getBonds()!=null) count++;
						}
					}
				}
			}
		}
		return count;
	}

	@Test
	public void testMultipleConformations()throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		// This needs MMCIF
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		StructureIO.setAtomCache(cache);

		Structure sCif = StructureIO.getStructure("3UCB");


		//Verify that we have all atoms from all conformations of the ligands

		assertEquals(ATOM_COUNT_3UCB, countBondedAtomsInLigandGroups(sCif));
	}

}
