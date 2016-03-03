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

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		StructureIO.setAtomCache(cache);

		Structure sCif = StructureIO.getStructure("4HHB");

		//Verify that we have all HEM atoms from the CIF file.
		assertEquals( HEM_COUNT_4HHB, countUniqueAtomsInLigandGroups(sCif) );
		
	}
	
	private int countUniqueAtomsInLigandGroups(Structure s){

		int count = 0;
		
		for (Chain c:s.getChains()) {
			for (Group g:c.getAtomGroups()) {
				if (!g.isWater() && !PolymerType.ALL_POLYMER_TYPES.contains(g.getChemComp().getPolymerType())) {
					System.out.println(g); 
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
		
		StructureIO.setAtomCache(cache);
		
		cache.setUseMmCif(true);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		 
		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		StructureIO.setAtomCache(cache);
		
		Structure sCif = StructureIO.getStructure("3UCB");

		
		//Verify that we have all atoms from all conformations of the ligands
		
		assertEquals(ATOM_COUNT_3UCB, countUniqueAtomsInLigandGroups(sCif));
	}

}