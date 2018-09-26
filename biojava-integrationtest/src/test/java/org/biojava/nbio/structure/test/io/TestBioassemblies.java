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
package org.biojava.nbio.structure.test.io;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestBioassemblies {

	
	
	/**
	 * A test for an NMR structure
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void test1E17() throws IOException, StructureException {

		AtomCache prevAtomCache = StructureIO.getAtomCache();

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true); 
		StructureIO.setAtomCache(cache);		

		List<Structure> multiModelBioAssemblies = StructureIO.getBiologicalAssemblies("1E17", true);
		
		List<Structure> flattenedBioAssemblies = StructureIO.getBiologicalAssemblies("1E17", false);
		
		// 1 bioassembly in this case
		assertEquals(1, multiModelBioAssemblies.size());
		assertEquals(1, flattenedBioAssemblies.size());
		
		// checking that we have 1 model only (the bioassemblies creation wipes out all models)
		assertEquals(1, multiModelBioAssemblies.get(0).nrModels());		
		assertEquals(1, flattenedBioAssemblies.get(0).nrModels());
		
		// the bioassembly is a monomer
		assertEquals(1, multiModelBioAssemblies.get(0).getPolyChains().size());
		assertEquals(1, flattenedBioAssemblies.get(0).getPolyChains().size());

		StructureIO.setAtomCache(prevAtomCache);
	}
	
	/**
	 * A test for an entry with a biounit that is a subset of the AU
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void test4TTX() throws IOException, StructureException {

		AtomCache prevAtomCache = StructureIO.getAtomCache();
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true); 
		StructureIO.setAtomCache(cache);		

		List<Structure> multiModelBioAssemblies = StructureIO.getBiologicalAssemblies("4TTX", true);
		
		List<Structure> flattenedBioAssemblies = StructureIO.getBiologicalAssemblies("4TTX", false);
		
		// 3 bioassemblies in this case
		assertEquals(3, multiModelBioAssemblies.size());
		assertEquals(3, flattenedBioAssemblies.size());
		
		// checking that we have 1 model only
		assertEquals(1, multiModelBioAssemblies.get(0).nrModels());		
		assertEquals(1, flattenedBioAssemblies.get(0).nrModels());
		
		// the 3 bioassemblies are dimers
		assertEquals(2, multiModelBioAssemblies.get(0).getPolyChains().size());
		assertEquals(2, flattenedBioAssemblies.get(0).getPolyChains().size());
		
		assertEquals(2, multiModelBioAssemblies.get(1).getPolyChains().size());
		assertEquals(2, flattenedBioAssemblies.get(1).getPolyChains().size());
		
		assertEquals(2, multiModelBioAssemblies.get(2).getPolyChains().size());
		assertEquals(2, flattenedBioAssemblies.get(2).getPolyChains().size());

		// all 3 flattened bioassemblies have operator id 1 in their new chain ids
		assertEquals("1", flattenedBioAssemblies.get(0).getPolyChains().get(0).getId().split("_")[1]);
		assertEquals("1", flattenedBioAssemblies.get(1).getPolyChains().get(0).getId().split("_")[1]);
		assertEquals("1", flattenedBioAssemblies.get(2).getPolyChains().get(0).getId().split("_")[1]);

		StructureIO.setAtomCache(prevAtomCache);

	}


	/**
	 * A test for an entry with cartesian product in assembly operators
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void test1M4X() throws IOException, StructureException {

		AtomCache prevAtomCache = StructureIO.getAtomCache();
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);

		Structure flattenedBioAssembly5 = StructureIO.getBiologicalAssembly("1M4X" , 5);

		// checking that we have 1 model only
		assertEquals(1, flattenedBioAssembly5.nrModels());

		// bioassembly 5 expands to 90 chains (3 chains x 5 operators x 6 operators), the expression is '(1-5)(61-88)'
		assertEquals(90, flattenedBioAssembly5.getPolyChains().size());

		// the operator ids are composed for this case, e.g. A_5x61
		assertTrue(flattenedBioAssembly5.getPolyChains().get(0).getId().contains("x"));
		assertEquals(2, flattenedBioAssembly5.getPolyChains().get(0).getId().split("_").length);

		StructureIO.setAtomCache(prevAtomCache);

	}
	
	/**
	 * A difficult case: see http://www.mail-archive.com/jmol-users@lists.sourceforge.net/msg25927.html
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test4OPJ() throws IOException, StructureException {

		AtomCache prevAtomCache = StructureIO.getAtomCache();
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true); 
		StructureIO.setAtomCache(cache);		

		List<Structure> multiModelBioAssemblies = StructureIO.getBiologicalAssemblies("4OPJ", true);
		
		List<Structure> flattenedBioAssemblies = StructureIO.getBiologicalAssemblies("4OPJ", false);
		
		// 2 bioassemblies in this case
		assertEquals(2, multiModelBioAssemblies.size());
		assertEquals(2, flattenedBioAssemblies.size());
		
		// checking number of models: 2 operators in each assembly
		assertEquals(2, multiModelBioAssemblies.get(0).nrModels());		
		assertEquals(1, flattenedBioAssemblies.get(0).nrModels());
		
		assertEquals(2, multiModelBioAssemblies.get(1).nrModels());		
		assertEquals(1, flattenedBioAssemblies.get(1).nrModels());

		// for multimodel bioassembly, we should have 2 models corresponding to 2 operators
		assertEquals(2, multiModelBioAssemblies.get(0).nrModels());
		// for flattened bioassembly we should have only 1 model
		assertEquals(1, flattenedBioAssemblies.get(0).nrModels());
		
		// 3 chains divided into 2 models in bioassembly 1
		assertEquals(3, multiModelBioAssemblies.get(0).getPolyChains(0).size() + multiModelBioAssemblies.get(0).getPolyChains(1).size());
		// 3 chains in flattened structure in bioassembly 1
		assertEquals(3, flattenedBioAssemblies.get(0).getPolyChains().size());

		// 3 chains divided into 2 models in bioassembly 2
		assertEquals(3, multiModelBioAssemblies.get(1).getPolyChains(0).size() + multiModelBioAssemblies.get(1).getPolyChains(1).size());
		// 3 chains in flattened structure in bioassembly 2
		assertEquals(3, flattenedBioAssemblies.get(1).getPolyChains().size());

		// chain ids and names don't contain underscores in multimodel
		for (int modelIdx = 0; modelIdx<multiModelBioAssemblies.get(0).nrModels(); modelIdx++) {
			List<Chain> model = multiModelBioAssemblies.get(0).getModel(modelIdx);
			for (Chain c:model) {
				System.out.println(c.getId()+" "+c.getName());
				assertTrue(!c.getId().contains("_"));
				assertTrue(!c.getName().contains("_"));
			}

		}

		// chain ids and names contain underscores in flattened
		for (Chain c:flattenedBioAssemblies.get(0).getChains()) {
			System.out.println(c.getId()+" "+c.getName());
			assertTrue(c.getId().contains("_"));
			assertTrue(c.getName().contains("_"));
		}

		StructureIO.setAtomCache(prevAtomCache);

	}

}
