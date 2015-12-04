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
package org.biojava.nbio.structure.symmetry.core;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TestQuatSymmetryDetection {

	private static final Logger logger = LoggerFactory.getLogger(TestQuatSymmetryDetection.class);
	
	@Test
	public void test1b4c() throws Exception {
		// an NMR multi-model entry		
		FileParsingParameters params = new FileParsingParameters();
		params.setParseBioAssembly(true);
		params.setAlignSeqRes(true);		
		AtomCache cache = new AtomCache();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		Structure pdb = StructureIO.getStructure("1b4c");
		
		String[] symmetries = getSymmetry(pdb, 1);
		
		// C2 symmetry
		assertEquals("C2",symmetries[0]);
		assertEquals("A2",symmetries[1]);
		// no pseudosymmetry
		assertNull(symmetries[2]);
		assertNull(symmetries[3]);
	}
	
	@Test
	public void test4hhb() throws Exception {
		// hemoglobin: has both symmetry and pseudosymmetry
		FileParsingParameters params = new FileParsingParameters();
		params.setParseBioAssembly(true);
		params.setAlignSeqRes(true);
		AtomCache cache = new AtomCache();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		Structure pdb = StructureIO.getStructure("4hhb");
		String[] symmetries = getSymmetry(pdb, 1);
		
		// C2 symmetry
		assertEquals("C2",symmetries[0]);
		assertEquals("A2B2",symmetries[1]);
		// D2 pseudo-symmetry
		assertEquals("D2",symmetries[2]);
		assertEquals("A4",symmetries[3]);		
	}
	
	@Test
	public void test4p2c() throws IOException, StructureException {
		// a structure with no global symmetry and some local symmetry
		FileParsingParameters params = new FileParsingParameters();
		params.setParseBioAssembly(true);
		params.setAlignSeqRes(true);
		AtomCache cache = new AtomCache();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		Structure pdb = StructureIO.getStructure("4p2c");
		String[] symmetries = getSymmetry(pdb, 1);
		
		// C1 global symmetry
		assertEquals("C1",symmetries[0]);
		assertEquals("A5B5C",symmetries[1]);
		
		// C3 local symmetry
		assertEquals("C5", symmetries[4]);
		assertEquals("A5B5", symmetries[5]);
	}

	/**
	 * Finds the symmetry of the biounit with the biojava quat symmetry algorithms
	 * @param bioUnitNumber
	 * @return an array of size 6 with members: symmetry, stoichiometry, pseudosymmetry, pseudoStoichiometry, localSymmetry, localStoichiometry
	 */
	private String[] getSymmetry(Structure pdb, int bioUnitNumber) {
		
		List<BiologicalAssemblyTransformation> transformations = pdb.getPDBHeader().getBioAssemblies().get(bioUnitNumber).getTransforms();
		
		if ( transformations == null || transformations.size() == 0){
			logger.warn("Could not load transformations for PDB biounit {}. Will not assign a symmetry value to it.", bioUnitNumber);
			return new String[]{null,null,null,null};
		}
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		Structure bioAssembly = builder.rebuildQuaternaryStructure(pdb, transformations);

		QuatSymmetryParameters parameters = new QuatSymmetryParameters();
        parameters.setOnTheFly(true);
		parameters.setVerbose(false);

		QuatSymmetryDetector detector = new QuatSymmetryDetector(bioAssembly, parameters);

		List<QuatSymmetryResults> globalResults = detector.getGlobalSymmetry(); 
		List<List<QuatSymmetryResults>> localResults = detector.getLocalSymmetries();
		
		if (globalResults.isEmpty()) {
			logger.warn("No global symmetry found for biounit {}. Will not assign a symmetry value to it.",  bioUnitNumber);
			return new String[]{null, null, null, null};
		}
		
		String symmetry = null;
		String stoichiometry = null;
		String pseudoSymmetry = null;
		String pseudoStoichiometry = null;
		String localSymmetry = null;
		String localStoichiometry = null;

		
		if (globalResults.size()>2) {
			StringBuilder sb = new StringBuilder();
			for (QuatSymmetryResults r:globalResults) {
				sb.append(r.getSymmetry()+" ");
			}
			logger.warn("More than 2 symmetry results found for biounit {}. The {} results are: {}", 
					bioUnitNumber, globalResults.size(), sb.toString());
		}
		
		for (QuatSymmetryResults r:globalResults) {
			
			if (r.getSubunits().isPseudoSymmetric()) {				
				pseudoSymmetry = r.getSymmetry();
				pseudoStoichiometry = r.getSubunits().getStoichiometry();
			} else {
				symmetry = r.getSymmetry();
				stoichiometry = r.getSubunits().getStoichiometry();
			}
			
		}
		// note: if there's no pseudosymmetry in the results then it remains null

		for (List<QuatSymmetryResults> r: localResults) {
			
			for (QuatSymmetryResults l:r) {
				localSymmetry = l.getSymmetry();
				localStoichiometry = l.getSubunits().getStoichiometry();
			}
		}

		
		return new String[]{symmetry, stoichiometry, pseudoSymmetry, pseudoStoichiometry, localSymmetry, localStoichiometry};
		
	}
}
