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
package org.biojava.nbio.structure.align.multiple;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.MultipleSuperimposer;
import org.biojava.nbio.structure.align.multiple.util.ReferenceSuperimposer;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.StructureIOFile;

/**
 * This class is a helper for all the other MultipleAlignment test classes.
 * It contains methods for generating sample MultipleAlignments with known
 * properties that can be used to test correctness of the calculation methods.
 * 
 * @author Aleix Lafita
 *
 */
public class TestSampleGenerator {

	/**
	 * Generate a MultipleAlignment of 3 identical structures slightly
	 * missaligned by a repetitive pattern. The resulting alignment contains
	 * 2 BlockSets with 2 Blocks each. Gaps and discontinuities are also
	 * included.<p>
	 * Fields left unfilled (null): distanceMatrices and ioTime.<p>
	 * Atoms are not downladed, but the structure (2gox) is obtained from the
	 * test/resources folder.
	 * 
	 * @return
	 * @throws StructureException
	 * @throws IOException
	 */
	public static MultipleAlignment testAlignment1()
			throws StructureException, IOException {

		//Obtain the structure atoms from resources
		StructureIOFile reader = new PDBFileReader();
		File f = new File("src/test/resources/2gox.pdb");
		Structure structure = null;
		try {
			structure = reader.getStructure(f);
		} catch (IOException e){
			AtomCache cache = new AtomCache();
			structure = cache.getStructure("2gox");
		}
		List<Atom[]> atomArrays = new ArrayList<Atom[]>(3);
		for (int str=0; str<3; str++){
			Atom[] atoms = StructureTools.
					getRepresentativeAtomArray(structure);
			atomArrays.add(StructureTools.cloneAtomArray(atoms));
		}

		//Set the ensemble properties
		MultipleAlignment msa = new MultipleAlignmentImpl();
		MultipleAlignmentEnsemble ensemble = msa.getEnsemble();
		ensemble.setAtomArrays(atomArrays);
		ensemble.setAlgorithmName("testAlignment");
		ensemble.setVersion("1.0");
		ensemble.setCalculationTime((long) 1000000000);
		ensemble.setStructureNames(Arrays.asList("2gox","2gox","2gox"));

		//Generate the MultipleAlignment - 2 blocks with 2 blocksets each
		int[] nextResidue = new int[3];
		for (int bs=0; bs<2; bs++){
			BlockSet blockSet = new BlockSetImpl(msa);
			for (int b=0; b<2; b++){
				List<List<Integer>> alnRes = new ArrayList<List<Integer>>(3);
				for (int str=0; str<3; str++){
					List<Integer> chain = new ArrayList<Integer>(50);
					for (int res=0; res<10; res++){
						//Introduce gaps and discontinuities
						if (nextResidue[str] % (2+str) == str)chain.add(null);
						else chain.add(nextResidue[str]);
						if (nextResidue[str] % (10) == str) nextResidue[str]++;
						nextResidue[str]++;
					}
					alnRes.add(chain);
					nextResidue[str]+=str;  //Spacing between Blocks
				}
				Block block = new BlockImpl(blockSet);
				block.setAlignRes(alnRes);
			}
		}

		//Superposition and scores
		MultipleSuperimposer imposer = new ReferenceSuperimposer();
		imposer.superimpose(msa);
		MultipleAlignmentScorer.calculateScores(msa);

		return msa;
	}

	/**
	 * Method that builds a MultipleAlignment of 4 structures with 2 BlockSets,
	 * with 2 and 1 Blocks respectively. The alignment contains gaps
	 * and non-consecutive residues, ideal to test for all possible cases.
	 * The alignment is manually constructed. All ensemble fields are filled,
	 * so that no one is null.<p>
	 * Atoms of four structures (globins) are downloaded.
	 * 
	 * @return MultipleAlignment the test multiple alignment
	 * @throws StructureException
	 * @throws IOException
	 */
	public static MultipleAlignment testAlignment2()
			throws StructureException, IOException {

		//Download the globin structures
		List<String> names = Arrays.asList("1mbc", "1hlb", "1thb.A", "1ith.A");
		AtomCache cache = new AtomCache();
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		for (String name:names) atomArrays.add(cache.getAtoms(name));

		MultipleAlignmentEnsemble ensemble = 
				new MultipleAlignmentEnsembleImpl();
		MultipleAlignment alignment = new MultipleAlignmentImpl(ensemble);

		//Set the ensemble properties (all filled)
		ensemble.setAtomArrays(atomArrays);
		ensemble.getDistanceMatrix();
		ensemble.setStructureNames(names);
		ensemble.setAlgorithmName("testAlignment");
		ensemble.setVersion("2.0");
		ensemble.setIoTime((long) 1000000000);
		ensemble.setCalculationTime((long) 1000000000);

		//Build the aligned positions: 2 BlockSets, 3 Blocks
		BlockSet blockSet1 = new BlockSetImpl(alignment);
		BlockSet blockSet2 = new BlockSetImpl(alignment);

		Block block1 = new BlockImpl(blockSet1);
		Block block2 = new BlockImpl(blockSet1);
		Block block3 = new BlockImpl(blockSet2);
		block1.setAlignRes(new ArrayList<List<Integer>>());
		block2.setAlignRes(new ArrayList<List<Integer>>());
		block3.setAlignRes(new ArrayList<List<Integer>>());

		List<Integer> aligned11 = Arrays.asList(
				0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21);
		List<Integer> aligned12 = Arrays.asList(
				29,30,31,32,33,34,35,36,38);
		List<Integer> aligned13 = Arrays.asList(
				123,124,125,126,127,128,129,130,131,132,133,134);

		List<Integer> aligned21 = Arrays.asList(10,11,12,13,null,15,16,17,
				null,19,20,21,22,23,24,25,null,27,28,29,30,null);
		List<Integer> aligned22 = Arrays.asList(39,40,41,42,43,44,45,46,48);
		List<Integer> aligned23 = Arrays.asList(133,134,135,136,137,138,139,
				140,141,142,143,144);

		List<Integer> aligned31 = Arrays.asList(0,1,2,3,null,5,6,7,8,9,10,11,
				12,13,14,15,16,17,18,19,20,21);
		List<Integer> aligned32 = Arrays.asList(29,30,31,32,33,34,35,36,38);
		List<Integer> aligned33 = Arrays.asList(117,118,119,120,121,122,123,
				124,125,126,127,128);

		List<Integer> aligned41 = Arrays.asList(0,1,2,3,null,5,6,7,8,9,10,11,
				12,13,14,15,null,17,18,19,20,21);
		List<Integer> aligned42 = Arrays.asList(30,31,32,33,34,35,36,37,39);
		List<Integer> aligned43 = Arrays.asList(121,122,123,124,125,126,127,
				128,129,130,131,132);

		block1.getAlignRes().add(aligned11);
		block1.getAlignRes().add(aligned21);
		block1.getAlignRes().add(aligned31);
		block1.getAlignRes().add(aligned41);

		block2.getAlignRes().add(aligned12);
		block2.getAlignRes().add(aligned22);
		block2.getAlignRes().add(aligned32);
		block2.getAlignRes().add(aligned42);

		block3.getAlignRes().add(aligned13);
		block3.getAlignRes().add(aligned23);
		block3.getAlignRes().add(aligned33);
		block3.getAlignRes().add(aligned43);

		//Superposition and scores
		MultipleSuperimposer imposer = new ReferenceSuperimposer();
		imposer.superimpose(alignment);
		MultipleAlignmentScorer.calculateScores(alignment);
		alignment.length();
		alignment.getCoreLength();

		return alignment;
	}
}
