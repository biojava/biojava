package org.biojava.nbio.structure.align;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class TestAlignmentConsistency {

	// Check that indices of the aligned residues are unique
	@Test
	public void testDuplicateIndices() throws IOException, StructureException {
		String[] algorithmIDs = {CeMain.algorithmName, FatCatRigid.algorithmName};

		AtomCache cache = new AtomCache();

		// 3j47 is a bunch of a-helices, so there are many valid ways to align chains
		// structurally between each other.
		List<Chain> chains = cache.getStructure("3j47").getChains();

		for(String algorithmID:algorithmIDs) {
			StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(algorithmID);
			System.out.println("Testing "+algorithmID);
			for (int c1 = 0; c1<chains.size()-1;c1++) {
				for (int c2=chains.size()-1;c2>c1;c2--) {
					Atom[] ca1 = StructureTools.getAtomCAArray(chains.get(c1));
					Atom[] ca2 = StructureTools.getAtomCAArray(chains.get(c2));

					AFPChain afpChain_fc = algorithm.align(ca1, ca2);
					assertEquals(1,afpChain_fc.getOptAln().length);

					int[][] optAln = afpChain_fc.getOptAln()[0];
					// two chains aligned
					assertEquals(2,optAln.length);

					//same number of aligned residues between the chains
					assertEquals(optAln[0].length,optAln[1].length);

					// no indices duplicated in the alignments
					for (int[] optAlnSeq : optAln) {
						long count_unique = Arrays.stream(optAlnSeq).distinct().count();
						long count_all = optAlnSeq.length;
						assertEquals(count_unique, count_all);
					}

				}
			}

		}

	}

}
