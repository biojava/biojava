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
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.StructureIOFile;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the correctness of various Text outputs for {@link MultipleAlignment}s.
 * <p>
 * Currently tested:
 * <ul><li>FASTA
 * <li>FatCat format
 * <li>Aligned Residues
 * </ul>
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentWriterTest {
	
	@Test
	public void testFASTA() throws StructureException, IOException{
		
		MultipleAlignment alignment = generateTestMultipleAlignment();
		String result = MultipleAlignmentWriter.toFASTA(alignment);
		System.out.println(result);
		
		StringBuffer expected = new StringBuffer();
		expected.append(">2gox\n");
		expected.append("---S-tDaErLkhl--IvTpSgAgeq----NmIgMtPt-viAv---HyL-dEt-eqWe\n");
		expected.append(">2gox\n");
		expected.append("GsrS-tDAeRLkh--LiVTpSGaGEqn---MiGMtPTviA-vh--YlDE-tEqwE-Kf\n");
		expected.append(">2gox\n");
		expected.append("GS-rsTDaERLkhl-IvTPSgAGEqnmig--MTPtVIavH-Yld-ETEqwEKf-G-LE\n");
		
		assertEquals(result,expected.toString());
		
	}
	
	@Test
	public void testFatCat(){
		
		
	}
	
	@Test
	public void testAlignedResidues(){
		
	}
	
	private MultipleAlignment generateTestMultipleAlignment() throws StructureException, IOException{
		
		//Obtain the structure atoms
		StructureIOFile reader = new PDBFileReader();
		File f = new File("src/main/resources/2gox.pdb");
		Structure structure = null;
		try {
			structure = reader.getStructure(f);
		} catch (IOException e){
			AtomCache cache = new AtomCache();
			structure = cache.getStructure("2gox");
		}
		List<Atom[]> atomArrays = new ArrayList<Atom[]>(3);
		for (int str=0; str<3; str++){
			Atom[] atoms = StructureTools.getRepresentativeAtomArray(structure);
			atomArrays.add(StructureTools.cloneAtomArray(atoms));
		}
		
		//Generate the MultipleAlignment - 2 blocks with 2 blocksets each
		MultipleAlignment msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setStructureNames(Arrays.asList("2gox","2gox","2gox"));
		msa.getEnsemble().setAtomArrays(atomArrays);
		int[] nextResidue = new int[3];
		for (int bs=0; bs<2; bs++){
			BlockSet blockSet = new BlockSetImpl(msa);
				for (int b=0; b<2; b++){
				List<List<Integer>> alnRes = new ArrayList<List<Integer>>(3);
				for (int str=0; str<3; str++){
					List<Integer> chain = new ArrayList<Integer>(50);
					for (int res=0; res<10; res++){
						//Introduce gaps and discontinuities to test for all cases
						if (nextResidue[str] % (2+str) == str)chain.add(null);
						else chain.add(nextResidue[str]);
						if (nextResidue[str] % (10) == str) nextResidue[str] ++;
						nextResidue[str]++;
					}
					alnRes.add(chain);
					nextResidue[str]+=str;  //Spacing between Blocks
				}
				Block block = new BlockImpl(blockSet);
				block.setAlignRes(alnRes);
			}
		}
		return msa;
	}
}
