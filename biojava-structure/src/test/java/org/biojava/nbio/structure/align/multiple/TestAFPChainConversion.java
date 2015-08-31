package org.biojava.nbio.structure.align.multiple;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.jama.Matrix;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test that all relevant information (except scores and cache) is correctly 
 * copied from the AFPChain to the generated MultipleAlignment object.
 * 
 * @author Aleix Lafita
 *
 */
public class TestAFPChainConversion {

	@Test
	public void testAFPconversion() throws Exception{

		//Fill an AFPChain with the general information
		AFPChain afp = new AFPChain();
		afp.setName1("name1");
		afp.setName2("name2");
		afp.setAlgorithmName("algorithm");
		afp.setVersion("1.0");
		afp.setCalculationTime(System.currentTimeMillis());
		//Generate a optimal alignment with 3 blocks and 5 residues per block
		int[][][] optAln = new int[3][][];
		for (int b=0; b<optAln.length; b++){
			int[][] block = new int[2][];
			for (int c=0; c<block.length; c++){
				int[] residues = {b+5,b+6,b+7,b+8,b+9};
				block[c] = residues;
			}
			optAln[b] = block;
		}
		afp.setOptAln(optAln);
		afp.setBlockNum(optAln.length);
		//Set the rotation matrix and shift to random numbers
		double[][] mat = {{0.13,1.5,0.84},{1.3,0.44,2.3},{1.0,1.2,2.03}};
		Matrix rot = new Matrix(mat);
		Atom shift = new AtomImpl();
		shift.setX(0.44);
		shift.setY(0.21);
		shift.setZ(0.89);
		Matrix[] blockRot = {rot,rot,rot};
		afp.setBlockRotationMatrix(blockRot);
		Atom[] blockShift = {shift,shift,shift};
		afp.setBlockShiftVector(blockShift);

		//Convert the AFPChain into a MultipleAlignment (without Atoms)
		MultipleAlignmentEnsemble ensemble = 
				new MultipleAlignmentEnsembleImpl(afp,null,null,true);
		MultipleAlignment msa = ensemble.getMultipleAlignment(0);

		//Test for all the information
		assertEquals(afp.getName1(),ensemble.getStructureNames().get(0));
		assertEquals(afp.getName2(), ensemble.getStructureNames().get(1));
		assertEquals(afp.getAlgorithmName(), ensemble.getAlgorithmName());
		assertEquals(afp.getVersion(),ensemble.getVersion());
		assertTrue(ensemble.getCalculationTime().equals(
				afp.getCalculationTime()));
		assertEquals(afp.getBlockNum(), msa.getBlockSets().size());
		for (int b = 0; b<afp.getBlockNum(); b++){
			assertEquals(Calc.getTransformation(
					afp.getBlockRotationMatrix()[b],
					afp.getBlockShiftVector()[b]),
					msa.getBlockSet(b).getTransformations().get(1));
		}

		//Test for the scores
		assertEquals(msa.getScore(MultipleAlignmentScorer.CE_SCORE),
				(Double) afp.getAlignScore());
		assertEquals(msa.getScore(MultipleAlignmentScorer.AVGTM_SCORE),
				(Double) afp.getTMScore());
		assertEquals(msa.getScore(MultipleAlignmentScorer.RMSD),
				(Double) afp.getTotalRmsdOpt());


		//Test for the optimal alignment
		for (int b=0; b<3; b++){
			for (int c=0; c<2; c++){
				for (int res=0; res<5; res++){
					Integer afpRes = afp.getOptAln()[b][c][res];
					assertEquals(afpRes, msa.getBlock(b).
							getAlignRes().get(c).get(res));
				}
			}
		}
	}
}