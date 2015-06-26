package org.biojava.nbio.structure.align.multiple;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Test that all relevant information (except scores and cache) is correctly copied
 * from the AFPChain to the generated MultipleAlignment object.
 * 
 * @author Aleix Lafita
 *
 */
public class AFPChainConversionTest {

	@Test
	public void testAFPconversion() throws Exception{
		
		//Fill an AFPChain with the general information
		AFPChain afp = new AFPChain();
		afp.setName1("name1");
		afp.setName2("name2");
		afp.setAlgorithmName("algorithm");
		afp.setVersion("1.0");
		afp.setCalculationTime(System.currentTimeMillis());
		//Generate a fake optimal alignment with three blocks and 5 residues per block
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
		//Set the rotation matrix to the identity and the shift vector to the origin
		Matrix rot = Matrix.identity(3, 3);
		Atom shift = new AtomImpl();
		shift.setX(0);
		shift.setY(0);
		shift.setZ(0);
		Matrix[] blockRot = {rot,rot,rot};
		afp.setBlockRotationMatrix(blockRot);
		Atom[] blockShift = {shift,shift,shift};
		afp.setBlockShiftVector(blockShift);
		
		//Convert the AFPChain into a MultipleAlignment (without Atoms)
		MultipleAlignmentEnsemble ensemble = new MultipleAlignmentEnsembleImpl(afp,null,null);
		MultipleAlignment msa = ensemble.getMultipleAlignments().get(0);
		
		//Test for all the information to be equal
		assertEquals(afp.getName1(),ensemble.getStructureNames().get(0));
		assertEquals(afp.getName2(), ensemble.getStructureNames().get(1));
		assertEquals(afp.getAlgorithmName(), ensemble.getAlgorithmName());
		assertEquals(afp.getVersion(),ensemble.getVersion());
		assertTrue(ensemble.getCalculationTime().equals(afp.getCalculationTime()));
		assertEquals(afp.getBlockNum(), msa.getBlockSets().size());
		assertEquals(Calc.getTransformation(afp.getBlockRotationMatrix()[0], afp.getBlockShiftVector()[0]), msa.getTransformations().get(1));
		
		//Test for the optimal alignment to be equal
		for (int b=0; b<3; b++){
			for (int c=0; c<2; c++){
				for (int res=0; res<5; res++){
					assert(afp.getOptAln()[b][c][res] == msa.getBlocks().get(b).getAlignRes().get(c).get(res));
				}
			}
		}
	}
		
}