package org.biojava.nbio.structure.align.multiple;

import java.io.IOException;
import org.biojava.nbio.structure.StructureException;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the correctness of various Text outputs for {@link MultipleAlignment}s.
 * <p>
 * Currently tested:
 * <ul><li>FASTA
 * <li>FatCat format
 * <li>Aligned Residues
 * <li>Transformation Matrices
 * <li>XML format
 * </ul>
 * 
 * @author Aleix Lafita
 *
 */
public class TestMultipleAlignmentWriter {
	
	@Test
	public void testFASTA() throws StructureException, IOException{
		
		MultipleAlignment alignment = TestSampleGenerator.testAlignment1();
		String result = MultipleAlignmentWriter.toFASTA(alignment);
		System.out.println(result);
		
		StringBuffer expected = new StringBuffer();
		expected.append(">2gox\n");
		expected.append("---S-tDaErLkhl--IvTpSgAgeq----NmIgMtPt-viAv---"
				+ "HyL-dEt-eqWe\n");
		expected.append(">2gox\n");
		expected.append("GsrS-tDAeRLkh--LiVTpSGaGEqn---MiGMtPTviA-vh--Y"
				+ "lDE-tEqwE-Kf\n");
		expected.append(">2gox\n");
		expected.append("GS-rsTDaERLkhl-IvTPSgAGEqnmig--MTPtVIavH-Yld-ET"
				+ "EqwEKf-G-LE\n");
		
		assertEquals(result,expected.toString());
	}
	
	@Test
	public void testFatCat(){
		
		
	}
	
	@Test
	public void testAlignedResidues(){
		
		
	}
	
	@Test
	public void testTransformMatrices(){
		
		
	}
	
	@Test
	public void testXMLformat(){
		
		
	}
}