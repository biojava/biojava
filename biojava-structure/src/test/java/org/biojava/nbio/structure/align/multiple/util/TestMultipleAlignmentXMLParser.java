package org.biojava.nbio.structure.align.multiple.util;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import javax.vecmath.Matrix4d;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.TestSampleGenerator;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.xml.MultipleAlignmentXMLParser;
import org.junit.Test;
import org.xml.sax.SAXException;

/**
 * Test the correctness of converting and parsing a MultipleAlignment
 * into an XML format. It checks that the alignments before converting
 * and after parsing are equivalent.
 * <p>
 * To test for the correctness of the XML format see the other test 
 * {@link TestMultipleAlignmentWriter}.
 * 
 * @author Aleix Lafita
 *
 */
public class TestMultipleAlignmentXMLParser {
	
	@Test
	public void testRecovery1() throws StructureException, IOException, 
			ParserConfigurationException, SAXException {
		
		MultipleAlignment sampleMSA = TestSampleGenerator.testAlignment1();
		MultipleAlignmentEnsemble before = sampleMSA.getEnsemble();
		
		String xml = MultipleAlignmentWriter.toXML(before);
		
		MultipleAlignmentEnsemble after = 
				MultipleAlignmentXMLParser.parseXMLfile(xml).get(0);
		
		after.getAtomArrays();
		after.getDistanceMatrix();
		
		assertTrue(equals(before,after));
	}
	
	@Test
	public void testRecovery2() throws StructureException, IOException, 
			ParserConfigurationException, SAXException {
		
		MultipleAlignment sampleMSA = TestSampleGenerator.testAlignment2();
		MultipleAlignmentEnsemble before = sampleMSA.getEnsemble();
		
		String xml = MultipleAlignmentWriter.toXML(before);
		
		MultipleAlignmentEnsemble after = 
				MultipleAlignmentXMLParser.parseXMLfile(xml).get(0);
		
		after.getAtomArrays();
		after.getDistanceMatrix();
		
		assertTrue(equals(before,after));
	}
	
	/**
	 * Compares the basic properties of two MultipleAlignmentEnsembles
	 * to determine if they are equal or not. Recursively checks if
	 * its member MultipleAlignments are equal.<p>
	 * Consider moving this code inside the MultipleAlignment DS.
	 * 
	 * @param a
	 * @param b
	 * @return true if they are equal, false otherwise
	 */
	private static boolean equals(MultipleAlignmentEnsemble a, 
			MultipleAlignmentEnsemble b){
		
		//Check creation properties
		if (a.getAlgorithmName()!=null){
			if (!a.getAlgorithmName().equals(b.getAlgorithmName())) 
				return false;
		} else if (b.getAlgorithmName()!=null) return false;
		
		if (a.getVersion()!=null){
			if (!a.getVersion().equals(b.getVersion())) 
				return false;
		} else if (b.getVersion()!=null) return false;
		
		if (a.getIoTime()!=null){
			if (!a.getIoTime().equals(b.getIoTime()))
				return false;
		} else if (b.getIoTime()!=null) return false;
		
		if (a.getCalculationTime()!=null){
			if (!a.getCalculationTime().equals(b.getCalculationTime()))
				return false;
		} else if (b.getCalculationTime()!=null) return false;
		
		if (a.getStructureNames()!=null){
			if (!a.getStructureNames().equals(b.getStructureNames()))
				return false;
		} else if (b.getStructureNames()!=null) return false;
		
		//Check sizes and lengths
		if (a.size() != b.size()) return false;
		if (a.getMultipleAlignments().size()!=b.getMultipleAlignments().size())
			return false;
		
		//Recursively check member alignments
		for (int i=0; i<a.getMultipleAlignments().size(); i++){
			MultipleAlignment msa1 = a.getMultipleAlignment(i);
			MultipleAlignment msa2 = b.getMultipleAlignment(i);
			if (!equals(msa1,msa2)) return false;
		}
		
		return true;
	}
	
	/**
	 * Compares the basic properties of two MultipleAlignments
	 * to determine if they are equal or not. Recursively checks if
	 * its member BlockSets are equal.<p>
	 * Consider moving this code inside the MultipleAlignment DS.
	 * 
	 * @param a
	 * @param b
	 * @return true if they are equal, false otherwise
	 */
	private static boolean equals(MultipleAlignment a, MultipleAlignment b){
		
		if (a.getBlockSets().size() != b.getBlockSets().size())
			return false;
		
		//Recursively check member alignments
		for (int i=0; i<a.getBlockSets().size(); i++){
			BlockSet bs1 = a.getBlockSets().get(i);
			BlockSet bs2 = b.getBlockSets().get(i);
			if (!equals(bs1,bs2)) return false;
		}
		return true;
	}
	
	/**
	 * Compares the basic properties of two BlockSets
	 * to determine if they are equal or not. Recursively checks if
	 * its member Blocks are equal.<p>
	 * Consider moving this code inside the MultipleAlignment DS.
	 * 
	 * @param a
	 * @param b
	 * @return true if they are equal, false otherwise
	 */
	private static boolean equals(BlockSet a, BlockSet b){
		
		if (a.getBlocks().size() != b.getBlocks().size())
			return false;
		
		//Recursively check member alignments
		for (int i=0; i<a.getBlocks().size(); i++){
			Block b1 = a.getBlocks().get(i);
			Block b2 = b.getBlocks().get(i);
			if (!equals(b1,b2)) return false;
		}
		//Check transformation matrices
		if (a.getTransformations() == null){
			if (b.getTransformations() != null) return false;
		} else if (b.getTransformations() == null) return false;
		else {
			for (int i=0; i<a.getTransformations().size(); i++){
				Matrix4d t1 = a.getTransformations().get(i);
				Matrix4d t2 = b.getTransformations().get(i);
				if (!t1.equals(t2)) return false;
			}
		}
		return true;
	}
	
	/**
	 * Compares the basic properties of two Blocks to determine if 
	 * they are equal or not.
	 * Consider moving this code inside the MultipleAlignment DS.
	 * 
	 * @param a
	 * @param b
	 * @return true if they are equal, false otherwise
	 */
	private static boolean equals(Block a, Block b){
		
		if (a.length() != b.length()) return false;
		if (a.getCoreLength() != b.getCoreLength()) return false;
		
		//Check all aligned residues
		for (int i=0; i<a.getAlignRes().size(); i++){
			List<Integer> chain1 = a.getAlignRes().get(i);
			List<Integer> chain2 = b.getAlignRes().get(i);
			if (!chain1.equals(chain2)) return false;
		}
		return true;
	}
}
