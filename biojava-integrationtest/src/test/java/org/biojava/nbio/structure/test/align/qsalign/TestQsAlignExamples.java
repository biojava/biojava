package org.biojava.nbio.structure.test.align.qsalign;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.quaternary.QsAlign;
import org.biojava.nbio.structure.align.quaternary.QsAlignParameters;
import org.biojava.nbio.structure.align.quaternary.QsAlignResult;
import org.biojava.nbio.structure.align.quaternary.QsRelation;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.junit.Test;

/**
 * Test the correctness of the {@link QsAlign} algorithm with some examples of
 * different levels of quaternary structure similarity.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class TestQsAlignExamples {

	/**
	 * Identity: test hemoglobin (4HHB) against itself.
	 */
	@Test
	public void testIdentity() throws StructureException, IOException {

		Structure s1 = StructureIO.getStructure("4hhb");
		Structure s2 = s1;

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QsAlignParameters alignParams = new QsAlignParameters();

		QsAlignResult result = QsAlign
				.align(s1, s2, clusterParams, alignParams);

		assertEquals(result.length(), 4);
		assertEquals(result.getRelation(), QsRelation.EQUIVALENT);
		assertEquals(result.getRmsd(), 0.0, 0.01);

	}
	
	/**
	 * Different: test two completely different proteins (4HHB, 3IFV).
	 */
	@Test
	public void testDifferent() throws StructureException, IOException {

		Structure s1 = StructureIO.getStructure("4hhb");
		Structure s2 = StructureIO.getStructure("3ifv");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QsAlignParameters alignParams = new QsAlignParameters();

		QsAlignResult result = QsAlign
				.align(s1, s2, clusterParams, alignParams);

		assertEquals(result.length(), 0);
		assertEquals(result.getRelation(), QsRelation.DIFFERENT);

	}

	/**
	 * Proliferating cell nuclear antigens (3IFV, 3HI8) are structurally
	 * equivalent C3 homotrimers.
	 */
	@Test
	public void testHomoEquivalent() throws StructureException, IOException {

		Structure s1 = StructureIO.getStructure("3ifv");
		Structure s2 = StructureIO.getStructure("BIO:3hi8:1");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QsAlignParameters alignParams = new QsAlignParameters();

		QsAlignResult result = QsAlign
				.align(s1, s2, clusterParams, alignParams);

		assertEquals(result.length(), 3);
		assertEquals(result.getRelation(), QsRelation.EQUIVALENT);
		assertTrue(result.getRmsd() < 10.0);

	}

	/**
	 * Phycocyanins (2VML, 2BV8) are equivalent D3 heterododecamers with A6B6
	 * stoichiometry.
	 */
	@Test
	public void testHeteroEquivalent() throws StructureException, IOException {

		Structure s1 = StructureIO.getStructure("2vml");
		Structure s2 = StructureIO.getStructure("2bv8");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QsAlignParameters alignParams = new QsAlignParameters();

		QsAlignResult result = QsAlign
				.align(s1, s2, clusterParams, alignParams);

		assertEquals(result.length(), 12);
		assertEquals(result.getRelation(), QsRelation.EQUIVALENT);
		assertTrue(result.getRmsd() < 10.0);
	}

	/**
	 * Hydratases (2B3M dimer, 1Q6W hexamer). The C2 dimer is
	 * triplicated into a D3 assembly.
	 */
	@Test
	public void testPartialComplete() throws StructureException,
			IOException {

		Structure s1 = StructureIO.getStructure("BIO:2b3m:1");
		Structure s2 = StructureIO.getStructure("BIO:1Q6W:7");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QsAlignParameters alignParams = new QsAlignParameters();

		QsAlignResult result = QsAlign
				.align(s1, s2, clusterParams, alignParams);

		assertEquals(result.length(), 2);
		assertEquals(result.getRelation(), QsRelation.PARTIAL_COMPLETE);
		assertTrue(result.getRmsd() < 10.0);

	}

	/**
	 * Cytochrome bc1 complexes (1BCC, 1KB9) have some equivalent Chains and
	 * some unmatched.
	 */
	@Test
	public void testPartialIncomplete() throws StructureException,
			IOException {

		Structure s1 = StructureIO.getStructure("1bcc");
		Structure s2 = StructureIO.getStructure("1kb9");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QsAlignParameters alignParams = new QsAlignParameters();

		QsAlignResult result = QsAlign
				.align(s1, s2, clusterParams, alignParams);

		assertEquals(result.length(), 8);
		assertEquals(result.getRelation(), QsRelation.PARTIAL_INCOMPLETE);
		assertTrue(result.getRmsd() < 10.0);
	}

}
