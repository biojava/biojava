package org.biojava.nbio.structure.align.quaternary;

import java.io.IOException;

import org.biojava.nbio.structure.StructureException;
import org.junit.Test;

/**
 * Test the correctness of the {@link QsAlign} algorithm with some examples of
 * different levels of quaternary structure similarity.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class TestQsAlign {

	/**
	 * Identity: test hemoglobin against itself.
	 */
	@Test
	public void testIdentity() throws StructureException, IOException {

	}

	/**
	 * Proliferating cell nuclear antigens (1PLR, 3HI8) are structurally
	 * equivalent C3 homotrimers.
	 */
	@Test
	public void testHomoEquivalent() throws StructureException, IOException {

	}

	/**
	 * Phycocyanins (2VML, 2BV8) are equivalent D3 heterododecamers with A6B6
	 * stoichiometry.
	 */
	@Test
	public void testHeteroEquivalent() throws StructureException, IOException {

	}

	/**
	 * Cytochrome bc1 complexes (1BCC, 1KB9) have some equivalent Chains and
	 * some unmatched.
	 */
	@Test
	public void testHeteroPartialComplete() throws StructureException,
			IOException {

	}

}
