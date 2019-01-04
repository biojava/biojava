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
package org.biojava.nbio.structure.align.multiple.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.TestSampleGenerator;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the correctness of various Text outputs for {@link MultipleAlignment}s.
 * <p>
 * Currently tested:
 * <ul>
 * <li>FASTA
 * <li>FatCat format
 * <li>Aligned Residues
 * </ul>
 *
 * @author Aleix Lafita
 *
 */
public class TestMultipleAlignmentWriter {

	private MultipleAlignment alignment1;
	private MultipleAlignment alignment2;

	/**
	 * Build the alignments in common for every writer output, so that they do
	 * not have to be created each time.
	 *
	 * @throws IOException
	 * @throws StructureException
	 */
	public TestMultipleAlignmentWriter() throws StructureException, IOException {

		alignment1 = TestSampleGenerator.testAlignment1();
		alignment2 = TestSampleGenerator.testAlignment2();
	}

	@Test
	public void testFASTA1() throws IOException {

		String result = MultipleAlignmentWriter.toFASTA(alignment1);

		FileReader file = new FileReader("src/test/resources/testMSTA1.fasta");
		BufferedReader reader = new BufferedReader(file);
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append("\n");
		}
		reader.close();

		String expected = stringBuilder.toString();
		assertEquals(expected, result);
	}

	@Test
	public void testFASTA2() throws IOException {

		String result = MultipleAlignmentWriter.toFASTA(alignment2);

		FileReader file = new FileReader("src/test/resources/testMSTA2.fasta");
		BufferedReader reader = new BufferedReader(file);
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append("\n");
		}
		reader.close();

		String expected = stringBuilder.toString();
		assertEquals(expected, result);
	}

	@Test
	public void testFatCat1() throws IOException {

		String result = MultipleAlignmentWriter.toFatCat(alignment1);

		FileReader file = new FileReader("src/test/resources/testMSTA1.fatcat");
		BufferedReader reader = new BufferedReader(file);
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append("\n");
		}
		reader.close();

		String expected = stringBuilder.toString();
		assertEquals(expected, result);
	}

	@Test
	public void testFatCat2() throws IOException {

		String result = MultipleAlignmentWriter.toFatCat(alignment2);

		FileReader file = new FileReader("src/test/resources/testMSTA2.fatcat");
		BufferedReader reader = new BufferedReader(file);
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append("\n");
		}
		reader.close();

		String expected = stringBuilder.toString();
		assertEquals(expected, result);
	}

	@Test
	public void testAlignedResidues1() throws IOException {

		String result = MultipleAlignmentWriter.toAlignedResidues(alignment1);

		FileReader file = new FileReader(
				"src/test/resources/testMSTA1_alnres.tsv");
		BufferedReader reader = new BufferedReader(file);
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append("\n");
		}
		reader.close();

		String expected = stringBuilder.toString();
		assertEquals(expected, result);
	}

	@Test
	public void testAlignedResidues2() throws IOException {

		String result = MultipleAlignmentWriter.toAlignedResidues(alignment2);

		FileReader file = new FileReader(
				"src/test/resources/testMSTA2_alnres.tsv");
		BufferedReader reader = new BufferedReader(file);
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append("\n");
		}
		reader.close();

		String expected = stringBuilder.toString();
		assertEquals(expected, result);
	}

}
