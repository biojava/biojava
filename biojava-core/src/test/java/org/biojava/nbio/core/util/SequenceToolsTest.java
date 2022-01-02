package org.biojava.nbio.core.util;

import static org.junit.Assert.assertThrows;
import static org.junit.jupiter.api.Assertions.assertAll;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.EmptySource;
import org.junit.jupiter.params.provider.NullAndEmptySource;
import org.junit.jupiter.params.provider.NullSource;

class SequenceToolsTest {

	String randomDNA(int n) {
		String[] nucs = new String[] { "A", "T", "C", "G" };
		Random r = new Random();
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < n; i++) {
			sb.append(nucs[r.nextInt(4)]);
		}
		return sb.toString();
	}

	@Nested
	class PermuteCyclic {

		@ParameterizedTest
		@CsvSource(value = { "ATCGT,1,TCGTA", "ATCGT,-1,TATCG", "ATCGT,0,ATCGT", "ATCGT,25,ATCGT","12345,1,23451" })
		void permuteCyclicBasic(String original, int n, String expected) {
			assertEquals(expected, SequenceTools.permuteCyclic(original, n));
		}
		
		@ParameterizedTest
		@CsvSource(value = { "ATCGT,CGTAT", "ATCGT,CGTAT" })
		@Disabled("fails with current implementation")
		void permuteCycleIntMaxMin(String original, String expected) {
			assertAll(
					()->assertEquals(expected, SequenceTools.permuteCyclic(original, Integer.MAX_VALUE)),
					()->assertEquals(expected, SequenceTools.permuteCyclic(original, Integer.MIN_VALUE))
				);
		}
		
		@ParameterizedTest
		@CsvSource(value = { "ATCGT,CGTAT", "ATCGT,CGTAT" })
		@DisplayName("Edge case fixed")
		void permuteCycleIntMaxMin2(String original, String expected) {
			assertAll(
					()->assertEquals(expected, SequenceTools.permuteCyclic(original, Integer.MAX_VALUE)),
					()->assertEquals(expected, SequenceTools.permuteCyclic(original, Integer.MIN_VALUE))
				);
		}

	}

	@Nested
	class PercentNucleotideContent {

		@ParameterizedTest
		@NullAndEmptySource
		@DisplayName("percent nucleotide sequence returns 0 for null "+
		 "or empty string")
		void nucleotideContentInvalidValues(String empty){
			assertEquals(0, SequenceTools.percentNucleotideSequence(empty));
		}

		@Test
		void nucleotideContentTest(){
			assertEquals(100, SequenceTools.percentNucleotideSequence("ATCGCAA"));
			assertEquals(100, SequenceTools.percentNucleotideSequence("UUACG"));
			assertEquals(100, SequenceTools.percentNucleotideSequence(randomDNA(1_000_000)));
			assertEquals(50, SequenceTools.percentNucleotideSequence("123CCG"));
			assertEquals(66, SequenceTools.percentNucleotideSequence("12TTAC"));				assertEquals(0, SequenceTools.percentNucleotideSequence("  HH"));
			assertEquals(0, SequenceTools.percentNucleotideSequence("actg"));
			}

		@Test
		void isNucleotideSequence () {
			assertTrue(SequenceTools.isNucleotideSequence("AACGAA"));
			assertFalse(SequenceTools.isNucleotideSequence("aacgaa"));
			assertFalse(SequenceTools.isNucleotideSequence("  HH"));
		}

		@ParameterizedTest
		@NullAndEmptySource
		@DisplayName("isNucleotide is false for null "+
			 "or empty string")
		void isnucleotideInvalidValues(String empty){
			assertFalse(SequenceTools.isNucleotideSequence(empty));
		}
	}
	@Nested
	@DisplayName("SequenceFromString")
	class SequenceFromString{
		SequenceTools tools = new SequenceTools();

		@Test
		void acceptsUpperCaseDNA() throws CompoundNotFoundException  {
			Sequence<?>nuc = tools.getSequenceFromString("ATCG");
			assertEquals(4, nuc.getLength());
		}

		@Test
		void acceptsLowerCaseDNA() throws CompoundNotFoundException  {
			Sequence<?>nuc = tools.getSequenceFromString("atcg");
			assertEquals(4, nuc.getLength());
		}

		@Test
		void rejectsRNA()throws CompoundNotFoundException {
			assertThrows(CompoundNotFoundException.class,
				 ()->tools.getSequenceFromString("AUCG"));
		}

		@Test
		void acceptsSingleLetterProtein()throws CompoundNotFoundException {
			Sequence<?> protein = tools.getSequenceFromString("HYDESS");
			assertEquals(6, protein.getLength());
		}

		@Test
		void interpets3LetterAACodeAsSingleLetter()throws CompoundNotFoundException {
			Sequence<?> protein = tools.getSequenceFromString("AlaGlySer");
			assertEquals(9, protein.getLength());
		}

		@EmptySource
		@ParameterizedTest
		@DisplayName("empty string return 0-length protein sequence")
		void emptyString(String empty) throws CompoundNotFoundException{
			Sequence<?> protein = tools.getSequenceFromString(empty);
			assertEquals(0, protein.getLength());
			assertTrue(protein instanceof ProteinSequence);
		}

		@NullSource
		@ParameterizedTest
		@DisplayName("null string throws NPE")
		void nullString(String nullStr) throws CompoundNotFoundException{
			assertThrows(NullPointerException.class,
			 ()-> tools.getSequenceFromString(nullStr));
		}
	}
}
