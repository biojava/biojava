package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertAll;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

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
					()->assertEquals(expected, SequenceTools.permuteCyclic2(original, Integer.MAX_VALUE)),
					()->assertEquals(expected, SequenceTools.permuteCyclic2(original, Integer.MIN_VALUE))
				);
		}
		
		@Test
		void permuteCyclicPerformance() {
			String dna = randomDNA(10_000_000);
			long start = System.currentTimeMillis();
			String rotated = SequenceTools.permuteCyclic(dna, 5_000_000);
			long end = System.currentTimeMillis();
			System.err.println(end-start);
			
			long start2 = System.currentTimeMillis();
			String rotated2 = SequenceTools.permuteCyclic2(dna, 5_000_000);
			long end2 = System.currentTimeMillis();
			System.err.println(end2-start2);
			assertTrue((end-start)/(end2-start2) > 5);
		}

	}

}
