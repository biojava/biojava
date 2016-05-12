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
package org.biojava.nbio.ronn;

import org.biojava.nbio.data.sequence.FastaSequence;
import org.biojava.nbio.data.sequence.SequenceUtil;
import org.biojava.nbio.ronn.Jronn.Range;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.fail;


public class JronnExample {

	private static final Logger logger = LoggerFactory.getLogger(JronnExample.class);

	/*
	@Test
	public void highPerformanceDisorderCalculation() {
		try {
			Jronn.calculateDisorder("src/test/resources/fasta.in", "src/test/resources/result.txt", 4, ResultLayout.HORIZONTAL);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	*/

	@Test
	public void rawScoresForMultipleSequences() {
		try {
		final List<FastaSequence> sequences = SequenceUtil.readFasta(new FileInputStream("src/test/resources/fasta.in"));
		@SuppressWarnings("unused")
		Map<FastaSequence, float[]>	rawProbabilityScores = Jronn.getDisorderScores(sequences);
		} catch (FileNotFoundException e) {
			logger.error("FileNotFoundException: ", e);
			fail(e.getMessage());
		} catch (IOException e) {
			logger.error("IOException: ", e);
			fail(e.getMessage());
		}
	}

	@Test
	public void rawScoresForSingleSequence() {
		@SuppressWarnings("unused")
		float[]	rawProbabilityScores = Jronn.getDisorderScores(new FastaSequence("name", "LLRGRHLMNGTMIMRPWNFLNDHHFPKFFPHLIEQQAIWLADWWRKKHC" +
				"RPLPTRAPTMDQWDHFALIQKHWTANLWFLTFPFNDKWGWIWFLKDWTPGSADQAQRACTWFFCHGHDTN" +
				"CQIIFEGRNAPERADPMWTGGLNKHIIARGHFFQSNKFHFLERKFCEMAEIERPNFTCRTLDCQKFPWDDP" ));
	}


	@Test
	public void disorderForMultipleSequences() {
		try {
			final List<FastaSequence> sequences = SequenceUtil.readFasta(new FileInputStream("src/test/resources/fasta.in"));
		 	@SuppressWarnings("unused")
			Map<FastaSequence, Range[]> ranges = Jronn.getDisorder(sequences);

		} catch (FileNotFoundException e) {
			logger.error("FileNotFoundException: ", e);
			fail(e.getMessage());
		} catch (IOException e) {
			logger.error("IOException: ", e);
			fail(e.getMessage());
		}
	}


	@Test
	public void disorderForSingleSequence() {
		@SuppressWarnings("unused")
		Range[]	ranges = Jronn.getDisorder(new FastaSequence("name", "LLRGRHLMNGTMIMRPWNFLNDHHFPKFFPHLIEQQAIWLADWWRKKHC" +
				"RPLPTRAPTMDQWDHFALIQKHWTANLWFLTFPFNDKWGWIWFLKDWTPGSADQAQRACTWFFCHGHDTN" +
				"CQIIFEGRNAPERADPMWTGGLNKHIIARGHFFQSNKFHFLERKFCEMAEIERPNFTCRTLDCQKFPWDDP" ));

	}




}
