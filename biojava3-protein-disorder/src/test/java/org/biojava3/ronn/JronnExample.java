package org.biojava3.ronn;

import static org.junit.Assert.fail;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.biojava3.data.sequence.FastaSequence;
import org.biojava3.data.sequence.SequenceUtil;
import org.biojava3.ronn.Jronn.Range;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


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
