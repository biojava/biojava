package org.biojava3.ronn;

import static org.junit.Assert.*;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.biojava3.data.sequence.FastaSequence;
import org.biojava3.data.sequence.SequenceUtil;
import org.biojava3.ronn.Jronn.Range;
import org.biojava3.ronn.ORonn.ResultLayout;
import org.junit.Test;


public class JronnExample {

	
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
		Map<FastaSequence, float[]>	rawProbabilityScores = Jronn.getDisorderScores(sequences); 
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	@Test
	public void rawScoresForSingleSequence() {
		float[]	rawProbabilityScores = Jronn.getDisorderScores(new FastaSequence("name", "LLRGRHLMNGTMIMRPWNFLNDHHFPKFFPHLIEQQAIWLADWWRKKHC" +
				"RPLPTRAPTMDQWDHFALIQKHWTANLWFLTFPFNDKWGWIWFLKDWTPGSADQAQRACTWFFCHGHDTN" +
				"CQIIFEGRNAPERADPMWTGGLNKHIIARGHFFQSNKFHFLERKFCEMAEIERPNFTCRTLDCQKFPWDDP" ));
	}
	
	
	@Test
	public void disorderForMultipleSequences() {
		try {
			final List<FastaSequence> sequences = SequenceUtil.readFasta(new FileInputStream("src/test/resources/fasta.in"));
		 	Map<FastaSequence, Range[]> ranges = Jronn.getDisorder(sequences);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	
	@Test
	public void disorderForSingleSequence() {
		Range[]	ranges = Jronn.getDisorder(new FastaSequence("name", "LLRGRHLMNGTMIMRPWNFLNDHHFPKFFPHLIEQQAIWLADWWRKKHC" +
				"RPLPTRAPTMDQWDHFALIQKHWTANLWFLTFPFNDKWGWIWFLKDWTPGSADQAQRACTWFFCHGHDTN" +
				"CQIIFEGRNAPERADPMWTGGLNKHIIARGHFFQSNKFHFLERKFCEMAEIERPNFTCRTLDCQKFPWDDP" ));
		
	}
	
	
	
	
}
