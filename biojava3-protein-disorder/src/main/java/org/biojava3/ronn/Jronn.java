package org.biojava3.ronn;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;



import org.biojava3.data.sequence.FastaSequence;
import org.biojava3.data.sequence.SequenceUtil;

public class Jronn {
	
	// Load models
	private static final ModelLoader loader = new ModelLoader();  
	static {
		try {
			loader.loadModels();
		} catch (NumberFormatException e) {
			throw new RuntimeException("Fails to load models!" + e.getMessage(), e); 
		} catch (IOException e) {
			throw new RuntimeException("Fails to load models!" + e.getMessage(), e);
		}
	}
	
	private static int availCpus = Runtime.getRuntime().availableProcessors(); 
	
	
	final int numberOfThreads;  
	
	public Jronn(int numberOfThreads) {
		if(numberOfThreads<=0) {
			throw new IllegalArgumentException("Number of threads must be greater then 0!");
		}
		if(numberOfThreads>availCpus*2) {
			System.err.println("You may be waisting your time by using "+availCpus+" threads! " +
					"You must know what you are doing!");
		}
		this.numberOfThreads = numberOfThreads; 
	}
	
	// Value class
	public static class Range {

		final int from; 
		final int to; 
	
		public Range(int from, int to) {
			assert from>=0; 
			assert from<to; 
			this.from = from; 
			this.to = to; 
		}

		@Override
		public String toString() { 
			return "Range" + " From:" + from + "\t" + "to: " + to + "\n";
		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + from;
			result = prime * result + to;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Range other = (Range) obj;
			if (from != other.from)
				return false;
			if (to != other.to)
				return false;
			return true;
		}
		
		
	}
	
	public float[] getDisorderScores(FastaSequence sequence) {
		    return predictSerial(sequence);
	}

	private static float[] predictSerial(FastaSequence fsequence) {
		ORonn.validateSequenceForRonn(fsequence);
		InputParameters in = new InputParameters();
		ORonn ronn;
		float[] disorder = null; 
		try {
			ronn = new ORonn(fsequence, loader, in);
			disorder = ronn.call().getMeanScores();
		} catch (NumberFormatException e) {
			throw new RuntimeException("Jronn fails to load models " + e.getLocalizedMessage(), e);
		} catch (IOException e) {
			throw new RuntimeException("Jronn fails to load models " + e.getLocalizedMessage(), e);
		}
		return disorder;  
	}

	public Range[] getDisorder(FastaSequence sequence) {
		float[] scores = getDisorderScores(sequence);
		return scoresToRanges(scores, RonnConstraint.DEFAULT_RANGE_PROBABILITY_THRESHOLD);
	}

	/**
	 * Convert raw scores to ranges. 
	 * @param scores
	 * @param probability
	 * @return
	 */
	public static Range[] scoresToRanges(float[] scores, float probability)  {
		assert scores!=null && scores.length>0;
		assert probability>0 && probability<1;
		
		int count=0;
		int regionLen=0;
		List<Range> ranges = new ArrayList<Range>();
		for(float score: scores) { 
			count++;
			// Round to 2 decimal points before comparison 
			score = (float) (Math.round(score*100.0)/100.0);
			if(score>probability) {
				regionLen++;
			} else {
				if(regionLen>0) {
					ranges.add(new Range(count-regionLen, count-1));
				}
				regionLen=0;
			}
		}
		// In case of the range to boundary runs to the very end of the sequence 
		if(regionLen>1) {
			ranges.add(new Range(count-regionLen+1, count));
		}
		return ranges.toArray(new Range[ranges.size()]); 		

	}
	public Map<FastaSequence,float[]> getDisorderScores(List<FastaSequence> sequences) {
		Map<FastaSequence,float[]> results = new TreeMap<FastaSequence, float[]>();
		for(FastaSequence fsequence : sequences) {
			results.put(fsequence, predictSerial(fsequence));
		}
		return results; 
	}
	
	
	public Map<FastaSequence,Range[]> getDisorder(List<FastaSequence> sequences) {
		Map<FastaSequence,Range[]> disorderRanges = new TreeMap<FastaSequence,Range[]>();
		for(FastaSequence fs: sequences) {
			disorderRanges.put(fs, getDisorder(fs));
		}
		return disorderRanges; 
	}
	
	public Map<FastaSequence,Range[]> getDisorder(String fastaFile) throws FileNotFoundException, IOException {
		final List<FastaSequence> sequences = SequenceUtil.readFasta(new FileInputStream(fastaFile));
		return getDisorder(sequences);
	}
	
	/*
	public void writeDisorder(String fastaFile, String outputFile, ResultLayout layout) throws FileNotFoundException, IOException {
		final List<FastaSequence> sequences = SequenceUtil.readFasta(new FileInputStream(fastaFile));
		ORonn.predictParallel(sequences, prms, loader); 
	} */
}
