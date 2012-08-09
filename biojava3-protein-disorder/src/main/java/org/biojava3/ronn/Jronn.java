 /*        BioJava development code
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
package org.biojava3.ronn;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;



import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.data.sequence.FastaSequence;
import org.biojava3.data.sequence.SequenceUtil;


/**
 * This class gives public API to RONN functions. 
 * It is build on top of the command line client. Due to this fact a few things 
 * could be improved and extended pending the refactoring of the command line client.  
 *
 * The input sequence limitations - the input sequence must not contain any ambiguous characters, 
 * and have a minimum length of 19 amino acids. 
 * 
 * @author Peter Troshin
 * @version 1.0
 * @since 3.0.2
 * 
 */
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
	
	
	/**
	 * Holder for the ranges, contain pointers to starting and ending position 
	 * on the sequence which comprises a disordered region. Immutable. 
	 * @author pvtroshin
	 */
	public static class Range {
		/**
		 * Range starting position counts from 1 (the first position on the sequence is 1)
		 */
		public final int from; 
		/**
		 * The range ending position includes the last residue. 
		 */
		public final int to; 
	
		public final float score;
		public Range(int from, int to, float score) {
			assert from>=0; 
			assert from<to; 
			this.from = from; 
			this.to = to; 
			this.score = score;
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
	
	/**
	 * Calculates the probability value for each residue in the protein sequence, 
	 * telling the probability that the residue belongs to disordered region. 
	 * In general, values greater than 0.5 considered to be in the disordered regions. 
	 *   
	 * @param sequence an instance of FastaSequence object, holding the name and the sequence. 
	 * @return the probability scores for each residue in the sequence
	 */
	public static float[] getDisorderScores(FastaSequence sequence) {
		    return predictSerial(sequence);
	}
	
	/**
	 * Calculates the probability value for each residue in the protein sequence, 
	 * telling the probability that the residue belongs to disordered region. 
	 * In general, values greater than 0.5 considered to be in the disordered regions. 
	 *   
	 * @param sequence an instance of FastaSequence object, holding the name and the sequence. 
	 * @return the probability scores for each residue in the sequence
	 */
	public static float[] getDisorderScores(ProteinSequence sequence) {
		
		FastaSequence seq = convertProteinSequencetoFasta(sequence);
		
		return predictSerial(seq);
	}
	
	/** Utility method to convert a BioJava ProteinSequence object to the FastaSequence 
	 *  object used internally in JRonn.
	 * 
	 * @param sequence
	 * @return
	 */
	public static FastaSequence convertProteinSequencetoFasta(ProteinSequence sequence){
		StringBuffer buf = new StringBuffer();
		for (AminoAcidCompound compound : sequence) {
			
			String c = compound.getShortName();
			
			if (! SequenceUtil.NON_AA.matcher(c).find()) {
				buf.append(c);
			} else {				
				buf.append("X");
			}									
		}
		
		return new FastaSequence(sequence.getAccession().getID(),buf.toString());
	}

	private static float[] predictSerial(FastaSequence fsequence) {
		ORonn.validateSequenceForRonn(fsequence);
		ORonn ronn;
		float[] disorder = null; 
		try {
			ronn = new ORonn(fsequence, loader);
			disorder = ronn.call().getMeanScores();
		} catch (NumberFormatException e) {
			throw new RuntimeException("Jronn fails to load models " + e.getLocalizedMessage(), e);
		} catch (IOException e) {
			throw new RuntimeException("Jronn fails to load models " + e.getLocalizedMessage(), e);
		}
		return disorder;  
	}

	/**
	 * Calculates the disordered regions of the sequence. More formally, the regions for which the 
	 * probability of disorder is greater then 0.50.  
	 *  
	 *   
	 * @param sequence an instance of FastaSequence object, holding the name and the sequence.
	 * @return the array of ranges if there are any residues predicted to have the 
	 * probability of disorder greater then 0.5, null otherwise. 
	 *
	 */
	public static Range[] getDisorder(FastaSequence sequence) {
		float[] scores = getDisorderScores(sequence);
		return scoresToRanges(scores, RonnConstraint.DEFAULT_RANGE_PROBABILITY_THRESHOLD);
	}

	/**
	 * Convert raw scores to ranges. Gives ranges for given probability of disorder value 
	 * @param scores the raw probability of disorder scores for each residue in the sequence.  
	 * @param probability the cut off threshold. Include all residues with the probability of disorder greater then this value
	 * @return the array of ranges if there are any residues predicted to have the 
	 * probability of disorder greater then {@code probability}, null otherwise.
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
					ranges.add(new Range(count-regionLen, count-1,score));
				}
				regionLen=0;
			}
		}
		// In case of the range to boundary runs to the very end of the sequence 
		if(regionLen>1) {
			ranges.add(new Range(count-regionLen+1, count,scores[scores.length-1]));
		}
		return ranges.toArray(new Range[ranges.size()]); 		

	}
	
	/**
	 * Calculates the probability of disorder scores for each residue in the sequence for 
	 * many sequences in the input.
	 * 
	 * @param sequences the list of the FastaSequence objects 
	 * @return the Map with key->FastaSequence, value->probability of disorder for each residue
	 * @see #getDisorder(FastaSequence)
	 */
	public static Map<FastaSequence,float[]> getDisorderScores(List<FastaSequence> sequences) {
		Map<FastaSequence,float[]> results = new TreeMap<FastaSequence, float[]>();
		for(FastaSequence fsequence : sequences) {
			results.put(fsequence, predictSerial(fsequence));
		}
		return results; 
	}
	
	/**
	 * Calculates the disordered regions of the sequence for many sequences in the input.
	 * 
	 * @param sequences sequences the list of the FastaSequence objects
	 * @return
	 * @see #getDisorder(FastaSequence)
	 */
	public static Map<FastaSequence,Range[]> getDisorder(List<FastaSequence> sequences) {
		Map<FastaSequence,Range[]> disorderRanges = new TreeMap<FastaSequence,Range[]>();
		for(FastaSequence fs: sequences) {
			disorderRanges.put(fs, getDisorder(fs));
		}
		return disorderRanges; 
	}
	
	/**
	 * Calculates the disordered regions of the protein sequence.
	 * @param fastaFile input file name containing the sequence in FASTA
	 * @return the Map with key->FastaSequence, value->the list of disordered regions for each sequence
	 * @throws FileNotFoundException if the input file cannot be found
	 * @throws IOException of the system cannot access or read from the input file 
	 * @see #getDisorder(FastaSequence)
	 * @see #Jronn.Range
	 */
	public static Map<FastaSequence,Range[]> getDisorder(String fastaFile) throws FileNotFoundException, IOException {
		final List<FastaSequence> sequences = SequenceUtil.readFasta(new FileInputStream(fastaFile));
		return getDisorder(sequences);
	}
	
	/**
	 * TODO 
	 * 
	 * High performance method for calculating disorder. Use multiple threads to achieve the speedup.
	 *  
	 * @param fastaFile  fully qualified path to the input FASTA file  
	 * @param outputFile file name of the file for the results 
	 * @param threadNumber the number of threads to use, default
	 * @param controls the format of the result file 
	 * @throws FileNotFoundException if input file in not found 
	 * @throws IOException if the input or the output files cannot be accessed  
	 * @see ORonn.ResultLayout
	 
	public static void calculateDisorder(String fastaFile, String outputFile, int threadNumber, ResultLayout layout) throws FileNotFoundException, IOException {
		final List<FastaSequence> sequences = SequenceUtil.readFasta(new FileInputStream(fastaFile));
		InputParameters in = new InputParameters(); 
		in.setFilePrm(fastaFile, InputParameters.inputKey);
		in.setFilePrm(outputFile, InputParameters.outputKey);
		//in.setThreadNum(Integer.toString(threadNumber)); 
		ORonn.predictParallel(sequences, in, loader); 
	}
	*/ 
}
