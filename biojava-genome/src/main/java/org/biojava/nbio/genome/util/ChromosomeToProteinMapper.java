package org.biojava.nbio.genome.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.biojava.nbio.genome.parsers.genename.GeneChromosomePosition;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;

import com.google.common.collect.Range;

public class ChromosomeToProteinMapper {
	
	private TwoBitParser parser;
	private String GENOME_URI;
	
	/** Sets a path to the genome.2bit file.
	 * 
	 *  @param path The path to the file containing the genome in .2bit format.
	 */
	public void setGenomeURI(String path) {
		GENOME_URI=path;
	}
	
	/** 
	 *  Reads a genome from a locally stored .2bit file.
	 */
	public void readGenome() throws Exception {
		File f = new File(GENOME_URI);
		this.parser = new TwoBitParser(f);
	}
	
	/** Sets a chromosome for TwoBitParser.
	 * 
	 * @param chr The chromosome name (e.g. chr21)
	 */
	public void setChromosome(String chr) throws Exception {
		parser.close();
		String[] names = parser.getSequenceNames();
		for(int i=0;i<names.length;i++) {
			if ( names[i].equals(chr) ) {
				parser.setCurrentSequence(names[i]);
				break;
			}
		}
	}
	
    /** Extracts the exons boundaries in CDS coordinates for genes living on the reverse DNA strand.
     *
     * @param exonStarts The list holding the genetic coordinates pointing to the start positions of the exons (including UTR regions)  
     * @param exonEnds The list holding the genetic coordinates pointing to the end positions of the exons (including UTR regions)
     * @param cdsStart The start position of a coding region
     * @param cdsEnd The end position of a coding region
     * 
     * @return the list of genetic positions corresponding to the exons boundaries in CDS coordinates
    */
    public static List<Range<Integer>> getCDSRegionsReverse(List<Integer> exonStarts, List<Integer> exonEnds,
            int cdsStart, int cdsEnd) {

        // remove exons that are fully landed in UTRs
        List<Integer> tmpS = new ArrayList<Integer>(exonStarts);
        List<Integer> tmpE = new ArrayList<Integer>(exonEnds);
        
        int j=0;
        for (int i = 0; i < tmpS.size(); i++) {
        	if ( ( tmpE.get(i) < cdsStart) || ( tmpS.get(i) > cdsEnd) ) {
        		exonStarts.remove(j);
        		exonEnds.remove(j);
        	}
        	else {
        		j++;
        	}
        }
        
        // remove untranslated regions from exons
        int nExons = exonStarts.size();
        exonStarts.remove(0);
        exonStarts.add(0, cdsStart);
        exonEnds.remove(nExons-1);
        exonEnds.add(cdsEnd);
        
        List<Range<Integer>> cdsRegion = new ArrayList<Range<Integer>>();
        for ( int i=0; i<nExons; i++ ) {
        	Range<Integer> r = Range.closed(exonStarts.get(i), exonEnds.get(i));
        	cdsRegion.add(r);
        }
		return cdsRegion;
    }
    
    /** Extracts the exons boundaries in CDS coordinates for genes living on the forward DNA strand.
    *
    * @param exonStarts The list holding the genetic coordinates pointing to the start positions of the exons (including UTR regions)  
    * @param exonEnds The list holding the genetic coordinates pointing to the end positions of the exons (including UTR regions)
    * @param cdsStart The start position of a coding region
    * @param cdsEnd The end position of a coding region
    * 
    * @return the list of genetic positions corresponding to the exons boundaries in CDS coordinates
   */
    public static List<Range<Integer>> getCDSRegionsForward(List<Integer> exonStarts, List<Integer> exonEnds,
            int cdsStart, int cdsEnd) {
    	
        // remove exons that are fully landed in UTRs
        List<Integer> tmpS = new ArrayList<Integer>(exonStarts);
        List<Integer> tmpE = new ArrayList<Integer>(exonEnds);
        
        int j=0;
        for (int i = 0; i < tmpS.size(); i++) {
        	if ( ( tmpE.get(i) < cdsStart) || ( tmpS.get(i) > cdsEnd) ) {
        		exonStarts.remove(j);
        		exonEnds.remove(j);
        	}
        	else {
        		j++;
        	}
        }
        
        // remove untranslated regions from exons
        int nExons = exonStarts.size();
        exonStarts.remove(0);
        exonStarts.add(0, cdsStart);
        exonEnds.remove(nExons-1);
        exonEnds.add(cdsEnd);
    	
        List<Range<Integer>> cdsRegion = new ArrayList<Range<Integer>>();
        for ( int i=0; i<nExons; i++ ) {
        	Range<Integer> r = Range.closed(exonStarts.get(i), exonEnds.get(i));
        	cdsRegion.add(r);
        }
		return cdsRegion;
    }
    
    /** Extracts the DNA sequence transcribed from the input genetic coordinates.
    *
    * @param gcp The container with chromosomal positions
    * 
    * @return the DNA sequence transcribed from the input genetic coordinates
   */
    public String getTranscriptSequence(GeneChromosomePosition gcp) throws IOException, CompoundNotFoundException {
    	return getTranscriptSequence(gcp.getExonStarts(), gcp.getExonEnds(), gcp.getCdsStart(), gcp.getCdsEnd(), gcp.getOrientation());
    }
    
    /** Extracts the DNA sequence transcribed from the input genetic coordinates.
    *
    * @param exonStarts The list holding the genetic coordinates pointing to the start positions of the exons (including UTR regions)  
    * @param exonEnds The list holding the genetic coordinates pointing to the end positions of the exons (including UTR regions)
    * @param cdsStart The start position of a coding region
    * @param cdsEnd The end position of a coding region
    * @param orientation The orientation of the strand where the gene is living
    * 
    * @return the DNA sequence transcribed from the input genetic coordinates
    */
	public String getTranscriptSequence(List<Integer> exonStarts, List<Integer> exonEnds, int codingStart, int codingEnd, Character orientation) throws IOException, CompoundNotFoundException {

		List<Range<Integer>> cdsRegion;
		if (orientation.equals("-")) {
			cdsRegion = getCDSRegionsReverse(exonStarts, exonEnds, codingStart, codingEnd);
		}
		else {
			cdsRegion = getCDSRegionsForward(exonStarts, exonEnds, codingStart, codingEnd);
		}

		String transcription = "";
		for (Range<Integer> range : cdsRegion) {
			int length = range.upperEndpoint() - range.lowerEndpoint();
			transcription += parser.loadFragment(range.lowerEndpoint(), length);
		}
		if (orientation.equals("-")) {
			transcription = new StringBuilder(transcription).reverse().toString();
			DNASequence dna = new DNASequence(transcription);
			SequenceView<NucleotideCompound> compliment = dna.getComplement();
			transcription = compliment.getSequenceAsString();
		}
		return transcription;
	}
	
    /** Converts the DNA sequence to protein sequence.
    *
    * @param sequence the DNA sequence
    * 
    * @return the protein sequence
    */
	public String convertDNAtoProteinSequence(String sequence) throws CompoundNotFoundException {
		DNASequence dna = new DNASequence(sequence);
		RNASequence mRNA = dna.getRNASequence();
		return mRNA.getProteinSequence().toString();
	}
}
