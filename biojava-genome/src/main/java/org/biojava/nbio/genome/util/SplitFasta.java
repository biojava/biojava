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

package org.biojava.nbio.genome.util;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;


/**
 * Utility to write each Fasta entry to a unique file
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SplitFasta {

	private static final Logger logger = LoggerFactory.getLogger(SplitFasta.class);

	public void processNucleotides(File fastaFileName,String uniqueid, File outputDirectory ) throws Exception{
		if(!outputDirectory.exists())
			outputDirectory.mkdirs();

		LinkedHashMap<String,DNASequence> dnaSequenceHashMap = FastaReaderHelper.readFastaDNASequence(fastaFileName);
		for(DNASequence dnaSequence : dnaSequenceHashMap.values()){
			String fileName = outputDirectory.getAbsolutePath() + File.separatorChar;
			if(uniqueid.length() > 0){
				fileName = fileName + dnaSequence.getAccession().getID() + ".fna";
			}else{
				fileName = fileName + uniqueid + dnaSequence.getAccession().getID() + ".fna";
			}
			ArrayList<DNASequence> dnaList = new ArrayList<DNASequence>();
			dnaList.add(dnaSequence);
			FastaWriterHelper.writeNucleotideSequence(new File(fileName), dnaList);
		}

	}

		public static void main( String[] args ){
		try{
			SplitFasta splitFasta = new SplitFasta();
			splitFasta.processNucleotides(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds.fna"), "", new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/individual"));
		}catch(Exception e){
			logger.error("Exception: ", e);
		}
	}

}
