/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.genome.util;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.io.FastaWriterHelper;


/**
 * Utility to write each Fasta entry to a unique file
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SplitFasta {

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
            e.printStackTrace();
        }
    }

}
