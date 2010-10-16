/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.genome.parsers.gff;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import org.biojava3.core.sequence.CDSSequence;
import org.biojava3.core.sequence.ChromosomeSequence;
import org.biojava3.core.sequence.GeneSequence;
import org.biojava3.core.sequence.SequenceComparator;
import org.biojava3.core.sequence.TranscriptSequence;
import org.biojava3.genome.GeneFeatureHelper;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GFF3Writer {

    /**
     * Output gff3 format for a DNA Sequence
     * @param fileName
     * @param chromosomeSequence
     * @throws Exception
     */
    public void write(OutputStream outputStream, LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList) throws Exception {

        outputStream.write("##gff-version 3\n".getBytes());
        for (String key : chromosomeSequenceList.keySet()) {
            ChromosomeSequence chromosomeSequence = chromosomeSequenceList.get(key);
            String gff3line = "";
   //         if(source.length() == 0){
   //             Collection<GeneSequence> genes = chromosomeSequence.getGeneSequences().values();
   //             for(GeneSequence gene : genes){
   //                 source = gene.getSource();
   //                 break;
   //             }
   //         }
   //         gff3line = key + "\t" + source + "\t" + "size" + "\t" + "1" + "\t" + chromosomeSequence.getBioEnd() + "\t.\t.\t.\tName=" + key + "\r\n";
   //         outputStream.write(gff3line.getBytes());

            for (GeneSequence geneSequence : chromosomeSequence.getGeneSequences().values()) {
                gff3line = key + "\t" + geneSequence.getSource() + "\t" + "gene" + "\t" + geneSequence.getBioBegin() + "\t" + geneSequence.getBioEnd() + "\t";
                Double score = geneSequence.getSequenceScore();
                if (score == null) {
                    gff3line = gff3line + ".\t";
                } else {
                    gff3line = gff3line + score + "\t";
                }
                gff3line = gff3line + geneSequence.getStrand().getStringRepresentation() + "\t";
                gff3line = gff3line + ".\t";
                gff3line = gff3line + "ID=" + geneSequence.getAccession().getID() + ";Name=" + geneSequence.getAccession().getID();
                gff3line = gff3line + getGFF3Note(geneSequence.getNotesList());
                gff3line = gff3line + "\n";
                outputStream.write(gff3line.getBytes());

                int transcriptIndex = 0;
                for (TranscriptSequence transcriptSequence : geneSequence.getTranscripts().values()) {
                    transcriptIndex++;

                    gff3line = key + "\t" + transcriptSequence.getSource() + "\t" + "mRNA" + "\t" + transcriptSequence.getBioBegin() + "\t" + transcriptSequence.getBioEnd() + "\t";
                    score = transcriptSequence.getSequenceScore();
                    if (score == null) {
                        gff3line = gff3line + ".\t";
                    } else {
                        gff3line = gff3line + score + "\t";
                    }
                    gff3line = gff3line + transcriptSequence.getStrand().getStringRepresentation() + "\t";
                    gff3line = gff3line + ".\t";
                    String id = geneSequence.getAccession().getID() + "." + transcriptIndex;
                    gff3line = gff3line + "ID=" + id + ";Parent=" + geneSequence.getAccession().getID() + ";Name=" + id;
                    gff3line = gff3line + getGFF3Note(transcriptSequence.getNotesList());

                    gff3line = gff3line + "\n";
                    outputStream.write(gff3line.getBytes());

                    String transcriptParentName = geneSequence.getAccession().getID() + "." + transcriptIndex;
                    ArrayList<CDSSequence> cdsSequenceList = new ArrayList(transcriptSequence.getCDSSequences().values());
                    Collections.sort(cdsSequenceList, new SequenceComparator());
                    for (CDSSequence cdsSequence : cdsSequenceList) {
                        gff3line = key + "\t" + cdsSequence.getSource() + "\t" + "CDS" + "\t" + cdsSequence.getBioBegin() + "\t" + cdsSequence.getBioEnd() + "\t";
                        score = cdsSequence.getSequenceScore();
                        if (score == null) {
                            gff3line = gff3line + ".\t";
                        } else {
                            gff3line = gff3line + score + "\t";
                        }
                        gff3line = gff3line + cdsSequence.getStrand().getStringRepresentation() + "\t";
                        gff3line = gff3line + cdsSequence.getPhase() + "\t";
                        gff3line = gff3line + "ID=" + cdsSequence.getAccession().getID() + ";Parent=" + transcriptParentName;
                        gff3line = gff3line + getGFF3Note(cdsSequence.getNotesList());

                        gff3line = gff3line + "\n";
                        outputStream.write(gff3line.getBytes());
                    }

                }
            }

        }


    }

    private String getGFF3Note(ArrayList<String> notesList) {
        String notes = "";

        if (notesList.size() > 0) {
            notes = ";Note=";
            int noteindex = 1;
            for (String note : notesList) {
                notes = notes + note;
                if (noteindex < notesList.size() - 1) {
                    notes = notes + " ";
                }
            } 

        }
        return notes;
    }

    public static void main(String args[]) throws Exception {

        if (true) {
            FileOutputStream fo = new FileOutputStream("/Users/Scooter/scripps/dyadic/geneid/geneid/c1-geneid.gff3");//-16
            LinkedHashMap<String, ChromosomeSequence> dnaSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGeneIDGFF2(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds.fna"), new File("/Users/Scooter/scripps/dyadic/geneid/geneid/c1_geneid.gff"));
            GFF3Writer gff3Writer = new GFF3Writer();
            gff3Writer.write(fo, dnaSequenceList);


     //       LinkedHashMap<String, ProteinSequence> proteinSequenceList = GeneFeatureHelper.getProteinSequences(chromosomeSequenceList.values());
     //       for(String id : proteinSequenceList.keySet()){
     //           ProteinSequence sequence = proteinSequenceList.get(id);
     //           System.out.println(id + " " + sequence.getSequenceAsString());

     //       }
            fo.close();
        }

        if (false) {
            FileOutputStream fo = new FileOutputStream("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/genemark_hmm.gff3");//-16
            LinkedHashMap<String, ChromosomeSequence> dnaSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGeneMarkGTF(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds.fna"), new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/genemark_hmm.gtf"));
            GFF3Writer gff3Writer = new GFF3Writer();
            gff3Writer.write(fo, dnaSequenceList);
            fo.close();
        }

        if (false) {
            LinkedHashMap<String, ChromosomeSequence> dnaSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGlimmerGFF3(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds-16.fna"), new File("/Users/Scooter/scripps/dyadic/GlimmerHMM/c1_glimmerhmm-16.gff"));
            GFF3Writer gff3Writer = new GFF3Writer();
            gff3Writer.write(System.out, dnaSequenceList);
        }
//        System.out.println(listGenes);
        //	GeneMarkGTF.write( list, args[1] );
    }
}
