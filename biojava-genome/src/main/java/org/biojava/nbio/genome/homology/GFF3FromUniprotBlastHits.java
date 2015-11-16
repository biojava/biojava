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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.genome.homology;

import org.biojava.nbio.genome.GeneFeatureHelper;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.*;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.features.DBReferenceInfo;
import org.biojava.nbio.core.sequence.features.DatabaseReferenceInterface;
import org.biojava.nbio.core.sequence.features.FeaturesKeyWordInterface;
import org.biojava.nbio.core.sequence.loader.UniprotProxySequenceReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 * @author Mark Chapman
 */
public class GFF3FromUniprotBlastHits {

	private static final Logger logger = LoggerFactory.getLogger(GFF3FromUniprotBlastHits.class);

    public void process(File xmlBlastHits, double ecutoff, LinkedHashMap<String, GeneSequence> geneSequenceHashMap, OutputStream gff3Output) throws Exception {
        LinkedHashMap<String, ArrayList<String>> hits = BlastHomologyHits.getMatches(xmlBlastHits, ecutoff);
        process(hits, geneSequenceHashMap, gff3Output);
    }

    public void process(LinkedHashMap<String, ArrayList<String>> hits, LinkedHashMap<String, GeneSequence> geneSequenceHashMap, OutputStream gff3Output) throws Exception {
        int size = hits.size();
        int index = 0;
 //       HashMap<String, String> scaffoldsReferencedHashMap = new HashMap<String, String>();
        for (String accessionid : hits.keySet()) {
            index++;
            if (index == 12) {
                index = 12;
            }
            logger.error(accessionid + " " + index + "/" + size);
            try {

                String[] data = accessionid.split(" ");
                String id = data[0];
                GeneSequence geneSequence = geneSequenceHashMap.get(id);
                if (geneSequence == null) {
                    logger.error("Not found " + id);
                    continue;
                }
                ArrayList<String> uniprotProteinHits = hits.get(accessionid);
                String uniprotBestHit = uniprotProteinHits.get(0);
                UniprotProxySequenceReader<AminoAcidCompound> uniprotSequence = new UniprotProxySequenceReader<AminoAcidCompound>(uniprotBestHit, AminoAcidCompoundSet.getAminoAcidCompoundSet());

                ProteinSequence proteinSequence = new ProteinSequence(uniprotSequence);
                String hitSequence = proteinSequence.getSequenceAsString();
                for (TranscriptSequence transcriptSequence : geneSequence.getTranscripts().values()) {


                    String predictedProteinSequence = transcriptSequence.getProteinSequence().getSequenceAsString();
                    ArrayList<ProteinSequence> cdsProteinList = transcriptSequence.getProteinCDSSequences();

                    ArrayList<CDSSequence> cdsSequenceList = new ArrayList<CDSSequence>(transcriptSequence.getCDSSequences().values());
                    String testSequence = "";
                    for (ProteinSequence cdsProteinSequence : cdsProteinList) {
                        testSequence = testSequence + cdsProteinSequence.getSequenceAsString();
                    }
                    if (!testSequence.equals(predictedProteinSequence) && (!predictedProteinSequence.equals(testSequence.substring(0, testSequence.length() - 1)))) {
                        DNASequence codingSequence = transcriptSequence.getDNACodingSequence();
                        logger.info("Coding Sequence: {}", codingSequence.getSequenceAsString());
                        logger.info("Sequence agreement error");
                        logger.info("CDS seq={}", testSequence);
                        logger.info("PRE seq={}", predictedProteinSequence);
                        logger.info("UNI seq={}", hitSequence);
                        //  throw new Exception("Protein Sequence compare error " + id);
                    }

                    SequencePair<ProteinSequence, AminoAcidCompound> alignment = Alignments.getPairwiseAlignment(
                            transcriptSequence.getProteinSequence(), proteinSequence,
                            PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(), 
                            SimpleSubstitutionMatrix.getBlosum62()
                            );
                    // System.out.println();
                    //    System.out.println(alignment.getSummary());
                    //   System.out.println(new Pair().format(alignment));
                    int proteinIndex = 0;
                    int gff3Index = 0;
                    for (int i = 0; i < cdsProteinList.size(); i++) {
                        ProteinSequence peptideSequence = cdsProteinList.get(i);
                        String seq = peptideSequence.getSequenceAsString();
                        Integer startIndex = null;
                        int offsetStartIndex = 0;
                        for (int s = 0; s < seq.length(); s++) {
                            startIndex = alignment.getIndexInTargetForQueryAt(proteinIndex + s);
                            if (startIndex != null) {
                                startIndex = startIndex + 1;
                                offsetStartIndex = s;
                                break;
                            }
                        }
                        Integer endIndex = null;

                        int offsetEndIndex = 0;
                        for (int e = 0; e < seq.length(); e++) {
                            endIndex = alignment.getIndexInTargetForQueryAt(proteinIndex + seq.length() - 1 - e);
                            if (endIndex != null) {
                                endIndex = endIndex + 1;
                                offsetEndIndex = e;
                                break;
                            }
                        }

                        proteinIndex = proteinIndex + seq.length();
                        if (startIndex != null && endIndex != null && startIndex != endIndex) {
                            CDSSequence cdsSequence = cdsSequenceList.get(i);
                            String hitLabel = "";
                            if (transcriptSequence.getStrand() == Strand.POSITIVE) {
                                hitLabel = uniprotBestHit + "_" + startIndex + "_" + endIndex;
                            } else {
                                hitLabel = uniprotBestHit + "_" + endIndex + "_" + startIndex;
                            }
                            int dnaBeginIndex = cdsSequence.getBioBegin() + (3 * offsetStartIndex);
                            int dnaEndIndex = cdsSequence.getBioEnd() - (3 * offsetEndIndex);
                            String scaffold = geneSequence.getParentChromosomeSequence().getAccession().getID();
                    //        if (scaffoldsReferencedHashMap.containsKey(scaffold) == false) {
                    //            String gff3line = scaffold + "\t" + geneSequence.getSource() + "\t" + "size" + "\t" + "1" + "\t" + geneSequence.getParentChromosomeSequence().getBioEnd() + "\t.\t.\t.\tName=" + scaffold + "\r\n";
                    //            gff3Output.write(gff3line.getBytes());
                    //            scaffoldsReferencedHashMap.put(scaffold, scaffold);
                    //        }

                            String line = scaffold + "\t" + geneSequence.getSource() + "_" + "UNIPROT\tmatch\t" + dnaBeginIndex + "\t" + dnaEndIndex + "\t.\t" + transcriptSequence.getStrand().getStringRepresentation() + "\t.\t";
                            if (gff3Index == 0) {
                                FeaturesKeyWordInterface featureKeyWords = proteinSequence.getFeaturesKeyWord();
                                String notes = "";
                                if (featureKeyWords != null) {
                                    ArrayList<String> keyWords = featureKeyWords.getKeyWords();
                                    if (keyWords.size() > 0) {
                                        notes = ";Note=";
                                        for (String note : keyWords) {
                                            if (note.equals("Complete proteome")) {
                                                continue;
                                            }
                                            if (note.equals("Direct protein sequencing")) {
                                                continue;
                                            }

                                            notes = notes + " " + note;
                                            geneSequence.addNote(note); // add note/keyword which can be output in fasta header if needed
                                        }
                                    }

                                }

                                DatabaseReferenceInterface databaseReferences = proteinSequence.getDatabaseReferences();
                                if (databaseReferences != null) {
                                    LinkedHashMap<String, ArrayList<DBReferenceInfo>> databaseReferenceHashMap = databaseReferences.getDatabaseReferences();
                                    ArrayList<DBReferenceInfo> pfamList = databaseReferenceHashMap.get("Pfam");
                                    ArrayList<DBReferenceInfo> cazyList = databaseReferenceHashMap.get("CAZy");
                                    ArrayList<DBReferenceInfo> goList = databaseReferenceHashMap.get("GO");
                                    ArrayList<DBReferenceInfo> eccList = databaseReferenceHashMap.get("BRENDA");
                                    if (pfamList != null && pfamList.size() > 0) {
                                        if (notes.length() == 0) {
                                            notes = ";Note=";
                                        }
                                        for (DBReferenceInfo note : pfamList) {
                                            notes = notes + " " + note.getId();
                                            geneSequence.addNote(note.getId()); // add note/keyword which can be output in fasta header if needed
                                        }
                                    }

                                    if (cazyList != null && cazyList.size() > 0) {
                                        if (notes.length() == 0) {
                                            notes = ";Note=";
                                        }
                                        for (DBReferenceInfo note : cazyList) {
                                            notes = notes + " " + note.getId();
                                            geneSequence.addNote(note.getId()); // add note/keyword which can be output in fasta header if needed
                                            // System.out.println("CAZy=" + note);
                                        }
                                    }

                                    if (eccList != null && eccList.size() > 0) {
                                        if (notes.length() == 0) {
                                            notes = ";Note=";
                                        }
                                        for (DBReferenceInfo note : eccList) {
                                            String dbid = note.getId();
                                            dbid = dbid.replace(".", "_"); //replace . with _ to facilitate searching in gbrowse
                                            notes = notes + " " + "EC:" + dbid;
                                            geneSequence.addNote("EC:" + dbid); // add note/keyword which can be output in fasta header if needed

                                        }
                                    }

                                    if (goList != null && goList.size() > 0) {
                                        if (notes.length() == 0) {
                                            notes = ";Note=";
                                        }
                                        for (DBReferenceInfo note : goList) {
                                            notes = notes + " " + note.getId();
                                            geneSequence.addNote(note.getId()); // add note/keyword which can be output in fasta header if needed
                                            LinkedHashMap<String, String> properties = note.getProperties();
                                            for (String propertytype : properties.keySet()) {
                                                if (propertytype.equals("evidence")) {
                                                    continue;
                                                }
                                                String property = properties.get(propertytype);

                                                if (property.startsWith("C:")) {
                                                    continue; // skip over the location
                                                }
                                                if (property.endsWith("...")) {
                                                    property = property.substring(0, property.length() - 3);
                                                }
                                                notes = notes + " " + property;
                                                geneSequence.addNote(property);
                                            }
                                        }
                                    }

                                }


                                line = line + "Name=" + hitLabel + ";Alias=" + uniprotBestHit + notes + "\n";
                            } else {
                                line = line + "Name=" + hitLabel + "\n";
                            }
                            gff3Index++;

                            gff3Output.write(line.getBytes());
                        }
                    }
                }
            } catch (Exception e) {
                logger.info("Accession Id: {}", accessionid, e);
            }
        }



    }


    

    public static void main(String[] args) {
        /*
            try {
                LinkedHashMap<String, ChromosomeSequence> dnaSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGeneMarkGTF(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds.fna"), new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/genemark_hmm.gtf"));
                LinkedHashMap<String, GeneSequence> geneSequenceList = GeneFeatureHelper.getGeneSequences(dnaSequenceList.values());
                FileOutputStream fo = new FileOutputStream("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/genemark_uniprot_match.gff3");

                GFF3FromUniprotBlastHits gff3FromUniprotBlastHits = new GFF3FromUniprotBlastHits();
                gff3FromUniprotBlastHits.process(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/c1-454Scaffolds-hits-uniprot_fungi.xml"), 1E-10, geneSequenceList, fo);
                fo.close();


            } catch (Exception e) {
                logger.error("Exception: ", e);


            }
        */

            try {
                LinkedHashMap<String, ChromosomeSequence> dnaSequenceHashMap = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGlimmerGFF3(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds-16.fna"), new File("/Users/Scooter/scripps/dyadic/GlimmerHMM/c1_glimmerhmm-16.gff"));
                LinkedHashMap<String, GeneSequence> geneSequenceList = GeneFeatureHelper.getGeneSequences(dnaSequenceHashMap.values());
                FileOutputStream fo = new FileOutputStream("/Users/Scooter/scripps/dyadic/outputGlimmer/genemark_uniprot_match-16.gff3");
                LinkedHashMap<String, ArrayList<String>> blasthits = BlastHomologyHits.getMatches(new File("/Users/Scooter/scripps/dyadic/blastresults/c1_glimmer_in_uniprot.xml"), 1E-10);
                logger.error("Number of uniprot hits " + blasthits.size());

                GFF3FromUniprotBlastHits gff3FromUniprotBlastHits = new GFF3FromUniprotBlastHits();
                gff3FromUniprotBlastHits.process(blasthits, geneSequenceList, fo);
                fo.close();
            } catch (Exception e) {
            	logger.error("Exception: ", e);
            }
    }
}
