/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.genome;

import java.io.File;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.logging.Logger;
import org.biojava3.genome.parsers.gff.Feature;

import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ExonSequence;
import org.biojava3.core.sequence.FastaReaderHelper;
import org.biojava3.core.sequence.FastaWriterHelper;
import org.biojava3.core.sequence.GeneSequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.TranscriptSequence;
import org.biojava3.genome.parsers.gff.GFF3;



import org.biojava3.genome.parsers.gff.GeneMarkGTF;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GeneFeatureHelper {

    private static final Logger log = Logger.getLogger(GeneFeatureHelper.class.getName());

    static public LinkedHashMap<String, DNASequence> loadFastaAddGeneFeaturesFromGFF3(File fastaSequenceFile, File gffFile) throws Exception {
        LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile);
        FeatureList listGenes = GFF3.read(gffFile.getAbsolutePath());
        dnaSequenceList = addGFF3GeneFeatures(dnaSequenceList, listGenes);
        return dnaSequenceList;
    }

    static public LinkedHashMap<String, DNASequence> addGFF3GeneFeatures(LinkedHashMap<String, DNASequence> dnaSequenceList, FeatureList listGenes) throws Exception {
        FeatureList mRNAFeatures = listGenes.selectByType("mRNA");
        for (FeatureI f : mRNAFeatures) {
            Feature mRNAFeature = (Feature) f;
            String geneid = mRNAFeature.getAttribute("ID");


            FeatureList gene = listGenes.selectByAttribute("Parent", geneid);
            FeatureI geneFeature = gene.get(0);
            DNASequence seq = dnaSequenceList.get(geneFeature.seqname());
            AccessionID geneAccessionID = new AccessionID(geneid);
            GeneSequence geneSequence = null;

            FeatureList cdsFeatures = gene.selectByType("CDS");
            FeatureI feature = cdsFeatures.get(0);
            Strand strand = Strand.POSITIVE;

            if (feature.location().isNegative()) {
                strand = strand.NEGATIVE;
            }
            cdsFeatures = cdsFeatures.sortByStart();







            String seqName = feature.seqname();
            FeatureI startCodon = null;
            FeatureI stopCodon = null;
            Integer startCodonBegin = null;
            Integer stopCodonEnd = null;
            String startCodonName = "";
            String stopCodonName = "";
            FeatureList startCodonList = gene.selectByAttribute("Note", "initial-exon");
            if (startCodonList != null && startCodonList.size() > 0) {
                startCodon = startCodonList.get(0);
                startCodonBegin = startCodon.location().bioStart();
                startCodonName = startCodon.getAttribute("ID");
            }

            FeatureList stopCodonList = gene.selectByAttribute("Note", "final-exon");

            if (stopCodonList != null && stopCodonList.size() > 0) {
                stopCodon = stopCodonList.get(0);
                stopCodonEnd = stopCodon.location().bioEnd();
                stopCodonName = stopCodon.getAttribute("ID");

            }




            if (startCodonBegin == null) {
                FeatureI firstFeature = cdsFeatures.get(0);
                startCodonBegin = firstFeature.location().bioStart();
            }

            if (stopCodonEnd == null) {
                FeatureI lastFeature = cdsFeatures.get(cdsFeatures.size() - 1);
                stopCodonEnd = lastFeature.location().bioEnd();
            }



            AccessionID transcriptAccessionID = new AccessionID(geneid);
            if (geneSequence == null) {
                geneSequence = seq.addGene(geneAccessionID, startCodonBegin, stopCodonEnd);
            } else {
                //if multiple transcripts for one gene make sure the gene is defined as the min and max start/end
                if (strand.equals(Strand.POSITIVE)) {
                    if (startCodonBegin < geneSequence.getBioBegin()) {
                        geneSequence.setBioBegin(startCodonBegin);
                    }
                    if (stopCodonEnd > geneSequence.getBioBegin()) {
                        geneSequence.setBioEnd(stopCodonEnd);
                    }
                } else {
                    if (startCodonBegin > geneSequence.getBioBegin()) {
                        geneSequence.setBioBegin(startCodonBegin);
                    }
                    if (stopCodonEnd < geneSequence.getBioBegin()) {
                        geneSequence.setBioEnd(stopCodonEnd);
                    }
                }
            }
            TranscriptSequence transcriptSequence = geneSequence.addTranscript(transcriptAccessionID, startCodonBegin, stopCodonEnd, strand);
            if (startCodon != null) {
                if (startCodonName == null || startCodonName.length() == 0) {
                    startCodonName = geneid + "-start_codon-" + startCodon.location().bioStart() + "-" + startCodon.location().bioEnd();
                }
                transcriptSequence.addStartCodonSequence(new AccessionID(startCodonName), startCodon.location().bioStart(), startCodon.location().bioEnd());
            }
            if (stopCodon != null) {
                if (stopCodonName == null || stopCodonName.length() == 0) {
                    stopCodonName = geneid + "-stop_codon-" + stopCodon.location().bioStart() + "-" + stopCodon.location().bioEnd();
                }
                transcriptSequence.addStopCodonSequence(new AccessionID(stopCodonName), stopCodon.location().bioStart(), stopCodon.location().bioEnd());
            }

            for (FeatureI cdsFeature : cdsFeatures) {
                Feature cds = (Feature) cdsFeature;
                String cdsName = cds.getAttribute("ID");
                if (cdsName == null || cdsName.length() == 0) {
                    cdsName = geneid + "-cds-" + cds.location().bioStart() + "-" + cds.location().bioEnd();
                }
                AccessionID cdsAccessionID = new AccessionID(cdsName);
                ExonSequence exonSequence = geneSequence.addExon(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd());
                transcriptSequence.addCDS(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd(), cds.frame());
            }

        }
        return dnaSequenceList;
    }

    static public LinkedHashMap<String, DNASequence> loadFastaAddGeneFeaturesFromGTF(File fastaSequenceFile, File gffFile) throws Exception {
        LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile);
        FeatureList listGenes = GeneMarkGTF.read(gffFile.getAbsolutePath());
        dnaSequenceList = addGTFGeneFeatures(dnaSequenceList, listGenes);
        return dnaSequenceList;
    }

    static public LinkedHashMap<String, DNASequence> addGTFGeneFeatures(LinkedHashMap<String, DNASequence> dnaSequenceList, FeatureList listGenes) throws Exception {
        Collection<String> geneIds = listGenes.attributeValues("gene_id");
        for (String geneid : geneIds) {
            FeatureList gene = listGenes.selectByAttribute("gene_id", geneid);
            FeatureI geneFeature = gene.get(0);
            DNASequence seq = dnaSequenceList.get(geneFeature.seqname());
            AccessionID geneAccessionID = new AccessionID(geneid);
            GeneSequence geneSequence = null;
            Collection<String> transcriptids = gene.attributeValues("transcript_id");
            for (String transcriptid : transcriptids) {
                // get all the individual features (exons, CDS regions, etc.) of this gene
                FeatureList transcriptFeature = listGenes.selectByAttribute("transcript_id", transcriptid);


                FeatureI feature = transcriptFeature.get(0);
                Strand strand = Strand.POSITIVE;

                if (feature.location().isNegative()) {
                    strand = strand.NEGATIVE;
                }

                String seqName = feature.seqname();
                FeatureI startCodon = null;
                FeatureI stopCodon = null;
                Integer startCodonBegin = null;
                Integer stopCodonEnd = null;
                String startCodonName = "";
                String stopCodonName = "";
                FeatureList startCodonList = transcriptFeature.selectByType("start_codon");
                if (startCodonList != null && startCodonList.size() > 0) {
                    startCodon = startCodonList.get(0);
                    startCodonBegin = startCodon.location().bioStart();
                    startCodonName = startCodon.getAttribute("transcript_name");
                }

                FeatureList stopCodonList = transcriptFeature.selectByType("stop_codon");

                if (stopCodonList != null && stopCodonList.size() > 0) {
                    stopCodon = stopCodonList.get(0);
                    stopCodonEnd = stopCodon.location().bioEnd();
                    stopCodonName = stopCodon.getAttribute("transcript_name");

                }


                // now select only the coding regions of this gene
                FeatureList cdsFeatures = transcriptFeature.selectByType("CDS");
                // sort them
                cdsFeatures = cdsFeatures.sortByStart();

                if (startCodonBegin == null) {
                    FeatureI firstFeature = cdsFeatures.get(0);
                    startCodonBegin = firstFeature.location().bioStart();
                }

                if (stopCodonEnd == null) {
                    FeatureI lastFeature = cdsFeatures.get(cdsFeatures.size() - 1);
                    stopCodonEnd = lastFeature.location().bioEnd();
                }



                AccessionID transcriptAccessionID = new AccessionID(transcriptid);
                if (geneSequence == null) {
                    geneSequence = seq.addGene(geneAccessionID, startCodonBegin, stopCodonEnd);
                } else {
                    //if multiple transcripts for one gene make sure the gene is defined as the min and max start/end
                    if (strand.equals(Strand.POSITIVE)) {
                        if (startCodonBegin < geneSequence.getBioBegin()) {
                            geneSequence.setBioBegin(startCodonBegin);
                        }
                        if (stopCodonEnd > geneSequence.getBioBegin()) {
                            geneSequence.setBioEnd(stopCodonEnd);
                        }
                    } else {
                        if (startCodonBegin > geneSequence.getBioBegin()) {
                            geneSequence.setBioBegin(startCodonBegin);
                        }
                        if (stopCodonEnd < geneSequence.getBioBegin()) {
                            geneSequence.setBioEnd(stopCodonEnd);
                        }
                    }
                }
                TranscriptSequence transcriptSequence = geneSequence.addTranscript(transcriptAccessionID, startCodonBegin, stopCodonEnd, strand);
                if (startCodon != null) {
                    if (startCodonName == null || startCodonName.length() == 0) {
                        startCodonName = transcriptid + "-start_codon-" + startCodon.location().bioStart() + "-" + startCodon.location().bioEnd();
                    }
                    transcriptSequence.addStartCodonSequence(new AccessionID(startCodonName), startCodon.location().bioStart(), startCodon.location().bioEnd());
                }
                if (stopCodon != null) {
                    if (stopCodonName == null || stopCodonName.length() == 0) {
                        stopCodonName = transcriptid + "-stop_codon-" + stopCodon.location().bioStart() + "-" + stopCodon.location().bioEnd();
                    }
                    transcriptSequence.addStopCodonSequence(new AccessionID(stopCodonName), stopCodon.location().bioStart(), stopCodon.location().bioEnd());
                }

                for (FeatureI cdsFeature : cdsFeatures) {
                    Feature cds = (Feature) cdsFeature;
                    String cdsName = cds.getAttribute("transcript_name");
                    if (cdsName == null || cdsName.length() == 0) {
                        cdsName = transcriptid + "-cds-" + cds.location().bioStart() + "-" + cds.location().bioEnd();
                    }
                    AccessionID cdsAccessionID = new AccessionID(cdsName);
                    ExonSequence exonSequence = geneSequence.addExon(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd());
                    transcriptSequence.addCDS(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd(), cds.frame());
                }
            }
        }
        return dnaSequenceList;
    }

    static public LinkedHashMap<String, ProteinSequence> getProteinSequences(Collection<DNASequence> dnaSequences) throws Exception {
        LinkedHashMap<String, ProteinSequence> proteinSequenceHashMap = new LinkedHashMap<String, ProteinSequence>();
        for (DNASequence dnaSequence : dnaSequences) {
            for (GeneSequence geneSequence : dnaSequence.getGeneSequences().values()) {
                for (TranscriptSequence transcriptSequence : geneSequence.getTranscripts().values()) {
                    DNASequence dnaCodingSequence = transcriptSequence.getDNACodingSequence();
                    System.out.println("CDS=" + dnaCodingSequence.getSequenceAsString());
                   
                    try {
                        ProteinSequence proteinSequence = transcriptSequence.getProteinSequence();
                        System.out.println(proteinSequence.getAccession().getID() + " " + proteinSequence);
                        if (proteinSequenceHashMap.containsKey(proteinSequence.getAccession().getID())) {
                            throw new Exception("Duplicate protein sequence id=" + proteinSequence.getAccession().getID() + " found at Gene id=" + geneSequence.getAccession().getID());
                        } else {
                            proteinSequenceHashMap.put(proteinSequence.getAccession().getID(), proteinSequence);
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                }

            }
        }
        return proteinSequenceHashMap;
    }

    public static void main(String args[]) throws Exception {
        if (false) {
            LinkedHashMap<String, DNASequence> dnaSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGTF(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds.fna"), new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/genemark_hmm.gtf"));
            LinkedHashMap<String, ProteinSequence> proteinSequenceList = GeneFeatureHelper.getProteinSequences(dnaSequenceList.values());
        }

        if (true) {
            LinkedHashMap<String, DNASequence> dnaSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGFF3(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/454Scaffolds.fna"), new File("/Users/Scooter/scripps/dyadic/GlimmerHMM/c1_glimmerhmm.gff"));
            LinkedHashMap<String, ProteinSequence> proteinSequenceList = GeneFeatureHelper.getProteinSequences(dnaSequenceList.values());
          //  for (ProteinSequence proteinSequence : proteinSequenceList.values()) {
          //      System.out.println(proteinSequence.getAccession().getID() + " " + proteinSequence);
          //  }
            FastaWriterHelper.writeProteinSequence(new File("/Users/Scooter/scripps/dyadic/GlimmerHMM/c1_predicted_glimmer.faa"), proteinSequenceList.values());

        }



    }
}
