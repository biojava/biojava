package org.biojava.nbio.genome;

import org.biojava.nbio.core.sequence.*;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GeneIDGFF2Reader;

import java.io.File;
import java.util.Collection;
import java.util.LinkedHashMap;

/**
 * @author Elizabeth James
 */
public class GeneGFF2FeatureHelper {

    /**
     * Loads Fasta file and GFF2 feature file generated from the geneid prediction algorithm
     *
     * @param fastaSequenceFile
     * @param gffFile
     * @return
     * @throws Exception
     */
    static public LinkedHashMap<String, ChromosomeSequence> loadFastaAddGeneFeaturesFromGeneIDGFF2(File fastaSequenceFile, File gffFile) throws Exception {
        LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile);
        LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneFeatureHelper.getChromosomeSequenceFromDNASequence(dnaSequenceList);
        FeatureList listGenes = GeneIDGFF2Reader.read(gffFile.getAbsolutePath());
        addGeneIDGFF2GeneFeatures(chromosomeSequenceList, listGenes);
        return chromosomeSequenceList;
    }
    /**
     * Load GFF2 feature file generated from the geneid prediction algorithm and map features onto the chromosome sequences
     *
     * @param chromosomeSequenceList
     * @param listGenes
     * @throws Exception
     */
    static public void addGeneIDGFF2GeneFeatures(LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList, FeatureList listGenes) throws Exception {
        Collection<String> geneIds = listGenes.attributeValues("gene_id");
        for (String geneid : geneIds) {
            FeatureList gene = listGenes.selectByAttribute("gene_id", geneid);
            FeatureI geneFeature = gene.get(0);
            ChromosomeSequence seq = chromosomeSequenceList.get(geneFeature.seqname());
            geneid = geneid.replaceAll("_", ".G");
            AccessionID geneAccessionID = new AccessionID(geneid);
            GeneSequence geneSequence = null;
            Collection<String> transcriptids = gene.attributeValues("gene_id");
            for (String transcriptid : transcriptids) {
                // get all the individual features (exons, CDS regions, etc.) of this gene
                FeatureList transcriptFeature = listGenes.selectByAttribute("gene_id", transcriptid);
                transcriptid = transcriptid.replaceAll("_", ".G");




                //      String seqName = feature.seqname();
                //FeatureI startCodon = null;
                //FeatureI stopCodon = null;
                Integer startCodonBegin = null;
                Integer stopCodonEnd = null;
                //String startCodonName = "";
                //String stopCodonName = "";


                // now select only the coding regions of this gene
                FeatureList firstFeatures = transcriptFeature.selectByType("First");
                FeatureList terminalFeatures = transcriptFeature.selectByType("Terminal");
                FeatureList internalFeatures = transcriptFeature.selectByType("Internal");
                FeatureList singleFeatures = transcriptFeature.selectByType("Single");
                FeatureList cdsFeatures = new FeatureList();
                cdsFeatures.add(firstFeatures);
                cdsFeatures.add(terminalFeatures);
                cdsFeatures.add(internalFeatures);
                cdsFeatures.add(singleFeatures);
                // sort them
                cdsFeatures = cdsFeatures.sortByStart();
                Strand strand = Strand.POSITIVE;
                FeatureI feature = cdsFeatures.get(0);
                if (feature.location().isNegative()) {
                    strand = Strand.NEGATIVE;
                }
                if (startCodonBegin == null) {
                    FeatureI firstFeature = cdsFeatures.get(0);
                    if (strand == Strand.NEGATIVE) {
                        startCodonBegin = firstFeature.location().bioEnd();
                    } else {
                        startCodonBegin = firstFeature.location().bioStart();
                    }
                }

                if (stopCodonEnd == null) {

                    FeatureI lastFeature = cdsFeatures.get(cdsFeatures.size() - 1);
                    if (strand == Strand.NEGATIVE) {
                        stopCodonEnd = lastFeature.location().bioStart();
                    } else {
                        stopCodonEnd = lastFeature.location().bioEnd();
                    }
                }
                //for gtf ordering can be strand based so first is last and last is first
                if (startCodonBegin > stopCodonEnd) {
                    int temp = startCodonBegin;
                    startCodonBegin = stopCodonEnd;
                    stopCodonEnd = temp;
                }

                AccessionID transcriptAccessionID = new AccessionID(transcriptid);
                if (geneSequence == null) {
                    geneSequence = seq.addGene(geneAccessionID, startCodonBegin, stopCodonEnd, strand);
                    geneSequence.setSource(((Feature) feature).source());
                } else {
                    //if multiple transcripts for one gene make sure the gene is defined as the min and max start/end

                    if (startCodonBegin < geneSequence.getBioBegin()) {
                        geneSequence.setBioBegin(startCodonBegin);
                    }
                    if (stopCodonEnd > geneSequence.getBioBegin()) {
                        geneSequence.setBioEnd(stopCodonEnd);
                    }

                }
                TranscriptSequence transcriptSequence = geneSequence.addTranscript(transcriptAccessionID, startCodonBegin, stopCodonEnd);
				/*
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
				*/

                for (FeatureI cdsFeature : cdsFeatures) {
                    Feature cds = (Feature) cdsFeature;
                    String cdsName = cds.getAttribute("transcript_name");
                    if (cdsName == null || cdsName.length() == 0) {
                        cdsName = transcriptid + "-cds-" + cds.location().bioStart() + "-" + cds.location().bioEnd();
                    }
                    AccessionID cdsAccessionID = new AccessionID(cdsName);
                    //ExonSequence exonSequence = geneSequence.addExon(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd());
                    CDSSequence cdsSequence = transcriptSequence.addCDS(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd(), cds.frame());
                    cdsSequence.setSequenceScore(cds.score());
                }
            }
        }

    }

}
