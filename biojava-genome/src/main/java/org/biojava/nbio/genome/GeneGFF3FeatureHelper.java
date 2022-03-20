package org.biojava.nbio.genome;

import org.biojava.nbio.core.sequence.*;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.genome.parsers.gff.*;

import java.io.File;
import java.io.FileWriter;
import java.util.LinkedHashMap;

/**
 * @author Elizabeth James
 */
public class GeneGFF3FeatureHelper {

    /**
     * Lots of variations in the ontology or descriptors that can be used in GFF3 which requires writing a custom parser to handle a GFF3 generated or used
     * by a specific application. Probably could be abstracted out but for now easier to handle with custom code to deal with gff3 elements that are not
     * included but can be extracted from other data elements.
     * @param fastaSequenceFile
     * @param gffFile
     * @param lazyloadsequences If set to true then the fasta file will be parsed for accession id but sequences will be read from disk when needed to save memory
     * @return
     * @throws Exception
     */
    static public LinkedHashMap<String, ChromosomeSequence> loadFastaAddGeneFeaturesFromGmodGFF3(File fastaSequenceFile, File gffFile, boolean lazyloadsequences) throws Exception {
        LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile,lazyloadsequences);
        LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneFeatureHelper.getChromosomeSequenceFromDNASequence(dnaSequenceList);
        FeatureList listGenes = GFF3Reader.read(gffFile.getAbsolutePath());
        addGmodGFF3GeneFeatures(chromosomeSequenceList, listGenes);
        return chromosomeSequenceList;
    }

    /**
     * Load GFF3 file using mRNA as the gene feature as not all GFF3 files are complete
     * @param chromosomeSequenceList
     * @param listGenes
     * @throws Exception
     */
    static public void addGmodGFF3GeneFeatures(LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList, FeatureList listGenes) throws Exception {


        // key off mRNA as being a known feature that may or may not have a parent gene


        FeatureList mRNAFeatures = listGenes.selectByType("mRNA");
        LinkedHashMap<String,FeatureList> featureIDHashMap = FeatureHelper.buildFeatureAtrributeIndex("ID", listGenes);
        LinkedHashMap<String,FeatureList> featureParentHashMap = FeatureHelper.buildFeatureAtrributeIndex("Parent", listGenes);

        for (FeatureI f : mRNAFeatures) {
            String geneID;
            String geneNote = null;
            String geneSource = null;
            String sequenceName = null;
            ChromosomeSequence seq = null;
            GeneSequence geneSequence = null;

            Feature mRNAFeature = (Feature) f;
            String mRNAID = mRNAFeature.getAttribute("ID");
            String mRNAsource = mRNAFeature.source();
            String mRNANote = mRNAFeature.getAttribute("Note");
            String mRNAParent = mRNAFeature.getAttribute("Parent");
            if (mRNAParent != null && mRNAParent.length() > 0) {
                // FeatureList geneFeatureList = listGenes.selectByAttribute("ID", mRNAParent);
                FeatureList geneFeatureList = featureIDHashMap.get(mRNAParent);
                Feature geneFeature = (Feature) geneFeatureList.get(0);
                geneID = geneFeature.getAttribute("ID");
                geneNote = geneFeature.getAttribute("Note");
                geneSource = geneFeature.source();
                sequenceName = geneFeature.seqname();

                //
            } else {
                //deal with cases where no parent gene is given
                geneID = mRNAID;
                geneSource = mRNAsource;
                sequenceName = mRNAFeature.seqname();
            }

            seq = chromosomeSequenceList.get(sequenceName);

            AccessionID geneAccessionID = new AccessionID(geneID);

            //  FeatureList mRNAChildren = listGenes.selectByAttribute("Parent", mRNAID);
            FeatureList mRNAChildren = featureParentHashMap.get(mRNAID);
            FeatureList cdsFeatures = mRNAChildren.selectByType("CDS");
            FeatureI feature = cdsFeatures.get(0);
            Strand strand = Strand.POSITIVE;

            if (feature.location().isNegative()) {
                strand = Strand.NEGATIVE;
            }
            cdsFeatures = cdsFeatures.sortByStart();







            //String seqName = feature.seqname();
            FeatureI startCodon = null;
            FeatureI stopCodon = null;
            Integer startCodonBegin = null;
            Integer stopCodonEnd = null;
            String startCodonName = "";
            String stopCodonName = "";
            FeatureList startCodonList = mRNAChildren.selectByType("five_prime_UTR");
            if (startCodonList != null && startCodonList.size() > 0) {
                startCodon = startCodonList.get(0);
                if (strand == Strand.NEGATIVE) {
                    startCodonBegin = startCodon.location().bioEnd();
                } else {
                    startCodonBegin = startCodon.location().bioStart();
                }
                startCodonName = startCodon.getAttribute("ID");
            }

            FeatureList stopCodonList = mRNAChildren.selectByType("three_prime_UTR");

            if (stopCodonList != null && stopCodonList.size() > 0) {
                stopCodon = stopCodonList.get(0);
                if (strand == Strand.NEGATIVE) {
                    stopCodonEnd = stopCodon.location().bioStart();
                } else {
                    stopCodonEnd = stopCodon.location().bioEnd();
                }
                stopCodonName = stopCodon.getAttribute("ID");

            }




            if (startCodonBegin == null) {
                if (strand == Strand.NEGATIVE) {
                    FeatureI firstFeature = cdsFeatures.get(0);

                    startCodonBegin = firstFeature.location().bioEnd();
                } else {
                    FeatureI firstFeature = cdsFeatures.get(0);

                    startCodonBegin = firstFeature.location().bioStart();
                }
            }

            if (stopCodonEnd == null) {
                if (strand == Strand.NEGATIVE) {
                    FeatureI lastFeature = cdsFeatures.get(cdsFeatures.size() - 1);
                    stopCodonEnd = lastFeature.location().bioStart();
                } else {
                    FeatureI lastFeature = cdsFeatures.get(cdsFeatures.size() - 1);
                    stopCodonEnd = lastFeature.location().bioEnd();
                }
            }
            //for gtf ordering can be strand based so first is last and last is first
            if (startCodonBegin > stopCodonEnd) {
                int temp = startCodonBegin;
                startCodonBegin = stopCodonEnd;
                stopCodonEnd = temp;
            }



            AccessionID transcriptAccessionID = new AccessionID(mRNAID);
            geneSequence = seq.getGene(geneID);
            if (geneSequence == null) {
                geneSequence = seq.addGene(geneAccessionID, startCodonBegin, stopCodonEnd, strand);
                geneSequence.setSource(geneSource);
                if (geneNote != null && geneNote.length() > 0) {
                    geneSequence.addNote(geneNote);
                }
            } else {

                if (startCodonBegin < geneSequence.getBioBegin()) {
                    geneSequence.setBioBegin(startCodonBegin);
                }
                if (stopCodonEnd > geneSequence.getBioBegin()) {
                    geneSequence.setBioEnd(stopCodonEnd);
                }

            }
            TranscriptSequence transcriptSequence = geneSequence.addTranscript(transcriptAccessionID, startCodonBegin, stopCodonEnd);
            transcriptSequence.setSource(mRNAsource);
            if (mRNANote != null && mRNANote.length() > 0) {
                transcriptSequence.addNote(mRNANote);

            }
            if (startCodon != null) {
                if (startCodonName == null || startCodonName.length() == 0) {
                    startCodonName = mRNAID + "-start_codon-" + startCodon.location().bioStart() + "-" + startCodon.location().bioEnd();
                }
                transcriptSequence.addStartCodonSequence(new AccessionID(startCodonName), startCodon.location().bioStart(), startCodon.location().bioEnd());
            }
            if (stopCodon != null) {
                if (stopCodonName == null || stopCodonName.length() == 0) {
                    stopCodonName = mRNAID + "-stop_codon-" + stopCodon.location().bioStart() + "-" + stopCodon.location().bioEnd();
                }
                transcriptSequence.addStopCodonSequence(new AccessionID(stopCodonName), stopCodon.location().bioStart(), stopCodon.location().bioEnd());
            }

            for (FeatureI cdsFeature : cdsFeatures) {
                Feature cds = (Feature) cdsFeature;
                String cdsNote = cdsFeature.getAttribute("Note");
                String cdsSource = cds.source();
                String cdsName = cds.getAttribute("ID");
                if (cdsName == null || cdsName.length() == 0) {
                    cdsName = mRNAID + "-cds-" + cds.location().bioStart() + "-" + cds.location().bioEnd();
                }
                AccessionID cdsAccessionID = new AccessionID(cdsName);
                ExonSequence exonSequence = geneSequence.addExon(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd());
                exonSequence.setSource(cdsSource);
                if (cdsNote != null && cdsNote.length() > 0) {
                    exonSequence.addNote(cdsNote);
                }
                transcriptSequence.addCDS(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd(), cds.frame());
            }
            geneSequence.addIntronsUsingExons();

        }

    }

    static public void addGlimmerGFF3GeneFeatures(LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList, FeatureList listGenes) throws Exception {
        FeatureList mRNAFeatures = listGenes.selectByType("mRNA");
        for (FeatureI f : mRNAFeatures) {
            Feature mRNAFeature = (Feature) f;
            String geneid = mRNAFeature.getAttribute("ID");
            String source = mRNAFeature.source();

            FeatureList gene = listGenes.selectByAttribute("Parent", geneid);
            FeatureI geneFeature = gene.get(0);
            ChromosomeSequence seq = chromosomeSequenceList.get(geneFeature.seqname());
            AccessionID geneAccessionID = new AccessionID(geneid);
            GeneSequence geneSequence = null;

            FeatureList cdsFeatures = gene.selectByType("CDS");
            FeatureI feature = cdsFeatures.get(0);
            Strand strand = Strand.POSITIVE;

            if (feature.location().isNegative()) {
                strand = Strand.NEGATIVE;
            }
            cdsFeatures = cdsFeatures.sortByStart();







            //String seqName = feature.seqname();
            FeatureI startCodon = null;
            FeatureI stopCodon = null;
            Integer startCodonBegin = null;
            Integer stopCodonEnd = null;
            String startCodonName = "";
            String stopCodonName = "";
            FeatureList startCodonList = gene.selectByAttribute("Note", "initial-exon");
            if (startCodonList != null && startCodonList.size() > 0) {
                startCodon = startCodonList.get(0);
                if (strand == Strand.NEGATIVE) {
                    startCodonBegin = startCodon.location().bioEnd();
                } else {
                    startCodonBegin = startCodon.location().bioStart();
                }
                startCodonName = startCodon.getAttribute("ID");
            }

            FeatureList stopCodonList = gene.selectByAttribute("Note", "final-exon");

            if (stopCodonList != null && stopCodonList.size() > 0) {
                stopCodon = stopCodonList.get(0);
                if (strand == Strand.NEGATIVE) {
                    stopCodonEnd = stopCodon.location().bioStart();
                } else {
                    stopCodonEnd = stopCodon.location().bioEnd();
                }
                stopCodonName = stopCodon.getAttribute("ID");

            }




            if (startCodonBegin == null) {
                if (strand == Strand.NEGATIVE) {
                    FeatureI firstFeature = cdsFeatures.get(0);

                    startCodonBegin = firstFeature.location().bioEnd();
                } else {
                    FeatureI firstFeature = cdsFeatures.get(0);

                    startCodonBegin = firstFeature.location().bioStart();
                }
            }

            if (stopCodonEnd == null) {
                if (strand == Strand.NEGATIVE) {
                    FeatureI lastFeature = cdsFeatures.get(cdsFeatures.size() - 1);
                    stopCodonEnd = lastFeature.location().bioStart();
                } else {
                    FeatureI lastFeature = cdsFeatures.get(cdsFeatures.size() - 1);
                    stopCodonEnd = lastFeature.location().bioEnd();
                }
            }
            //for gtf ordering can be strand based so first is last and last is first
            if (startCodonBegin > stopCodonEnd) {
                int temp = startCodonBegin;
                startCodonBegin = stopCodonEnd;
                stopCodonEnd = temp;
            }



            AccessionID transcriptAccessionID = new AccessionID(geneid);
            if (geneSequence == null) {
                geneSequence = seq.addGene(geneAccessionID, startCodonBegin, stopCodonEnd, strand);
                geneSequence.setSource(source);
            }
			/*
			else {

				if (startCodonBegin < geneSequence.getBioBegin()) {
					geneSequence.setBioBegin(startCodonBegin);
				}
				if (stopCodonEnd > geneSequence.getBioBegin()) {
					geneSequence.setBioEnd(stopCodonEnd);
				}
			}
			*/
            TranscriptSequence transcriptSequence = geneSequence.addTranscript(transcriptAccessionID, startCodonBegin, stopCodonEnd);
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
                //ExonSequence exonSequence = geneSequence.addExon(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd());
                transcriptSequence.addCDS(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd(), cds.frame());
            }

        }

    }

    static public LinkedHashMap<String, ChromosomeSequence> loadFastaAddGeneFeaturesFromGlimmerGFF3(File fastaSequenceFile, File gffFile) throws Exception {
        LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile);
        LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneFeatureHelper.getChromosomeSequenceFromDNASequence(dnaSequenceList);
        FeatureList listGenes = GFF3Reader.read(gffFile.getAbsolutePath());
        addGlimmerGFF3GeneFeatures(chromosomeSequenceList, listGenes);
        return chromosomeSequenceList;
    }
    /**
     * Output a gff3 feature file that will give the length of each scaffold/chromosome in the fasta file.
     * Used for gbrowse so it knows length.
     * @param fastaSequenceFile
     * @param gffFile
     * @throws Exception
     */
    static public void outputFastaSequenceLengthGFF3(File fastaSequenceFile, File gffFile) throws Exception {
        LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile);
        String fileName = fastaSequenceFile.getName();
        FileWriter fw = new FileWriter(gffFile);
        String newLine = System.getProperty("line.separator");
        fw.write("##gff-version 3" + newLine);
        for (DNASequence dnaSequence : dnaSequenceList.values()) {
            String gff3line = dnaSequence.getAccession().getID() + "\t" + fileName + "\t" + "contig" + "\t" + "1" + "\t" + dnaSequence.getBioEnd() + "\t.\t.\t.\tName=" + dnaSequence.getAccession().getID() + newLine;
            fw.write(gff3line);
        }
        fw.close();
    }

}
