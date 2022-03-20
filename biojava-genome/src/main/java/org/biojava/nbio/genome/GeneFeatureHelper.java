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
package org.biojava.nbio.genome;

import org.biojava.nbio.genome.parsers.gff.*;
import org.biojava.nbio.core.sequence.*;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GeneFeatureHelper {

	private static final Logger logger = LoggerFactory.getLogger(GeneFeatureHelper.class);

	static public LinkedHashMap<String, ChromosomeSequence> loadFastaAddGeneFeaturesFromUpperCaseExonFastaFile(File fastaSequenceFile, File uppercaseFastaFile, boolean throwExceptionGeneNotFound) throws Exception {
		LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = new LinkedHashMap<String, ChromosomeSequence>();
		LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile);
		for (String accession : dnaSequenceList.keySet()) {
			DNASequence contigSequence = dnaSequenceList.get(accession);
			ChromosomeSequence chromsomeSequence = new ChromosomeSequence(contigSequence.getSequenceAsString());
			chromsomeSequence.setAccession(contigSequence.getAccession());
			chromosomeSequenceList.put(accession, chromsomeSequence);
		}


		LinkedHashMap<String, DNASequence> geneSequenceList = FastaReaderHelper.readFastaDNASequence(uppercaseFastaFile);
		for (DNASequence dnaSequence : geneSequenceList.values()) {
			String geneSequence = dnaSequence.getSequenceAsString();
			String lcGeneSequence = geneSequence.toLowerCase();
			String reverseGeneSequence = dnaSequence.getReverse().getSequenceAsString();
			String lcReverseGeneSequence = reverseGeneSequence.toLowerCase();
			Integer bioStart = null;
			Integer bioEnd = null;
			Strand strand = Strand.POSITIVE;
			boolean geneFound = false;
			String accession = "";
			DNASequence contigDNASequence = null;
			for (String id : dnaSequenceList.keySet()) {
				accession = id;
				contigDNASequence = dnaSequenceList.get(id);
				String contigSequence = contigDNASequence.getSequenceAsString().toLowerCase();
				bioStart = contigSequence.indexOf(lcGeneSequence);
				if (bioStart != -1) {
					bioStart = bioStart + 1;
					bioEnd = bioStart + geneSequence.length() - 1;
					geneFound = true;
					break;
				} else {
					bioStart = contigSequence.indexOf(lcReverseGeneSequence);
					if (bioStart != -1) {
						bioStart = bioStart + 1;
						bioEnd = bioStart - geneSequence.length() - 1;
						strand = Strand.NEGATIVE;
						geneFound = true;
						break;
					}
				}
			}

			if (geneFound) {
				logger.info("Gene {} found at {} {} {} {}",
						dnaSequence.getAccession().toString(), contigDNASequence.getAccession().toString(), bioStart, bioEnd, strand);
				ChromosomeSequence chromosomeSequence = chromosomeSequenceList.get(accession);

				ArrayList<Integer> exonBoundries = new ArrayList<Integer>();

				//look for transitions from lowercase to upper case
				for (int i = 0; i < geneSequence.length(); i++) {
					if (i == 0 && Character.isUpperCase(geneSequence.charAt(i))) {
						exonBoundries.add(i);
					} else if (i == geneSequence.length() - 1) {
						exonBoundries.add(i);
					} else if (Character.isUpperCase(geneSequence.charAt(i)) && Character.isLowerCase(geneSequence.charAt(i - 1))) {
						exonBoundries.add(i);
					} else if (Character.isUpperCase(geneSequence.charAt(i)) && Character.isLowerCase(geneSequence.charAt(i + 1))) {
						exonBoundries.add(i);
					}
				}
				if (strand == Strand.NEGATIVE) {
					Collections.reverse(exonBoundries);
				}


				String geneaccession = dnaSequence.getAccession().getID();
				String note = geneaccession;
				String[] values = geneaccession.split(" ");
				geneaccession = values[0];



				GeneSequence geneSeq = chromosomeSequence.addGene(new AccessionID(geneaccession), bioStart, bioEnd, strand);
				geneSeq.addNote(note);
				geneSeq.setSource(uppercaseFastaFile.getName());
				//String transcriptName = geneaccession + "-transcript";
				//TranscriptSequence transcriptSequence = geneSeq.addTranscript(new AccessionID(transcriptName), bioStart, bioEnd);

				int runningFrameLength = 0;
				for (int i = 0; i < exonBoundries.size() - 1; i = i + 2) {
					int cdsBioStart = exonBoundries.get(i) + bioStart;
					int cdsBioEnd = exonBoundries.get(i + 1) + bioStart;
					runningFrameLength = runningFrameLength + Math.abs(cdsBioEnd - cdsBioStart) + 1;
					//String cdsName = transcriptName + "-cds-" + cdsBioStart + "-" + cdsBioEnd;

					//AccessionID cdsAccessionID = new AccessionID(cdsName);
					//ExonSequence exonSequence = geneSeq.addExon(cdsAccessionID, cdsBioStart, cdsBioEnd);
					int remainder = runningFrameLength % 3;
					//int frame = 0;
					if (remainder == 1) {
						//frame = 2; // borrow 2 from next CDS region
					} else if (remainder == 2) {
						//frame = 1;
					}
					//CDSSequence cdsSequence = transcriptSequence.addCDS(cdsAccessionID, cdsBioStart, cdsBioEnd, frame);
				}


			} else {
				if (throwExceptionGeneNotFound) {
					throw new Exception(dnaSequence.getAccession().toString() + " not found");
				}
				logger.info("Gene not found {}", dnaSequence.getAccession().toString());
			}

		}
		return chromosomeSequenceList;
	}



	static public LinkedHashMap<String, ChromosomeSequence> getChromosomeSequenceFromDNASequence(LinkedHashMap<String, DNASequence> dnaSequenceList) {
		LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = new LinkedHashMap<String, ChromosomeSequence>();
		for (String key : dnaSequenceList.keySet()) {
			DNASequence dnaSequence = dnaSequenceList.get(key);
			ChromosomeSequence chromosomeSequence = new ChromosomeSequence(dnaSequence.getProxySequenceReader()); //we want the underlying sequence but don't need storage
			chromosomeSequence.setAccession(dnaSequence.getAccession());
			chromosomeSequenceList.put(key, chromosomeSequence);
		}
		return chromosomeSequenceList;
	}


	static public LinkedHashMap<String, ChromosomeSequence> loadFastaAddGeneFeaturesFromGeneMarkGTF(File fastaSequenceFile, File gffFile) throws Exception {
		LinkedHashMap<String, DNASequence> dnaSequenceList = FastaReaderHelper.readFastaDNASequence(fastaSequenceFile);
		LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneFeatureHelper.getChromosomeSequenceFromDNASequence(dnaSequenceList);
		FeatureList listGenes = GeneMarkGTFReader.read(gffFile.getAbsolutePath());
		addGeneMarkGTFGeneFeatures(chromosomeSequenceList, listGenes);
		return chromosomeSequenceList;
	}

	static public void addGeneMarkGTFGeneFeatures(LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList, FeatureList listGenes) throws Exception {
		Collection<String> geneIds = listGenes.attributeValues("gene_id");
		for (String geneid : geneIds) {
			//       if(geneid.equals("45_g")){
			//           int dummy =1;
			//       }
			FeatureList gene = listGenes.selectByAttribute("gene_id", geneid);
			FeatureI geneFeature = gene.get(0);
			ChromosomeSequence seq = chromosomeSequenceList.get(geneFeature.seqname());
			AccessionID geneAccessionID = new AccessionID(geneid);
			GeneSequence geneSequence = null;
			Collection<String> transcriptids = gene.attributeValues("transcript_id");
			for (String transcriptid : transcriptids) {
				// get all the individual features (exons, CDS regions, etc.) of this gene


				FeatureList transcriptFeature = listGenes.selectByAttribute("transcript_id", transcriptid);
				// now select only the coding regions of this gene
				FeatureList cdsFeatures = transcriptFeature.selectByType("CDS");
				// sort them
				cdsFeatures = cdsFeatures.sortByStart();

				FeatureI feature = cdsFeatures.get(0);
				Strand strand = Strand.POSITIVE;

				if (feature.location().isNegative()) {
					strand = Strand.NEGATIVE;
				}

				//String seqName = feature.seqname();
				FeatureI startCodon = null;
				FeatureI stopCodon = null;
				Integer startCodonBegin = null;
				Integer stopCodonEnd = null;
				String startCodonName = "";
				String stopCodonName = "";
				FeatureList startCodonList = transcriptFeature.selectByType("start_codon");
				if (startCodonList != null && startCodonList.size() > 0) {
					startCodon = startCodonList.get(0);
					if (strand == Strand.POSITIVE) {
						startCodonBegin = startCodon.location().bioStart();
					} else {
						startCodonBegin = startCodon.location().bioEnd();
					}
					startCodonName = startCodon.getAttribute("transcript_name");
				}

				FeatureList stopCodonList = transcriptFeature.selectByType("stop_codon");

				if (stopCodonList != null && stopCodonList.size() > 0) {
					stopCodon = stopCodonList.get(0);
					if (strand == Strand.POSITIVE) {
						stopCodonEnd = stopCodon.location().bioEnd();
					} else {
						stopCodonEnd = stopCodon.location().bioStart();
					}

					stopCodonName = stopCodon.getAttribute("transcript_name");

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
					// for genemark it appears frame of 2 =1 and frame of 1 = 2
					// doesn't matter when you string cds regions together as one block
					// but does make a difference when you try to make a protein sequence for each CDS region where
					// you give up or borrow based on the frame value
					// compared with gff like files and docs for geneid and glimmer where geneid and glimmer both do it the same
					// way that appears to match the gff3 docs.
					int frame = cds.frame();
					if (frame == 1) {
						frame = 2;
					} else if (frame == 2) {
						frame = 1;
					} else {
						frame = 0;
					}
					String cdsName = cds.getAttribute("transcript_name");
					if (cdsName == null || cdsName.length() == 0) {
						cdsName = transcriptid + "-cds-" + cds.location().bioStart() + "-" + cds.location().bioEnd();
					}
					AccessionID cdsAccessionID = new AccessionID(cdsName);
					//ExonSequence exonSequence = geneSequence.addExon(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd());
					transcriptSequence.addCDS(cdsAccessionID, cdsFeature.location().bioStart(), cdsFeature.location().bioEnd(), frame);
				}
			}
		}

	}

	static public LinkedHashMap<String, ProteinSequence> getProteinSequences(Collection<ChromosomeSequence> chromosomeSequences) throws Exception {
		LinkedHashMap<String, ProteinSequence> proteinSequenceHashMap = new LinkedHashMap<String, ProteinSequence>();
		for (ChromosomeSequence dnaSequence : chromosomeSequences) {
			for (GeneSequence geneSequence : dnaSequence.getGeneSequences().values()) {
				for (TranscriptSequence transcriptSequence : geneSequence.getTranscripts().values()) {
					//TODO remove?
//                    DNASequence dnaCodingSequence = transcriptSequence.getDNACodingSequence();
//                    logger.info("CDS={}", dnaCodingSequence.getSequenceAsString());

					try {
						ProteinSequence proteinSequence = transcriptSequence.getProteinSequence();

//                        logger.info("{} {}", proteinSequence.getAccession().getID(), proteinSequence);
						if (proteinSequenceHashMap.containsKey(proteinSequence.getAccession().getID())) {
							throw new Exception("Duplicate protein sequence id=" + proteinSequence.getAccession().getID() + " found at Gene id=" + geneSequence.getAccession().getID());
						} else {
							proteinSequenceHashMap.put(proteinSequence.getAccession().getID(), proteinSequence);
						}
					} catch (Exception e) {
						logger.error("Exception: ", e);
					}

				}

			}
		}
		return proteinSequenceHashMap;
	}

	static public LinkedHashMap<String, GeneSequence> getGeneSequences(Collection<ChromosomeSequence> chromosomeSequences) throws Exception {
		LinkedHashMap<String, GeneSequence> geneSequenceHashMap = new LinkedHashMap<String, GeneSequence>();
		for (ChromosomeSequence chromosomeSequence : chromosomeSequences) {
			for (GeneSequence geneSequence : chromosomeSequence.getGeneSequences().values()) {
				geneSequenceHashMap.put(geneSequence.getAccession().getID(), geneSequence);
			}
		}

		return geneSequenceHashMap;
	}

	public static void main(String[] args) throws Exception {
/*        if (false) {
			LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGeneMarkGTF(new File("Scaffolds.fna"), new File("genemark_hmm.gtf"));
			LinkedHashMap<String, ProteinSequence> proteinSequenceList = GeneFeatureHelper.getProteinSequences(chromosomeSequenceList.values());
		}

		if (false) {
			LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceList = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGlimmerGFF3(new File("Scaffolds.fna"), new File("glimmerhmm.gff"));
			LinkedHashMap<String, ProteinSequence> proteinSequenceList = GeneFeatureHelper.getProteinSequences(chromosomeSequenceList.values());
			//  for (ProteinSequence proteinSequence : proteinSequenceList.values()) {
			//      logger.info(proteinSequence.getAccession().getID() + " " + proteinSequence);
			//  }
			FastaWriterHelper.writeProteinSequence(new File("predicted_glimmer.faa"), proteinSequenceList.values());

		}
		if (false) {
			GeneFeatureHelper.outputFastaSequenceLengthGFF3(new File("Scaffolds.fna"), new File("scaffolds.gff3"));
		}

 */

		try {

			if (true) {
				//File fastaSequenceFile = new File("Scaffolds.fna");
				//File gff3File = new File("geneid-v6.gff3");
				//LinkedHashMap<String, ChromosomeSequence> chromosomeSequenceHashMap = GeneFeatureHelper.loadFastaAddGeneFeaturesFromGmodGFF3(fastaSequenceFile, gff3File,true);
			}

		} catch (Exception e) {
			logger.error("Exception: ", e);
		}

	}
}
