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
package org.biojava.nbio.genome.query;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GeneMarkGTFReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class OutputHitsGFF {

	private static final Logger logger = LoggerFactory.getLogger(OutputHitsGFF.class);

	public void process(File blastXMLFile, File gffFile, File gffOutputFile, double maxEScore, double percentageAligned, boolean includeFrameShift, boolean includeNegativeStrand) throws Exception {
		BlastXMLQuery blastXMLQuery = new BlastXMLQuery(blastXMLFile.getAbsolutePath());
		LinkedHashMap<String, ArrayList<String>> hits = blastXMLQuery.getHitsQueryDef(maxEScore);
		FeatureList listGenes = GeneMarkGTFReader.read(gffFile.getAbsolutePath());
		FeatureList hitGenes = new FeatureList();
		for (String id : hits.keySet()) {
			String[] values = id.split(" ");
			String gene_id = values[0];
			FeatureList gene = listGenes.selectByAttribute("gene_id", gene_id);
			for (FeatureI geneFeature : gene) {

				if (!includeNegativeStrand && geneFeature.location().isNegative()) {
					continue;
				}
				if (!includeFrameShift) {
					boolean frameShift = false;
					FeatureList cdsList = gene.selectByType("CDS");
					for(FeatureI cdsFeature : cdsList){
						int frame = ((Feature)cdsFeature).frame();
						if(frame != 0){
							frameShift = true;
							break;
						}
					}
					if(frameShift)
						continue;
				}
				hitGenes.add(geneFeature);
			}
		}

	//    GeneMarkGTFReader.write(hitGenes, gffOutputFile.getAbsolutePath());
	}


		public static void main(String[] args) {
		try {
			OutputHitsGFF outputHitsGFF = new OutputHitsGFF();
			outputHitsGFF.process(new File("hits-uniprot_fungi.xml"),
					new File("genemark_hmm.gtf"),
					new File("genemark_hits_hmm.gtf"), 0, 100, true, true);


		} catch (Exception e) {
			logger.error("Execution: ", e);
		}
	}
}
