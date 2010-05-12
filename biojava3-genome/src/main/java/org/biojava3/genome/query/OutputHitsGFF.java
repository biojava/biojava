/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.genome.query;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import org.biojava3.genome.parsers.gff.Feature;
import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GeneMarkGTF;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class OutputHitsGFF {

    public void process(File blastXMLFile, File gffFile, File gffOutputFile, double maxEScore, double percentageAligned, boolean includeFrameShift, boolean includeNegativeStrand) throws Exception {
        BlastXMLQuery blastXMLQuery = new BlastXMLQuery(blastXMLFile.getAbsolutePath());
        LinkedHashMap<String, ArrayList<String>> hits = blastXMLQuery.getHitsQueryDef(maxEScore);
        FeatureList listGenes = GeneMarkGTF.read(gffFile.getAbsolutePath());
        FeatureList hitGenes = new FeatureList();
        for (String id : hits.keySet()) {
            String[] values = id.split(" ");
            String gene_id = values[0];
            FeatureList gene = listGenes.selectByAttribute("gene_id", gene_id);
            for (FeatureI geneFeature : gene) {

                if (includeNegativeStrand == false && geneFeature.location().isNegative()) {
                    continue;
                }
                if (includeFrameShift == false) {
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

        GeneMarkGTF.write(hitGenes, gffOutputFile.getAbsolutePath());
    }


        public static void main(String[] args) {
        try {
            OutputHitsGFF outputHitsGFF = new OutputHitsGFF();
            outputHitsGFF.process(new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/c1-454Scaffolds-hits-uniprot_fungi.xml"),
                    new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/genemark_hmm.gtf"),
                    new File("/Users/Scooter/scripps/dyadic/analysis/454Scaffolds/genemark_hits_hmm.gtf"), 0, 100, true, true);


        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
