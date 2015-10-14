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
package org.biojava.nbio.structure.symmetry.core;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;

import java.util.*;

/**
 * Represents a set of non-identical protein sequences.
 */
public class ProteinSequenceClusterer {
	private Structure structure = null;
	private Structure structure2 = null;
	private QuatSymmetryParameters parameters = null;
	
	private List<Atom[]> caUnaligned = new ArrayList<Atom[]>();
	private List<String> chainIds = new ArrayList<String>();
	private List<Integer> modelNumbers = new ArrayList<Integer>();
	private List<String> sequences = new ArrayList<String>();	
	private List<SequenceAlignmentCluster> seqClusters = new ArrayList<SequenceAlignmentCluster>();
	private int nucleicAcidChainCount = 0;
	private boolean modified = true;

	public ProteinSequenceClusterer(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
	}
	
	public ProteinSequenceClusterer(Structure structure1, Structure structure2,  QuatSymmetryParameters parameters) {
		this.structure = structure1;
		this.structure2 = structure2;
		this.parameters = parameters;
	}
	
	public List<SequenceAlignmentCluster> getSequenceAlignmentClusters() {
		run();
		return seqClusters;
	}
	
	public int getProteinChainCount() {
		run();
		return sequences.size();
	}
	
	/**
	 * @return the nucleicAcidChainCount
	 */
	public int getNucleicAcidChainCount() {
		run();
		return nucleicAcidChainCount;
	}

	public static void sortSequenceClustersBySize(List<SequenceAlignmentCluster> clusters) {
		Collections.sort(clusters, new Comparator<SequenceAlignmentCluster>() {
			@Override
			public int compare(SequenceAlignmentCluster c1, SequenceAlignmentCluster c2) {
				int sign = Math.round(Math.signum(c2.getSequenceCount() - c1.getSequenceCount()));
				if (sign != 0) {
					return sign;
				}
				return Math.round(Math.signum(c2.getSequenceAlignmentLength() - c1.getSequenceAlignmentLength()));
			}
		});
	}
	
	private void run() {
		if (modified) {
			extractProteinChains();
			clusterChains();
			modified = false;
		}
	}
	/**
	 * Populate all fields. If two structres are give, concatenate their chains.
	 */
	private void extractProteinChains() {
		ProteinChainExtractor extractor = new ProteinChainExtractor(structure,  parameters);
		caUnaligned = extractor.getCalphaTraces();
		chainIds  = extractor.getChainIds();
		sequences = extractor.getSequences();
		modelNumbers = extractor.getModelNumbers();
		nucleicAcidChainCount = extractor.getNucleicAcidChainCount();
		
		if (structure2 != null) {
			extractor = new ProteinChainExtractor(structure2,  parameters);
			caUnaligned.addAll(extractor.getCalphaTraces());
			chainIds.addAll(extractor.getChainIds());
			sequences.addAll(extractor.getSequences());
			modelNumbers.addAll(extractor.getModelNumbers());
		}
	}
	
	/**
	 * Cluster chains based on their sequence. Initializes seqClusters to the set
	 * of non-identical sequences.
	 */
	private void clusterChains() {
		boolean[] processed = new boolean[caUnaligned.size()];
		Arrays.fill(processed, false);
	
		for (int i = 0; i < caUnaligned.size(); i++) {
			if (processed[i]) {
				continue;
			}
			processed[i] = true;
			// create new sequence cluster
            UniqueSequenceList seqList = new UniqueSequenceList(caUnaligned.get(i), chainIds.get(i), modelNumbers.get(i), 0, sequences.get(i));
            SequenceAlignmentCluster seqCluster = new SequenceAlignmentCluster(parameters);
            seqCluster.addUniqueSequenceList(seqList);	
            seqClusters.add(seqCluster);
			
            for (int j = i + 1; j < caUnaligned.size(); j++) {
            	if (processed[j]) {
            		continue;
            	}
            	// Mark any future identical sequences as processed
            	for (SequenceAlignmentCluster c: seqClusters) {
            			if (c.identityMatch(caUnaligned.get(j), chainIds.get(j), modelNumbers.get(j), 0, sequences.get(j))) {
            				processed[j] = true;
            				//System.out.println("found identity match: " + i + " - " + j);
            				break;
            			}
            	}
            }
		}
		sortSequenceClustersBySize(seqClusters);
	}
}
