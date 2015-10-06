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

import org.biojava.nbio.structure.Structure;

import java.util.Collections;
import java.util.List;

/**
 * Clusters the chains of one or two structures by sequence.
 */
public class ClusterProteinChains {
	private Structure structure = null;
	private Structure structure2 = null;
	private QuatSymmetryParameters parameters = null;
	private ClusterMerger merger = null;
	private int proteinChainCount = 0;
	private int nucleicAcidChainCount = 0;

	public ClusterProteinChains(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
		run();
	}
	
	public ClusterProteinChains(Structure structure1, Structure structure2, QuatSymmetryParameters parameters) {
		this.structure = structure1;
		this.structure2 = structure2;
		this.parameters = parameters;
		run();
	}
	
	/**
	 * Get a non-redundent set of clusters for a given sequence cutoff
	 * @param sequenceIdentityThreshold
	 * @return
	 */
	public List<SequenceAlignmentCluster> getSequenceAlignmentClusters(double sequenceIdentityThreshold) {
		if (merger == null) {
			return Collections.emptyList();
		}
		return merger.getMergedClusters(sequenceIdentityThreshold);
	}
	
	/**
	 * @return the proteinChainCount
	 */
	public int getProteinChainCount() {
		return proteinChainCount;
	}

	/**
	 * @return the nucleicAcidChainCount
	 */
	public int getNucleicAcidChainCount() {
		return nucleicAcidChainCount;
	}

	private void run () {
		// cluster protein entities
		List<SequenceAlignmentCluster> seqClusters = null;
		
		if (structure2 == null) {
			ProteinSequenceClusterer clusterer = new ProteinSequenceClusterer(structure, parameters);
			seqClusters = clusterer.getSequenceAlignmentClusters();
			proteinChainCount = clusterer.getProteinChainCount();
			nucleicAcidChainCount = clusterer.getNucleicAcidChainCount();
		} else if (structure !=null && structure2 != null) {
			ProteinSequenceClusterer clusterer = new ProteinSequenceClusterer(structure, structure2, parameters);
			seqClusters = clusterer.getSequenceAlignmentClusters();
		}
		if (seqClusters == null  || seqClusters.size() == 0) {
			return;
		}
	
		// calculate pairwise aligment between protein clusters
		merger = new ClusterMerger(seqClusters, parameters);
		merger.calcPairwiseAlignments();
	}
}
