package org.biojava.bio.structure.symmetry.core;

import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.Structure;

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
