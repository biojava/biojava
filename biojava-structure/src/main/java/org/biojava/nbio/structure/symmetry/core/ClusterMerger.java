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

import java.util.*;

/**
 * Merges clusters based on their sequence identity. This class does the actual
 * agglomerative clustering calculation, while {@link SequenceAlignmentCluster}
 * stores the results.
 */
public class ClusterMerger {
	private List<SequenceAlignmentCluster> clusters = null;
    private QuatSymmetryParameters parameters = null;
    
    List<PairwiseAlignment> pairwiseAlignments = Collections.emptyList();

	public ClusterMerger(List<SequenceAlignmentCluster> clusters, QuatSymmetryParameters parameters) {
		this.clusters = clusters;
		this.parameters = parameters;
	}
	
	/**
	 * Aligns all pairs of input clusters, calculating their pairwise alignments
	 */
	public void calcPairwiseAlignments() {
		pairwiseAlignments = new ArrayList<PairwiseAlignment>();

		boolean[] merged = new boolean[clusters.size()];
		Arrays.fill(merged, false);

		for (int i = 0, n = clusters.size(); i < n-1; i++) {
			if (! merged[i]) {
				SequenceAlignmentCluster c1 = clusters.get(i);
				for (int j = i + 1; j < n; j++) {
					SequenceAlignmentCluster c2 = clusters.get(j);
					PairwiseAlignment alignment = c1.getPairwiseAlignment(c2);
					if (alignment != null && 
							alignment.getAlignmentLengthFraction() >= parameters.getAlignmentFractionThreshold() &&
							alignment.getRmsd() <= parameters.getRmsdThreshold()) {
						merged[j] = true;
						pairwiseAlignments.add(alignment);
						if (parameters.isVerbose()) {
							System.out.println("ClusterMerger: pairwise cluster alignment: " + i + "-" + j + " seq. identity: " + alignment.getSequenceIdentity());
							System.out.println(alignment);
						}
					}
				}
			}
		}
	}
	
	/**
	 * Combine clusters based on the given sequence identity
	 * @param sequenceIdentityCutoff
	 * @return
	 */
	public List<SequenceAlignmentCluster> getMergedClusters(double sequenceIdentityCutoff) {		
		List<SequenceAlignmentCluster> mergedClusters = new ArrayList<SequenceAlignmentCluster>();
		Map<SequenceAlignmentCluster, Integer> map = getClusterMap();
		
		boolean[] processed = new boolean[clusters.size()];
		Arrays.fill(processed, false);
		
		for (int i = 0, n = clusters.size(); i < n; i++) {
			SequenceAlignmentCluster cluster = clusters.get(i);
			SequenceAlignmentCluster clone = null;
			if (! processed[i]) {
				clone = (SequenceAlignmentCluster) cluster.clone();
				mergedClusters.add(clone);
				processed[i] = true;
			}
			
			for (PairwiseAlignment alignment: pairwiseAlignments) {
				if (alignment.getCluster1() == cluster && alignment.getSequenceIdentity() >= sequenceIdentityCutoff) {
					clone.setMinSequenceIdentity(Math.min(clone.getMinSequenceIdentity(), alignment.getSequenceIdentity()));
					clone.setMaxSequenceIdentity(Math.max(clone.getMaxSequenceIdentity(), alignment.getSequenceIdentity()));
					combineClusters(clone, alignment);
					int index = map.get(alignment.getCluster2());
					processed[index] = true;
				}
			}
		}
		
		ProteinSequenceClusterer.sortSequenceClustersBySize(mergedClusters);
		return mergedClusters;
	}
	
	
	private Map<SequenceAlignmentCluster, Integer> getClusterMap() {
		 Map<SequenceAlignmentCluster, Integer> map = new HashMap<SequenceAlignmentCluster, Integer>();
		 for (int i = 0, n = clusters.size(); i < n; i++) {
			 map.put(clusters.get(i), i);
		 }
		 return map;
	}
	
	private void combineClusters(SequenceAlignmentCluster c1, PairwiseAlignment alignment) {
		SequenceAlignmentCluster c2 = (SequenceAlignmentCluster) alignment.getCluster2().clone();
		int[][][] align = alignment.getAlignment();

		// add alignment for reference sequence
		UniqueSequenceList u =c2.getUniqueSequenceList().get(0);
		// set new sequence alignment
		List<Integer> align1 = new ArrayList<Integer>(align[0][0].length);
		for (Integer a1: align[0][0]) {
			align1.add(a1);
		}
		u.setAlignment1(align1);

		List<Integer> align2 = new ArrayList<Integer>(align[0][1].length);
		for (Integer a2: align[0][1]) {
			align2.add(a2);
		}
		u.setAlignment2(align2);	
		c1.addUniqueSequenceList(u);
		
		// note, i starts at 1 (ONE), since i = 0 corresponds to reference sequence, 
		// which has already been processed above
		for (int i = 1; i < c2.getUniqueSequenceList().size(); i++) {
			u =c2.getUniqueSequenceList().get(i);
			List<Integer> oldAlign1 = u.getAlignment1();
			List<Integer> oldAlign2 = u.getAlignment2();
			List<Integer> newAlign1 = new ArrayList<Integer>();
			List<Integer> newAlign2 = new ArrayList<Integer>();
			for (int j = 0; j < align2.size(); j++) {
				Integer element = align2.get(j);
				Integer index = oldAlign1.indexOf(element);
				// map alignment to first reference alignment
				if (index != null && index >= 0) {
					newAlign1.add(align1.get(j));
					newAlign2.add(oldAlign2.get(index));
				}
			}
			u.setAlignment1(newAlign1);
			u.setAlignment2(newAlign2);	
			c1.addUniqueSequenceList(u);
		}
	}
}
