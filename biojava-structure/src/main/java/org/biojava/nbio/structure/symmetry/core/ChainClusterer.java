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

import javax.vecmath.Point3d;
import java.util.*;

/**
 * Wraps a sequence clustering with structural information
 */
public class ChainClusterer  {	
	private List<SequenceAlignmentCluster> seqClusters = new ArrayList<SequenceAlignmentCluster>();	
	private boolean modified = true;

	private List<Atom[]> caAligned = new ArrayList<Atom[]>();
	private List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
	
	public ChainClusterer(List<SequenceAlignmentCluster> seqClusters) {
		this.seqClusters = seqClusters;
		this.modified = true;
	}
	
	public List<Point3d[]> getCalphaCoordinates() {
        run();
		return caCoords;
	}
	
	public List<Atom[]> getCalphaTraces() {
		run();
		return caAligned;
	}
	
	public List<String> getChainIds() {
		run();
		List<String> chainIdList = new ArrayList<String>();

		for (int i = 0; i < seqClusters.size(); i++) {
	        SequenceAlignmentCluster cluster = seqClusters.get(i);
	        for (String chainId: cluster.getChainIds()) {
	        	chainIdList.add(chainId);
	        }
		}
		return chainIdList;
	}
	
	
	public List<Integer> getModelNumbers() {
		run();
		List<Integer> modNumbers = new ArrayList<Integer>();

		for (int i = 0; i < seqClusters.size(); i++) {
	        SequenceAlignmentCluster cluster = seqClusters.get(i);
	        for (Integer number: cluster.getModelNumbers()) {
	        	modNumbers.add(number);
	        }
		}
		return modNumbers;
	}
	
	public String getStoichiometry() {
		run();
		StringBuilder formula = new StringBuilder();
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

		for (int i = 0; i < seqClusters.size(); i++) {
			String c = "?";
			if (i < alpha.length()) {
				c = alpha.substring(i, i+1);
			}
			formula.append(c);
			int multiplier = seqClusters.get(i).getSequenceCount();
			if (multiplier > 1) {
				formula.append(multiplier);
			}
		}
		return formula.toString();
	}

	/**
	 * Get valid symmetry order for this stoichiometry.
	 * @return
	 */
	public List<Integer> getFolds() {
		run();
		List<Integer> stoichiometry = new ArrayList<Integer>(seqClusters.size());
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			stoichiometry.add(seqCount);
		}
		return getValidFolds(stoichiometry);
	}
	/**
	 * Find valid symmetry orders for a given stoichiometry. For instance,
	 * an A6B4 protein would give [1,2] because (A6B4)1 and (A3B2)2 are valid
	 * decompositions.
	 * @param stoichiometry List giving the number of copies in each chain cluster
	 * @return The common factors of the stoichiometry
	 */
	public static List<Integer> getValidFolds(List<Integer> stoichiometry){
		List<Integer> denominators = new ArrayList<Integer>();

		int nChains = Collections.max(stoichiometry);
		
		// Remove duplicate stoichiometries
		Set<Integer> nominators = new TreeSet<Integer>(stoichiometry);

		// find common denominators
		for (int d = 1; d <= nChains; d++) {
			boolean isDivisable=true;
			for (Integer n : nominators) {
				if (n % d != 0) {
					isDivisable = false;
					break;
				}
			}
			if(isDivisable) {
				denominators.add(d);
			}
		}
		return denominators;
	}
	
	public List<Integer> getSequenceClusterIds() {
		run();
		List<Integer> list = new ArrayList<Integer>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			for (int i = 0; i < seqCount; i++) {
				list.add(id);
			}
		}
		return list;
	}
	
	
	public int getSequenceClusterCount() {
		run();
		return seqClusters.size();
	}
	
	public List<SequenceAlignmentCluster> getSequenceAlignmentClusters() {
		return seqClusters;
	}
	
	public List<Boolean> getPseudoStoichiometry() {
		run();
		List<Boolean> list = new ArrayList<Boolean>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			Boolean pseudo = seqClusters.get(id).isPseudoStoichiometric();
			for (int i = 0; i < seqCount; i++) {
				list.add(pseudo);
			}
		}
		return list;
	}
	
	public List<Double> getMinSequenceIdentity() {
		run();
		List<Double> list = new ArrayList<Double>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			double minSequenceIdentity = seqClusters.get(id).getMinSequenceIdentity();
			for (int i = 0; i < seqCount; i++) {
				list.add(minSequenceIdentity);
			}
		}
		return list;
	}
	
	public List<Double> getMaxSequenceIdentity() {
		run();
		List<Double> list = new ArrayList<Double>();
		
		for (int id = 0; id < seqClusters.size(); id++) {
			int seqCount = seqClusters.get(id).getSequenceCount();
			double maxSequenceIdentity = seqClusters.get(id).getMaxSequenceIdentity();
			for (int i = 0; i < seqCount; i++) {
				list.add(maxSequenceIdentity);
			}
		}
		return list;
	}
	
	@Override
	public String toString() {
		run();
		StringBuilder builder = new StringBuilder();
		builder.append("Sequence alignment clusters: " + seqClusters.size());
		builder.append("\n");
		for (SequenceAlignmentCluster s: seqClusters) {
			builder.append("# seq: ");
			builder.append(s.getSequenceCount());
			builder.append(" alignment length: ");
			builder.append(s.getSequenceAlignmentLength());
			builder.append("\n");
		}
		return builder.toString();
	}
	
	private void run() {
		if (modified) {
			modified = false;
			calcAlignedSequences();
			createCalphaTraces();
		}
	}
	

	private void calcAlignedSequences() {
		caAligned = new ArrayList<Atom[]>();
		for (SequenceAlignmentCluster cluster: seqClusters) {
			caAligned.addAll(cluster.getAlignedCalphaAtoms());	
		}
	}
	
	private void createCalphaTraces() {
		for (Atom[] atoms: caAligned) {
			Point3d[] trace = new Point3d[atoms.length];
			for (int j = 0; j < atoms.length; j++) {
				trace[j] = new Point3d(atoms[j].getCoords());
			}
			caCoords.add(trace);
		}
	}
}
