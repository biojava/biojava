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
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.ce.CeParameters;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.seq.SmithWaterman3Daligner;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

/**
 * Represents a cluster of equivalent sequences
 *
 */
public class SequenceAlignmentCluster implements Cloneable {
	
	private static final Logger logger = LoggerFactory.getLogger(SequenceAlignmentCluster.class);
	
	private QuatSymmetryParameters parameters = null;
	private List<UniqueSequenceList> uniqueSequenceList = new ArrayList<UniqueSequenceList>();
	private List<Atom[]> alignedCAlphaAtoms = null;

	private int alignmentLength = 0;
	private double minSequenceIdentity = 1.0;
	private double maxSequenceIdentity = 0.0;
	
	private boolean modified = true;

	public SequenceAlignmentCluster (QuatSymmetryParameters parameters) {
		this.parameters = parameters;
	}
	
	public boolean isPseudoStoichiometric() {
		return minSequenceIdentity < parameters.getSequencePseudoSymmetryThreshold();
	}
	
	public double getMinSequenceIdentity() {
		if (! isPseudoStoichiometric()) {
			return 1.0;
		}
		return minSequenceIdentity;
	}
	
	public void setMinSequenceIdentity(double minSequenceIdentity) {
		this.minSequenceIdentity = minSequenceIdentity;
	}

	public double getMaxSequenceIdentity() {
		if (! isPseudoStoichiometric()) {
			return 1.0;
		}
		return maxSequenceIdentity;
	}

	public void setMaxSequenceIdentity(double maxSequenceIdentity) {
		this.maxSequenceIdentity = maxSequenceIdentity;
	}
	
	public void addUniqueSequenceList(UniqueSequenceList sequenceList) {
		uniqueSequenceList.add(sequenceList);
		modified = true;
	}

	/**
	 * @return the number of sequences which have been added to this cluster
	 */
	public int getSequenceCount() {
		return uniqueSequenceList.size();
	}
	
	public int getSequenceAlignmentLength() {
		run();
		return alignmentLength;
	}
	
	public List<UniqueSequenceList> getUniqueSequenceList() {
		return uniqueSequenceList;
	}

	public List<String> getChainIds() {
		List<String> ids = new ArrayList<String>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			ids.add(list.getChainId());
		}
		return ids;
	}
	
	public List<Integer> getModelNumbers() {
		List<Integer> numbers = new ArrayList<Integer>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			numbers.add(list.getModelNumber());
		}
		return numbers;
	}
	
	public List<Integer> getStructureIds() {
		List<Integer> numbers = new ArrayList<Integer>();
		for (UniqueSequenceList list: uniqueSequenceList) {
			numbers.add(list.getStructureId());
		}
		return numbers;
	}
	
	public List<Atom[]> getAlignedCalphaAtoms() {
		run();
		return alignedCAlphaAtoms;
	}
	
	/**
	 * Match a sequence to this cluster at 100% identity.
	 * 
	 * If the given sequence matches the cluster seed (100%), then performs an
	 * alignment to the seed and adds it to the {@link #getUniqueSequenceList()
	 * unique sequence list}.
	 *  
	 * @param cAlphaAtoms
	 * @param chainId
	 * @param modelNumber
	 * @param structureId
	 * @param sequence
	 * @return
	 */
	public boolean identityMatch(Atom[] cAlphaAtoms, String chainId, int modelNumber, int structureId, String sequence) {
		UniqueSequenceList u = uniqueSequenceList.get(0);

		// check for 100% identity match of reference sequence
		String refSequence = u.getSeqResSequence();
		boolean seqMatch = refSequence.equals(sequence);

		// if reference (SEQRES) sequences match 100%, 
		// find alignment of atom sequences by Smith-Waterman alignment
		if (seqMatch) {
			List<Integer> alig1 = new ArrayList<Integer>();
			List<Integer> alig2 = new ArrayList<Integer>();
			Atom[] referenceAtoms = u.getCalphaAtoms();
			int inCommon = 0;
			try {
				inCommon = alignIdenticalSequence(referenceAtoms, cAlphaAtoms, alig1, alig2);
			} catch (StructureException e) {
				// this happens in some cases like all X sequences, e.g. 1s1o or 1s4a
				logger.warn("Could not align identical sequences {}: {} and {}: {}. Chains won't be clustered together.",
						u.getChainId(),refSequence,chainId,sequence);
			}

			if (inCommon > 0) {
				UniqueSequenceList seqList = new UniqueSequenceList(cAlphaAtoms, chainId, modelNumber, structureId, sequence);
				seqList.setAlignment1(alig1);
				seqList.setAlignment2(alig2);
				//			System.out.println(alig1);
				//			System.out.println(alig2);
				addUniqueSequenceList(seqList);
				return true;
			}
		}

		return false;
	}
	
	public PairwiseAlignment getPairwiseAlignment(SequenceAlignmentCluster other) {
		PairwiseAlignment alignment = new PairwiseAlignment(this, other);
		
		Atom[] referenceAtoms1 = this.getUniqueSequenceList().get(0).getCalphaAtoms();
		Atom[] referenceAtoms2 = other.getUniqueSequenceList().get(0).getCalphaAtoms();
		
		double alignmentLengthFraction = (double)Math.min(referenceAtoms1.length, referenceAtoms2.length) /
				Math.max(referenceAtoms1.length, referenceAtoms2.length);
	
		if (alignmentLengthFraction < parameters.getAlignmentFractionThreshold()) {
			return null;
		}
		
		AFPChain afp = alignPairByStructure(referenceAtoms1, referenceAtoms2,parameters.isVerbose());
		if (afp == null) {
			return null;
		}
		
		if (! afp.isSignificantResult()) {
			return null;

    		// alternative: tmSCore:
    		// double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
    		// if ( tmScore < 0.35) {
    		// return null ...
		}
		
		int[][][] align = afp.getOptAln();
		if (align == null) {
			return null;
		}
			
    	alignmentLengthFraction = (double)afp.getOptLength()/Math.max(referenceAtoms1.length, referenceAtoms2.length);
    	alignment.setAlignmentLengthFraction(alignmentLengthFraction);
    	alignment.setRmsd(afp.getChainRmsd());
    	alignment.setSequenceIdentity(afp.getIdentity());
    	alignment.setAlignment(afp.getOptAln());
    	
		return alignment;
	}
	
	@Override
	public Object clone() {
	    SequenceAlignmentCluster copy = null;
		try {
			copy = (SequenceAlignmentCluster) super.clone();
		} catch (CloneNotSupportedException e) {
			logger.error("CloneNotSupportedException caught",e);
		}
		// deep copy sequences
		copy.uniqueSequenceList = new ArrayList<UniqueSequenceList>();
		for (UniqueSequenceList seq: this.getUniqueSequenceList()) {
			copy.addUniqueSequenceList((UniqueSequenceList) seq.clone());
		}
		return copy;
	}
	
	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for (UniqueSequenceList u: uniqueSequenceList) {
			builder.append(u.toString());
			builder.append("\n");
		}
		return builder.toString();
	}
	
	private void run() {
		if (modified) {
			alignedCAlphaAtoms = null;
			createAlignedCAlphaAtoms();
			modified = false;
		}
	}

	private static AFPChain alignPairBySequence(Atom[] ca1Seq, Atom[] ca2Seq) throws StructureException { 
		SmithWaterman3Daligner aligner = new SmithWaterman3Daligner();
		return aligner.align(ca1Seq, ca2Seq);
	}
	
	private static AFPChain alignPairByStructure(Atom[] ca1Seq, Atom[] ca2Seq, boolean verbose) {
       CeParameters params = new CeParameters();

        AFPChain afp = null;
		try {
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			afp = algorithm.align(ca1Seq,ca2Seq,params);
			if (verbose) {
				System.out.println(afp.toFatcat(ca1Seq, ca2Seq));
			}
		} catch (StructureException e) {
			logger.error("StructureException caught",e);
		}            
		return afp;
	}
	
	
	private static int alignIdenticalSequence(Atom[] ca1Seq, Atom[] ca2Seq, List<Integer> align1, List<Integer> align2) throws StructureException {
		AFPChain afp = alignPairBySequence(ca1Seq, ca2Seq);
		int[][][] align = afp.getOptAln();
		if (align == null) {
			return 0;
		}
		int len =  afp.getOptLength();

		List<Integer> delta = new ArrayList<Integer>();
		Set<Integer> unique = new HashSet<Integer>();

		for (int i = 0; i < len; i++) {
			Atom a1 = ca1Seq[align[0][0][i]];
			String residueName1 = a1.getGroup().getPDBName();
			Atom a2 = ca2Seq[align[0][1][i]];
			String residueName2 = a2.getGroup().getPDBName();
			if (residueName1.equals(residueName2)) {
			    int n1 = a1.getGroup().getResidueNumber().getSeqNum();
			    int n2 = a2.getGroup().getResidueNumber().getSeqNum();
			    delta.add(n2-n1);
			    unique.add(n2-n1);
			}
		}
		
		int offset = 0;
		int frequency = 0;
        for (Integer i: unique) {
        	int freq = Collections.frequency(delta, i);
        	if (freq > frequency) {
        		offset = i;
        		frequency = freq;
        	}
        }
        
        for (int i = 0; i < len; i++) {
        	Atom a1 = ca1Seq[align[0][0][i]];
			int n1 = a1.getGroup().getResidueNumber().getSeqNum();
			Atom a2 = ca2Seq[align[0][1][i]];
			int n2 = a2.getGroup().getResidueNumber().getSeqNum();
			if (n2 - offset == n1) {
				align1.add(align[0][0][i]);
				align2.add(align[0][1][i]);
			}
        }
//        System.out.println("PDB alignment: ");
//        System.out.println(align1);
//        System.out.println(align2);
        return align1.size();
	}
	
	private void createAlignedCAlphaAtoms() {
		List<Integer> indices = getReferenceResidueIndices();
		alignmentLength = indices.size();
		alignedCAlphaAtoms = new ArrayList<Atom[]>();
		for (UniqueSequenceList u: uniqueSequenceList) {
			List<Integer> alignment1 = u.getAlignment1();
			List<Integer> alignment2 = u.getAlignment2();
			List<Integer> alignmentIndices = new ArrayList<Integer>();
			for (int i = 0; i < alignment1.size(); i++) {
				int a1 = alignment1.get(i);
				if (indices.contains(a1)) {
					alignmentIndices.add(alignment2.get(i));
				}
			}
			Atom[] unalignedAtoms = u.getCalphaAtoms();
			Atom[] alignedAtoms = new Atom[alignmentIndices.size()];
			for (int j = 0; j < alignedAtoms.length; j++) {
				alignedAtoms[j] = unalignedAtoms[alignmentIndices.get(j)];
			}
			alignedCAlphaAtoms.add(alignedAtoms);
		}
	}
	
	private List<Integer> getReferenceResidueIndices() {
		List<Integer> indices = new ArrayList<Integer>(uniqueSequenceList.get(0).getAlignment1());
		for (UniqueSequenceList u: uniqueSequenceList) {
           indices.retainAll(u.getAlignment1());
		}
		return indices;
	}
}
