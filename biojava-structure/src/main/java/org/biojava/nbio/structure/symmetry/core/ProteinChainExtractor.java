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

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Extracts information about all the chains in a structure, including
 * chain Ids, sequences, and atoms. Includes both protein and nucleic acid chains.
 */
public class ProteinChainExtractor  {
	
	private static final Logger logger = LoggerFactory.getLogger(ProteinChainExtractor.class);
	
	private Structure structure = null;
	private QuatSymmetryParameters parameters = null;
	private boolean modified = true;
	private int adjustedMinimumSequenceLength = 0;

	private List<Atom[]> cAlphaTrace = new ArrayList<Atom[]>();	
	private List<String> chainIds = new ArrayList<String>();
	private List<Integer> modelNumbers = new ArrayList<Integer>();
	private List<String> sequences = new ArrayList<String>();
	private int nucleicAcidChainCount = 0;

	public ProteinChainExtractor(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
		modified = true;
	}

	public List<Atom[]> getCalphaTraces() {
		run();
		return cAlphaTrace;
	}

	public List<String> getChainIds() {
		run();
		return chainIds;
	}

	public List<Integer> getModelNumbers() {
		run();
		return modelNumbers;
	}

	public List<String> getSequences() {
		run();
		return sequences;
	}

	/**
	 * @return the nucleicAcidChainCount
	 */
	public int getNucleicAcidChainCount() {
		run();
		return nucleicAcidChainCount;
	}

	public int getAdjustedMinimumSequenceLength() {
		run();
		return adjustedMinimumSequenceLength;
	}

	private void run() {
		if (modified) {
			calcAdjustedMinimumSequenceLength();
			extractProteinChains();
			modified = false;
		}
	}

	private void extractProteinChains() {
		int models = 1;
		if (structure.isBiologicalAssembly()) {
			models = structure.nrModels();
		}

		if (parameters.isVerbose()) {
			System.out.println("Protein chains used in calculation:");
			System.out.println("Adjusted minimum sequence length: " + adjustedMinimumSequenceLength);
		}

		for (int i = 0; i < models; i++) {
			for (Chain c : structure.getChains(i)) {
				if (isNucleicAcidChain(c)) {
					nucleicAcidChainCount++;
				} //TODO Should we break here for DNA? If "CA" atoms are present, could cause bugs. -Spencer 9-2015 
				Atom[] ca = StructureTools.getAtomCAArray(c);
				ca = retainStandardAminoAcidResidues(ca);

				if (ca.length >= adjustedMinimumSequenceLength) {
					if (parameters.isVerbose()) {
						System.out.println("Chain " + c.getChainID() + " Calpha atoms: " + ca.length + " seqres: " + c.getSeqResSequence());
					}

					cAlphaTrace.add(ca);
					chainIds.add(c.getChainID());
					modelNumbers.add(i);
					sequences.add(replaceQuestionMarks(c.getSeqResSequence()));
				}
			}
		}
	}

	/**
	 * Returns an adapted minimum sequence length. This method ensure that
	 * structure that only have short chains are not excluded by the
	 * minimumSequenceLength cutoff value;
	 * @return
	 */
	private void calcAdjustedMinimumSequenceLength() {
		int models = 1;
		if (structure.isBiologicalAssembly()) {
			models = structure.nrModels();
		}

		int maxLength = Integer.MIN_VALUE;
		int minLength = Integer.MAX_VALUE;

		List<Integer> lengths = new ArrayList<Integer>();
		for (int i = 0; i < models; i++) {
			for (Chain c : structure.getChains(i)) {
				Atom[] ca = StructureTools.getAtomCAArray(c);
				ca = retainStandardAminoAcidResidues(ca);
				if (ca.length >= parameters.getAbsoluteMinimumSequenceLength()) {
					maxLength = Math.max(ca.length, maxLength);
					minLength = Math.min(ca.length, minLength);
					lengths.add(ca.length);
				}
			}
		}

		adjustedMinimumSequenceLength = parameters.getMinimumSequenceLength();

		if (lengths.size() < 2) {
			return;
		}

		double median = 0;
		Collections.sort(lengths);
		if (lengths.size() %2 == 1) {
			int middle = (lengths.size()-1) / 2;
			median = lengths.get(middle);
		} else {
			int middle2 = lengths.size()/2;
			int middle1 = middle2-1;
			median = 0.5 * (lengths.get(middle1) + lengths.get(middle2));
		}

		if (minLength >= median * parameters.getMinimumSequenceLengthFraction()) {
			adjustedMinimumSequenceLength = Math.min(minLength, parameters.getMinimumSequenceLength());
		}
	}
	
	private boolean isNucleicAcidChain(Chain chain) {
		int count = 0;
		for (Group group: chain.getAtomGroups()) {
			PolymerType type = group.getChemComp().getPolymerType();
			if (type != null && (type.equals(PolymerType.dna) || type.equals(PolymerType.rna) || type.equals(PolymerType.dnarna))) {
				count++;
			}
		}
		return count > 3;
	}

	// In some cases "?" are in the sequence string. Example 2WS1. This is caused
	// because the chemical component file YNM doesn't contain a one-letter code.
	private String replaceQuestionMarks(String sequence) {
		return sequence.replaceAll("\\?", "X");
	}

	private Atom[] retainStandardAminoAcidResidues(Atom[] atoms) {
		List<Atom> atomList = new ArrayList<Atom>(atoms.length);
		for (Atom atom: atoms) {
			Group group = atom.getGroup();
			if (group.getPDBName().equalsIgnoreCase("UNK")) {
				continue;
			}
			if (! isAminoAcid(group)) {
				continue;
			}
			atomList.add(atom);
		}
		return atomList.toArray(new Atom[atomList.size()]);
	}

	private boolean isAminoAcid(Group group) {
		ChemComp cc= group.getChemComp();
		if (cc.getResidueType() == null) {
			logger.warn("null residue type for: " + group.getPDBName());
			return false;
		}
		return (cc.getResidueType().equals(ResidueType.lPeptideLinking) || cc.getResidueType().equals(ResidueType.glycine));
	}
}
