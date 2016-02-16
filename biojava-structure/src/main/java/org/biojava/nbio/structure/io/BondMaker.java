/*
 *                  BioJava development code
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
 * Created on Mar. 6, 2014
 *
 */
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompBond;

import java.util.ArrayList;
import java.util.List;

/**
 * Adds polymer bonds for peptides and nucleotides based on distance cutoffs and
 * intra-group (residue) bonds based on data from the Chemical Component Dictionary
 * to the Structure object.
 * 
 * TODO the current implementation adds bonds to the first model only. This
 * should be sufficient for homogeneous models, but here are a few inhomogeneous models
 * in the PDB. A better handling of models should be considered in the future.
 * 
 * @author Peter Rose
 * @author Ulysse Carion
 *
 */
public class BondMaker {
	/**
	 * Maximum peptide (C - N) bond length considered for bond formation
	 */
	private static final double MAX_PEPTIDE_BOND_LENGTH = 1.8;
	/**
	 * Maximum nucleotide (P - O3') bond length considered for bond formation
	 */
	private static final double MAX_NUCLEOTIDE_BOND_LENGTH = 2.1;

	private Structure structure = null;

	public BondMaker(Structure structure) {
		this.structure = structure;
	}

	public void makeBonds() {
		formPeptideBonds();
		formNucleotideBonds();
		formIntraResidueBonds();
		trimBondLists();
	}

	private void formPeptideBonds() {
		for (Chain chain : structure.getChains()) {
			List<Group> groups = chain.getSeqResGroups();

			for (int i = 0; i < groups.size() - 1; i++) {
				if (!(groups.get(i) instanceof AminoAcidImpl)
						|| !(groups.get(i + 1) instanceof AminoAcidImpl))
					continue;

				AminoAcidImpl tail = (AminoAcidImpl) groups.get(i);
				AminoAcidImpl head = (AminoAcidImpl) groups.get(i + 1);

				// atoms with no residue number don't have atom information
				if (tail.getResidueNumber() == null
						|| head.getResidueNumber() == null) {
					continue;
				}

				Atom carboxylC;
				Atom aminoN;

				carboxylC = tail.getC();
				aminoN = head.getN();


				if (carboxylC == null || aminoN == null) {
					// some structures may be incomplete and not store info
					// about all of their atoms

					continue;
				}


				if (Calc.getDistance(carboxylC, aminoN) < MAX_PEPTIDE_BOND_LENGTH) {
					new BondImpl(carboxylC, aminoN, 1);
				}

			}
		}
	}

	private void formNucleotideBonds() {
		for (Chain chain : structure.getChains()) {
			List<Group> groups = chain.getSeqResGroups();

			for (int i = 0; i < groups.size() - 1; i++) {
				if (!(groups.get(i) instanceof NucleotideImpl)
						|| !(groups.get(i + 1) instanceof NucleotideImpl))
					continue;

				NucleotideImpl tail = (NucleotideImpl) groups.get(i);
				NucleotideImpl head = (NucleotideImpl) groups.get(i + 1);

				// atoms with no residue number don't have atom information
				if (tail.getResidueNumber() == null
						|| head.getResidueNumber() == null) {
					continue;
				}

				Atom phosphorous = tail.getP();
				Atom oThreePrime = head.getO3Prime();

				if (phosphorous == null || oThreePrime == null) {
					continue;
				}


				if (Calc.getDistance(phosphorous, oThreePrime) < MAX_NUCLEOTIDE_BOND_LENGTH) {
					new BondImpl(phosphorous, oThreePrime, 1);
				}

			}
		}
	}

	private void formIntraResidueBonds() {
		for (Chain chain : structure.getChains()) {
			List<Group> groups = chain.getAtomGroups();

			for (Group group : groups) {
				// atoms with no residue number don't have atom information
				if (group.getResidueNumber() == null) {
					continue;
				}

				ChemComp aminoChemComp = ChemCompGroupFactory.getChemComp(group
						.getPDBName());

				for (ChemCompBond chemCompBond : aminoChemComp.getBonds()) {

					Atom a = group.getAtom(chemCompBond.getAtom_id_1());
					Atom b = group.getAtom(chemCompBond.getAtom_id_2());
					if ( a != null && b != null){

						int bondOrder = chemCompBond.getNumericalBondOrder();

						new BondImpl(a, b, bondOrder);
					} else  {

						// Some of the atoms were missing. That's fine, there's
						// nothing to do in this case.
					}
				}
			}
		}
	}

	private void trimBondLists() {
		for (Chain chain : structure.getChains()) {
			for (Group group : chain.getAtomGroups()) {
				for (Atom atom : group.getAtoms()) {
					if (atom.getBonds().size() > 0) {
						((ArrayList<Bond>) atom.getBonds()).trimToSize();
					}
				}
			}
		}
	}

}
