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
 * Created on Aug 2, 2010
 * Author: Jianjiong Gao
 *
 */

package org.biojava.nbio.protmod.structure;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.mmcif.MetalBondParser;
import org.biojava.nbio.structure.io.mmcif.chem.MetalBondDistance;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public final class StructureUtil {
	private StructureUtil() {
		throw new AssertionError();
	}

	/**
	 *
	 * @param group a {@link Group} in structure.
	 * @param isAminoAcid true if it is an amino acid.
	 * @return the {@link StructureGroup} of the group.
	 */
	public static StructureGroup getStructureGroup(Group group, boolean isAminoAcid) {
		ResidueNumber resNum = group.getResidueNumber();
		return new StructureGroup(resNum, group.getPDBName(), isAminoAcid);
	}

	/**
	 *
	 * @param atom a {@link Atom} in structure.
	 * @param isParentAminoAcid true if the containing group is an amino acid.
	 * @return the {@link StructureAtom} of the atom.
	 */
	public static StructureAtom getStructureAtom(Atom atom, boolean isParentAminoAcid) {

		Group g = atom.getGroup();
		String chainId = g.getChainId();
		StructureGroup strucGroup = getStructureGroup(g, isParentAminoAcid);
		strucGroup.setChainId(chainId);
		return new StructureAtom(strucGroup, atom.getName());
	}

	/**
	 *
	 * @param atom1 the first {@link Atom} in structure.
	 * @param isParentAminoAcid1 true if the first containing group is an amino acid..
	 * @param atom2 the second {@link Atom} in structure.
	 * @param isParentAminoAcid2 true if the second containing group is an amino acid..
	 * @return the {@link StructureAtomLinkage} of the two atoms.
	 */
	public static StructureAtomLinkage getStructureAtomLinkage(Atom atom1,
			boolean isParentAminoAcid1, Atom atom2, boolean isParentAminoAcid2) {
		StructureAtom strucAtom1 = getStructureAtom(atom1, isParentAminoAcid1);
		StructureAtom strucAtom2 = getStructureAtom(atom2, isParentAminoAcid2);
		double distance = getAtomDistance(atom1, atom2);
		return new StructureAtomLinkage(strucAtom1, strucAtom2, distance);
	}

	/**
	 *
	 * @param atom1 the first {@link Atom} in structure.
	 * @param atom2 the second {@link Atom} in structure.
	 * @return the distance between the two atoms in Angstrom.
	 */
	public static double getAtomDistance(Atom atom1, Atom atom2) {
		double distance;

		distance = Calc.getDistance(atom1, atom2);

		return distance;
	}

	/**
	 * Find a linkage between two groups within tolerance of bond length,
	 * from potential atoms.
	 * @param group1 the first {@link Group}.
	 * @param group2 the second {@link Group}.
	 * @param potentialNamesOfAtomOnGroup1 potential names of the atom on the first group.
	 * 		  If null, search all atoms on the first group.
	 * @param potentialNamesOfAtomOnGroup2 potential names of the atom on the second group.
	 * 		  If null, search all atoms on the second group.
	 * @param ignoreNCLinkage true to ignore all N-C linkages
	 * @param bondLengthTolerance bond length error tolerance.
	 * @return an array of two Atoms that form bond between each other
	 *  if found; null, otherwise.
	 */
	public static Atom[] findNearestAtomLinkage(final Group group1, final Group group2,
			List<String> potentialNamesOfAtomOnGroup1, List<String> potentialNamesOfAtomOnGroup2,
			final boolean ignoreNCLinkage, double bondLengthTolerance) {


		List<Atom[]> linkages = findAtomLinkages(group1, group2,
				potentialNamesOfAtomOnGroup1, potentialNamesOfAtomOnGroup2,
				ignoreNCLinkage, bondLengthTolerance);

		Atom[] ret = null;
		double minDistance = Double.POSITIVE_INFINITY;

		for (Atom[] linkage : linkages) {
			double distance;

			distance = Calc.getDistance(linkage[0], linkage[1]);


			if (distance < minDistance) {
				minDistance = distance;
				ret = linkage;
			}
		}

		return ret;
	}

	/**
	 * Find linkages between two groups within tolerance of bond length,
	 * from potential atoms.
	 * @param group1 the first {@link Group}.
	 * @param group2 the second {@link Group}.
	 * @param ignoreNCLinkage true to ignore all N-C linkages
	 * @param bondLengthTolerance bond length error tolerance.
	 * @return a list, each element of which is an array of two Atoms that form bond
	 * between each other.
	 */
	public static List<Atom[]> findAtomLinkages(final Group group1,
			final Group group2, final boolean ignoreNCLinkage,
			final double bondLengthTolerance) {
		return findAtomLinkages(group1, group2,
				null, null, ignoreNCLinkage, bondLengthTolerance);
	}

	/**
	 * Find linkages between two groups within tolerance of bond length,
	 * from potential atoms.
	 * @param group1 the first {@link Group}.
	 * @param group2 the second {@link Group}.
	 * @param potentialNamesOfAtomOnGroup1 potential names of the atom on the first group.
	 * 		  If null, search all atoms on the first group.
	 * @param potentialNamesOfAtomOnGroup2 potential names of the atom on the second group.
	 * 		  If null, search all atoms on the second group.
	 * @param ignoreNCLinkage true to ignore all N-C linkages
	 * @param bondLengthTolerance bond length error tolerance.
	 * @return a list, each element of which is an array of two Atoms that form bond
	 * between each other.
	 */
	public static List<Atom[]> findAtomLinkages(final Group group1,
			final Group group2,
			List<String> potentialNamesOfAtomOnGroup1,
			List<String> potentialNamesOfAtomOnGroup2,
			final boolean ignoreNCLinkage,
			final double bondLengthTolerance) {
		if (group1==null || group2==null) {
			throw new IllegalArgumentException("Null group(s).");
		}

		if (bondLengthTolerance<0) {
			throw new IllegalArgumentException("bondLengthTolerance cannot be negative.");
		}

		List<Atom[]> ret = new ArrayList<Atom[]>();

		if (potentialNamesOfAtomOnGroup1 == null) {
			// if empty name, search for all atoms
			potentialNamesOfAtomOnGroup1 = getAtomNames(group1);
		}

		if (potentialNamesOfAtomOnGroup2 == null) {
			// if empty name, search for all atoms
			potentialNamesOfAtomOnGroup2 = getAtomNames(group2);
		}

		for (String namesOfAtomOnGroup1 : potentialNamesOfAtomOnGroup1) {
			for (String namesOfAtomOnGroup2 : potentialNamesOfAtomOnGroup2) {
				Atom[] atoms = findLinkage(group1, group2, namesOfAtomOnGroup1,
						namesOfAtomOnGroup2, bondLengthTolerance);
				if (atoms != null) {
					if (ignoreNCLinkage &&
							((atoms[0].getName().equals("N") && atoms[1].getName().equals("C"))
									|| (atoms[0].getName().equals("C") && atoms[1].getName().equals("N")))
							) {
						continue;
					}

					ret.add(atoms);
				}
			}
		}

		return ret;
	}

	/**
	 * Find a linkage between two groups within tolerance of bond length.
	 * @param group1 the first {@link Group}.
	 * @param group2 the second {@link Group}.
	 * @param nameOfAtomOnGroup1 atom name of the first group.
	 * @param nameOfAtomOnGroup2 atom name of the second group.
	 * @param bondLengthTolerance bond length error tolerance.
	 * @return an array of two Atoms that form bond between each other
	 *  if found; null, otherwise.
	 */
	public static Atom[] findLinkage(final Group group1, final Group group2,
			String nameOfAtomOnGroup1, String nameOfAtomOnGroup2,
			double bondLengthTolerance) {

		Atom[] ret = new Atom[2];

		ret[0] = group1.getAtom(nameOfAtomOnGroup1);
		ret[1] = group2.getAtom(nameOfAtomOnGroup2);

		if (ret[0]==null || ret[1]==null) {
			return null;
		}


		Atom a1 = ret[0];
		Atom a2 = ret[1];

		boolean hasBond =  a1.hasBond(a2);

		if ( hasBond ) {

			return ret;
		}

		// is it a metal ?

		if ( a1.getElement().isMetal() || a2.getElement().isMetal()){

			MetalBondDistance defined = getMetalDistanceCutoff(a1.getElement().name(),a2.getElement().name());

			if ( defined != null) {

				if (hasMetalBond(a1, a2, defined))
					return ret;
				else
					return null;
			}

		}

			// not a metal

			double distance = Calc.getDistance(a1, a2);

			float radiusOfAtom1 = ret[0].getElement().getCovalentRadius();
			float radiusOfAtom2 = ret[1].getElement().getCovalentRadius();

			if (Math.abs(distance - radiusOfAtom1 - radiusOfAtom2)
					> bondLengthTolerance) {
				return null;
			}


		return ret;
	}

	private static boolean hasMetalBond(Atom a1, Atom a2, MetalBondDistance definition) {

		double distance = Calc.getDistance(a1,a2);

		Float min = definition.getLowerLimit();
		Float max = definition.getUpperLimit();

		return ( min < distance && max > distance);

	}

	private static MetalBondDistance getMetalDistanceCutoff(String name1, String name2) {
		Map<String,List<MetalBondDistance>> defs= MetalBondParser.getMetalBondDefinitions();

		List<MetalBondDistance> distances = defs.get(name1);

		if ( distances == null){

			distances = defs.get(name2);
			String tmp = name1;
			 name2 = name1;
			 name1 = tmp;
		}
		if ( distances == null){
			return null;
		}

		for  ( MetalBondDistance d : distances){
			if ( name1.equalsIgnoreCase(d.getAtomType1()) && name2.equalsIgnoreCase(d.getAtomType2()) )
				return d;
		}

		// no matching atom definitions found
		return null;

	}

	/**
	 *
	 * @param group a {@link Group}.
	 * @return all atom names of the group.
	 */
	public static List<String> getAtomNames(Group group) {
		List<Atom> atoms = group.getAtoms();
		if (atoms == null) {
			return Collections.emptyList();
		}

		int n = atoms.size();
		List<String> ret = new ArrayList<String>(n);
		for (int i=0; i<n; i++) {
			ret.add(atoms.get(i).getName());
		}

		return ret;
	}

	/**
	 * Get all amino acids in a chain.
	 * @param chain
	 * @return
	 */
	public static List<Group> getAminoAcids(Chain chain) {

		List<Group> gs = new ArrayList<>();

		for ( Group g : chain.getAtomGroups()){

			if ( g.isAminoAcid())
				gs.add(g);

		}

		return gs;

	}

}
