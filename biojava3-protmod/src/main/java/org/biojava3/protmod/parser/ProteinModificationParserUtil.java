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
 * Created on Jun 24, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;

final class ProteinModificationParserUtil {
	/**
	 * Utility class. All methods are static.
	 */
	private ProteinModificationParserUtil() {
		throw new AssertionError();
	}
	
	/**
	 * Find a linkage between two groups within tolerance of bond length,
	 * from potential atoms.
	 * @param group1 the first {@link Group}.
	 * @param isGroup1AminoAcid true if group1 is an amino acid.
	 * @param group2 the second {@link Group}.
	 * @param isGroup2AminoAcid true if group2 is an amino acid.
	 * @param potentialNamesOfAtomOnGroup1 potential names of the atom on the first group.
	 * 		  If null, search all atoms on the first group.
	 * @param potentialNamesOfAtomOnGroup2 potential names of the atom on the second group.
	 * 		  If null, search all atoms on the second group.
	 * @param bondLengthTolerance bond length error tolerance.
	 * @return an array of two Atoms that form bond between each other
	 *  if found; null, otherwise.
	 */
	public static Atom[] findNearestNonNCAtomLinkage(final Group group1,
			final boolean isGroup1AminoAcid, final Group group2, final boolean isGroup2AminoAcid,
			List<String> potentialNamesOfAtomOnGroup1, List<String> potentialNamesOfAtomOnGroup2,
			double bondLengthTolerance) {
		List<Atom[]> linkages = findNonNCAtomLinkages(group1, isGroup1AminoAcid, group2, isGroup2AminoAcid,
				potentialNamesOfAtomOnGroup1, potentialNamesOfAtomOnGroup2, bondLengthTolerance);
		
		Atom[] ret = null;
		double minDistance = Double.POSITIVE_INFINITY;

		for (Atom[] linkage : linkages) {
			double distance;
			try {
				distance = Calc.getDistance(linkage[0], linkage[1]);
			} catch (StructureException e) {
				continue;
			}
			
			if (distance < minDistance) {
				minDistance = distance;
				ret = linkage;
			}
		}
		
		return ret;
	}
	
	/**
	 * Find non-N-C linkages between two groups within tolerance of bond length,
	 * from potential atoms.
	 * @param group1 the first {@link Group}.
	 * @param isGroup1AminoAcid true if group1 is an amino acid.
	 * @param group2 the second {@link Group}.
	 * @param isGroup2AminoAcid true if group2 is an amino acid.
	 * @param bondLengthTolerance bond length error tolerance.
	 * @return a list, each element of which is an array of two Atoms that form bond 
	 * between each other.
	 */
	public static List<Atom[]> findNonNCAtomLinkages(final Group group1, final boolean isGroup1AminoAcid,
			final Group group2, final boolean isGroup2AminoAcid, final double bondLengthTolerance) {
		return findNonNCAtomLinkages(group1, isGroup1AminoAcid, group2,
				isGroup2AminoAcid, null, null, bondLengthTolerance);
	}
	
	/**
	 * Find non-N-C linkages between two groups within tolerance of bond length,
	 * from potential atoms.
	 * @param group1 the first {@link Group}.
	 * @param isGroup1AminoAcid true if group1 is an amino acid.
	 * @param group2 the second {@link Group}.
	 * @param isGroup2AminoAcid true if group2 is an amino acid.
	 * @param potentialNamesOfAtomOnGroup1 potential names of the atom on the first group.
	 * 		  If null, search all atoms on the first group.
	 * @param potentialNamesOfAtomOnGroup2 potential names of the atom on the second group.
	 * 		  If null, search all atoms on the second group.
	 * @param bondLengthTolerance bond length error tolerance.
	 * @return a list, each element of which is an array of two Atoms that form bond 
	 * between each other.
	 */
	public static List<Atom[]> findNonNCAtomLinkages(final Group group1, 
			final boolean isGroup1AminoAcid,  final Group group2,
			final boolean isGroup2AminoAcid,
			List<String> potentialNamesOfAtomOnGroup1,
			List<String> potentialNamesOfAtomOnGroup2,
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
					if (isGroup1AminoAcid && isGroup2AminoAcid &&
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
		double distance;
		
		try {
			ret[0] = group1.getAtom(nameOfAtomOnGroup1);
			ret[1] = group2.getAtom(nameOfAtomOnGroup2);
			distance = Calc.getDistance(ret[0], ret[1]);
		} catch (StructureException e) {
			return null;
		}
		
		if (ret[0]==null || ret[1]==null) {
			return null;
		}
		
		float radiusOfAtom1 = ret[0].getElement().getCovalentRadius();
		float radiusOfAtom2 = ret[1].getElement().getCovalentRadius();
		
		if (Math.abs(distance-radiusOfAtom1 -radiusOfAtom2)
				> bondLengthTolerance) {
			return null;
		}
		
		return ret;
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
}
