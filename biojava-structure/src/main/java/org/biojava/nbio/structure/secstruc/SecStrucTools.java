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
package org.biojava.nbio.structure.secstruc;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupIterator;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;

/**
 * This class contains methods for obtaining and converting secondary structure
 * information from BioJava {@link Structure}s.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucTools {

	/**
	 * Obtain the List of secondary structure information (SecStrucInfo) of a
	 * Structure.
	 * 
	 * @param s
	 *            Structure with SS assignments
	 * @return List of SecStrucInfo objects
	 */
	public static List<SecStrucInfo> getSecStrucInfo(Structure s) {

		List<SecStrucInfo> listSSI = new ArrayList<SecStrucInfo>();
		GroupIterator iter = new GroupIterator(s);

		while (iter.hasNext()) {
			Group g = iter.next();
			if (g.hasAminoAtoms()) {
				Object p = g.getProperty(Group.SEC_STRUC);
				if (!(p == null)) {
					SecStrucInfo ss = (SecStrucInfo) p;
					listSSI.add(ss);
				}
			}
		}

		return listSSI;
	}

	/**
	 * Obtain the List of secondary structure elements (SecStrucElement) of a
	 * Structure.
	 * 
	 * @param s
	 *            Structure with SS assignments
	 * @return List of SecStrucElement objects
	 */
	public static List<SecStrucElement> getSecStrucElements(Structure s) {

		List<SecStrucElement> listSSE = new ArrayList<SecStrucElement>();
		GroupIterator iter = new GroupIterator(s);

		// SecStruc information - initialize
		SecStrucType type = SecStrucType.coil;
		ResidueNumber previous = new ResidueNumber();
		ResidueNumber start = new ResidueNumber();
		String chainId = "";
		int count = 0; // counts the number of residues in SSE

		// Create a map for the IDs of the SSE in the structure
		Map<SecStrucType, Integer> ids = new TreeMap<SecStrucType, Integer>();
		for (SecStrucType t : SecStrucType.values())
			ids.put(t, 1);

		while (iter.hasNext()) {
			Group g = iter.next();

			if (g.hasAminoAtoms()) {
				Object p = g.getProperty(Group.SEC_STRUC);
				if (p == null)
					continue;
				SecStrucInfo ss = (SecStrucInfo) p;

				if (count > 0) {
					// If chain and type are equal increment counter
					if (ss.type == type && chainId == g.getChainId()) {
						previous = g.getResidueNumber();
						count++;
						continue;
					} else {
						// Save the current SSE if chain or type change
						SecStrucElement sse = new SecStrucElement(type, start,
								previous, count, ids.get(type), chainId);
						listSSE.add(sse);
						ids.put(type, ids.get(type) + 1);
						count = 0;

						// Initialize a new SSE one
						if (ss.type != SecStrucType.coil) {
							type = ss.type;
							start = g.getResidueNumber();
							previous = start;
							chainId = g.getChainId();
							count = 1;
						}
					}
				} else {
					// This is for the first residue only
					if (ss.type != SecStrucType.coil) {
						type = ss.type;
						start = g.getResidueNumber();
						previous = start;
						chainId = g.getChainId();
						count = 1;
					}
				}
			}
		}
		return listSSE;
	}

	/**
	 * Obtain the List of secondary structure elements (SecStrucElement) of a
	 * List of Groups (assumed to be sequential, this is, connected in the
	 * original Structure).
	 * 
	 * @param groups
	 *            Structure with SS assignments
	 * @return List of SecStrucElement objects
	 */
	public static List<SecStrucElement> getSecStrucElements(List<Group> groups) {

		List<SecStrucElement> listSSE = new ArrayList<SecStrucElement>();

		// SecStruc information - initialize
		SecStrucType type = SecStrucType.coil;
		ResidueNumber previous = new ResidueNumber();
		ResidueNumber start = new ResidueNumber();
		String chainId = "";
		int count = 0; // counts the number of residues in SSE

		// Create a map for the IDs of the SSE in the structure
		Map<SecStrucType, Integer> ids = new TreeMap<SecStrucType, Integer>();
		for (SecStrucType t : SecStrucType.values())
			ids.put(t, 1);

		for (Group g : groups) {

			if (g.hasAminoAtoms()) {
				Object p = g.getProperty(Group.SEC_STRUC);
				if (p == null)
					continue;
				SecStrucInfo ss = (SecStrucInfo) p;

				if (count > 0) {
					// If chain and type are equal increment counter
					if (ss.type == type && chainId == g.getChainId()) {
						previous = g.getResidueNumber();
						count++;
						continue;
					} else {
						// Save the current SSE if chain or type change
						SecStrucElement sse = new SecStrucElement(type, start,
								previous, count, ids.get(type), chainId);
						listSSE.add(sse);
						ids.put(type, ids.get(type) + 1);
						count = 0;

						// Initialize a new SSE one
						if (ss.type != SecStrucType.coil) {
							type = ss.type;
							start = g.getResidueNumber();
							previous = start;
							chainId = g.getChainId();
							count = 1;
						}
					}
				} else {
					// This is for the first residue only
					if (ss.type != SecStrucType.coil) {
						type = ss.type;
						start = g.getResidueNumber();
						previous = start;
						chainId = g.getChainId();
						count = 1;
					}
				}
			}
		}
		return listSSE;
	}
}
