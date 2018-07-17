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
 * Created on Jan 4, 2006
 *
 */
package org.biojava.nbio.structure;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.Grid;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.core.util.FileDownloadUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A class that provides some tool methods.
 *
 * @author Andreas Prlic
 * @author Jules Jacobsen
 * @since 1.0
 */
public class StructureTools {

	private static final Logger logger = LoggerFactory
			.getLogger(StructureTools.class);

	// Amino Acid backbone
	/**
	 * The atom name of the backbone C-alpha atom. Note that this can be
	 * ambiguous depending on the context since Calcium atoms use the same name
	 * in PDB.
	 */
	public static final String CA_ATOM_NAME = "CA";

	/**
	 * The atom name for the backbone amide nitrogen
	 */
	public static final String N_ATOM_NAME = "N";

	/**
	 * The atom name for the backbone carbonyl
	 */
	public static final String C_ATOM_NAME = "C";

	/**
	 * The atom name for the backbone carbonyl oxygen
	 */
	public static final String O_ATOM_NAME = "O";

	/**
	 * The atom name of the side-chain C-beta atom
	 */
	public static final String CB_ATOM_NAME = "CB";

	// Nucleotide backbone
	/**
	 * The atom name of the backbone C1' in RNA
	 */
	public static final String C1_ATOM_NAME = "C1'";
	/**
	 * The atom name of the backbone C2' in RNA
	 */
	public static final String C2_ATOM_NAME = "C2'";
	/**
	 * The atom name of the backbone C3' in RNA
	 */
	public static final String C3_ATOM_NAME = "C3'";
	/**
	 * The atom name of the backbone C4' in RNA
	 */
	public static final String C4_ATOM_NAME = "C4'";
	/**
	 * The atom name of the backbone O2' in RNA
	 */
	public static final String O2_ATOM_NAME = "O2'";
	/**
	 * The atom name of the backbone O3' in RNA
	 */
	public static final String O3_ATOM_NAME = "O3'";
	/**
	 * The atom name of the backbone O4' in RNA
	 */
	public static final String O4_ATOM_NAME = "O4'";
	/**
	 * The atom name of the backbone O4' in RNA
	 */
	public static final String O5_ATOM_NAME = "O5'";
	/**
	 * The atom name of the backbone O4' in RNA
	 */
	public static final String OP1_ATOM_NAME = "OP1";
	/**
	 * The atom name of the backbone O4' in RNA
	 */
	public static final String OP2_ATOM_NAME = "OP2";
	/**
	 * The atom name of the backbone phosphate in RNA
	 */
	public static final String P_ATOM_NAME = "P";

	/**
	 * The atom used as representative for nucleotides, equivalent to
	 * {@link #CA_ATOM_NAME} for proteins
	 */
	public static final String NUCLEOTIDE_REPRESENTATIVE = P_ATOM_NAME;

	/**
	 * The character to use for unknown compounds in sequence strings
	 */
	public static final char UNKNOWN_GROUP_LABEL = 'X';

	/**
	 * Below this ratio of aminoacid/nucleotide residues to the sequence total,
	 * we use simple majority of aminoacid/nucleotide residues to decide the
	 * character of the chain (protein/nucleotide)
	 */
	public static final double RATIO_RESIDUES_TO_TOTAL = 0.95;

	/**
	 * Threshold for plausible binding of a ligand to the selected substructure
	 */
	public static final double DEFAULT_LIGAND_PROXIMITY_CUTOFF = 5;

	// there is a file format change in PDB 3.0 and nucleotides are being
	// renamed
	private static final Map<String, Character> nucleotides30;
	private static final Map<String, Character> nucleotides23;

	// amino acid 3 and 1 letter code definitions
	private static final Map<String, Character> aminoAcids;

	private static final Set<Element> hBondDonorAcceptors;

	static {
		nucleotides30 = new HashMap<String, Character>();
		nucleotides30.put("DA", 'A');
		nucleotides30.put("DC", 'C');
		nucleotides30.put("DG", 'G');
		nucleotides30.put("DT", 'T');
		nucleotides30.put("DI", 'I');
		nucleotides30.put("A", 'A');
		nucleotides30.put("G", 'G');
		nucleotides30.put("C", 'C');
		nucleotides30.put("U", 'U');
		nucleotides30.put("I", 'I');

		// the DNA linkers - the +C , +G, +A +T +U and +I have been replaced
		// with these:
		nucleotides30.put("TAF", UNKNOWN_GROUP_LABEL); // Fluorinated Thymine
		nucleotides30.put("TC1", UNKNOWN_GROUP_LABEL); // Furanosyl
		nucleotides30.put("TFE", UNKNOWN_GROUP_LABEL); // Fluorinated Thymine
		nucleotides30.put("TFO", UNKNOWN_GROUP_LABEL); // Tenofovir (3'
		// terminator)
		nucleotides30.put("TGP", UNKNOWN_GROUP_LABEL); // Guanine variant
		nucleotides30.put("THX", UNKNOWN_GROUP_LABEL); // 5' terminator
		nucleotides30.put("TLC", UNKNOWN_GROUP_LABEL); // Thymine with dicyclic
		// sugar
		nucleotides30.put("TLN", UNKNOWN_GROUP_LABEL); // locked Thymine
		nucleotides30.put("LCG", UNKNOWN_GROUP_LABEL); // locked Guanine
		nucleotides30.put("TP1", UNKNOWN_GROUP_LABEL); // Thymine peptide
		// nucleic acid, with
		// added methyl
		nucleotides30.put("CP1", UNKNOWN_GROUP_LABEL); // Cytidine peptide
		// nucleic acid, with
		// added methyl
		nucleotides30.put("TPN", UNKNOWN_GROUP_LABEL); // Thymine peptide
		// nucleic acid
		nucleotides30.put("CPN", UNKNOWN_GROUP_LABEL); // Cytidine peptide
		// nucleic acid
		nucleotides30.put("GPN", UNKNOWN_GROUP_LABEL); // Guanine peptide
		// nucleic acid
		nucleotides30.put("APN", UNKNOWN_GROUP_LABEL); // Adenosine peptide
		// nucleic acid
		nucleotides30.put("TPC", UNKNOWN_GROUP_LABEL); // Thymine variant

		// store nucleic acids (C, G, A, T, U, and I), and
		// the modified versions of nucleic acids (+C, +G, +A, +T, +U, and +I),
		// and
		nucleotides23 = new HashMap<String, Character>();
		String[] names = { "C", "G", "A", "T", "U", "I", "+C", "+G", "+A",
				"+T", "+U", "+I" };
		for (String n : names) {
			nucleotides23.put(n, n.charAt(n.length() - 1));
		}

		aminoAcids = new HashMap<String, Character>();
		aminoAcids.put("GLY", 'G');
		aminoAcids.put("ALA", 'A');
		aminoAcids.put("VAL", 'V');
		aminoAcids.put("LEU", 'L');
		aminoAcids.put("ILE", 'I');
		aminoAcids.put("PHE", 'F');
		aminoAcids.put("TYR", 'Y');
		aminoAcids.put("TRP", 'W');
		aminoAcids.put("PRO", 'P');
		aminoAcids.put("HIS", 'H');
		aminoAcids.put("LYS", 'K');
		aminoAcids.put("ARG", 'R');
		aminoAcids.put("SER", 'S');
		aminoAcids.put("THR", 'T');
		aminoAcids.put("GLU", 'E');
		aminoAcids.put("GLN", 'Q');
		aminoAcids.put("ASP", 'D');
		aminoAcids.put("ASN", 'N');
		aminoAcids.put("CYS", 'C');
		aminoAcids.put("MET", 'M');
		// MSE is only found as a molecular replacement for MET
		aminoAcids.put("MSE", 'M');
		// 'non-standard', genetically encoded
		// http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html
		// IUBMB recommended name is 'SEC' but the wwPDB currently use 'CSE'
		// likewise 'PYL' (IUBMB) and 'PYH' (PDB)
		aminoAcids.put("CSE", 'U');
		aminoAcids.put("SEC", 'U');
		aminoAcids.put("PYH", 'O');
		aminoAcids.put("PYL", 'O');

		hBondDonorAcceptors = new HashSet<Element>();
		hBondDonorAcceptors.add(Element.N);
		hBondDonorAcceptors.add(Element.O);
		hBondDonorAcceptors.add(Element.S);

	}

	/**
	 * Count how many Atoms are contained within a Structure object.
	 *
	 * @param s
	 *            the structure object
	 * @return the number of Atoms in this Structure
	 */
	public static final int getNrAtoms(Structure s) {

		int nrAtoms = 0;

		Iterator<Group> iter = new GroupIterator(s);

		while (iter.hasNext()) {
			Group g = iter.next();
			nrAtoms += g.size();
		}

		return nrAtoms;
	}

	/**
	 * Count how many groups are contained within a structure object.
	 *
	 * @param s
	 *            the structure object
	 * @return the number of groups in the structure
	 */
	public static final int getNrGroups(Structure s) {
		int nrGroups = 0;

		List<Chain> chains = s.getChains(0);
		for (Chain c : chains) {
			nrGroups += c.getAtomLength();
		}
		return nrGroups;
	}

	/**
	 * Returns an array of the requested Atoms from the Structure object.
	 * Iterates over all groups and checks if the requested atoms are in this
	 * group, no matter if this is a {@link AminoAcid} or {@link HetatomImpl}
	 * group. If the group does not contain all requested atoms then no atoms
	 * are added for that group. For structures with more than one model, only
	 * model 0 will be used.
	 *
	 * @param s
	 *            the structure to get the atoms from
	 *
	 * @param atomNames
	 *            contains the atom names to be used.
	 * @return an Atom[] array
	 */
	public static final Atom[] getAtomArray(Structure s, String[] atomNames) {
		List<Chain> chains = s.getModel(0);

		List<Atom> atoms = new ArrayList<Atom>();

		extractAtoms(atomNames, chains, atoms);
		return atoms.toArray(new Atom[atoms.size()]);

	}

	/**
	 * Returns an array of the requested Atoms from the Structure object. In
	 * contrast to {@link #getAtomArray(Structure, String[])} this method
	 * iterates over all chains. Iterates over all chains and groups and checks
	 * if the requested atoms are in this group, no matter if this is a
	 * {@link AminoAcid} or {@link HetatomImpl} group. If the group does not
	 * contain all requested atoms then no atoms are added for that group. For
	 * structures with more than one model, only model 0 will be used.
	 *
	 * @param s
	 *            the structure to get the atoms from
	 *
	 * @param atomNames
	 *            contains the atom names to be used.
	 * @return an Atom[] array
	 */
	public static final Atom[] getAtomArrayAllModels(Structure s,
			String[] atomNames) {

		List<Atom> atoms = new ArrayList<Atom>();

		for (int i = 0; i < s.nrModels(); i++) {
			List<Chain> chains = s.getModel(i);
			extractAtoms(atomNames, chains, atoms);
		}
		return atoms.toArray(new Atom[atoms.size()]);

	}

	/**
	 * Convert all atoms of the structure (all models) into an Atom array
	 *
	 * @param s
	 *            input structure
	 * @return all atom array
	 */
	public static final Atom[] getAllAtomArray(Structure s) {
		List<Atom> atoms = new ArrayList<Atom>();

		AtomIterator iter = new AtomIterator(s);
		while (iter.hasNext()) {
			Atom a = iter.next();
			atoms.add(a);
		}
		return atoms.toArray(new Atom[atoms.size()]);
	}
	/**
	 * Convert all atoms of the structure (specified model) into an Atom array
	 *
	 * @param s
	 *            input structure
	 * @return all atom array
	 */
	public static final Atom[] getAllAtomArray(Structure s, int model) {
		List<Atom> atoms = new ArrayList<Atom>();

		AtomIterator iter = new AtomIterator(s,model);
		while (iter.hasNext()) {
			Atom a = iter.next();
			atoms.add(a);
		}
		return atoms.toArray(new Atom[atoms.size()]);

	}

	/**
	 * Returns and array of all atoms of the chain, including
	 * Hydrogens (if present) and all HETATOMs. Waters are not included.
	 *
	 * @param c
	 *            input chain
	 * @return all atom array
	 */
	public static final Atom[] getAllAtomArray(Chain c) {
		List<Atom> atoms = new ArrayList<Atom>();

		for (Group g : c.getAtomGroups()) {
			if (g.isWater())
				continue;
			for (Atom a : g.getAtoms()) {
				atoms.add(a);
			}
		}
		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * List of groups from the structure not included in ca (e.g. ligands).
	 *
	 * Unaligned groups are searched from all chains referenced in ca, as well
	 * as any chains in the first model of the structure from ca[0], if any.
	 *
	 * @param ca an array of atoms
	 * @return
	 */
	public static List<Group> getUnalignedGroups(Atom[] ca) {
		Set<Chain> chains = new HashSet<Chain>();
		Set<Group> caGroups = new HashSet<Group>();

		// Create list of all chains in this structure
		Structure s = null;
		if (ca.length > 0) {
			Group g = ca[0].getGroup();
			if (g != null) {
				Chain c = g.getChain();
				if (c != null) {
					s = c.getStructure();
				}
			}
		}
		if (s != null) {
			// Add all chains from the structure
			for (Chain c : s.getChains(0)) {
				chains.add(c);
			}
		}

		// Add groups and chains from ca
		for (Atom a : ca) {
			Group g = a.getGroup();
			if (g != null) {
				caGroups.add(g);

				Chain c = g.getChain();
				if (c != null) {
					chains.add(c);
				}
			}
		}

		// Iterate through all chains, finding groups not in ca
		List<Group> unadded = new ArrayList<Group>();
		for (Chain c : chains) {
			for (Group g : c.getAtomGroups()) {
				if (!caGroups.contains(g)) {
					unadded.add(g);
				}
			}
		}
		return unadded;
	}

	/**
	 * Finds all ligand groups from the target which fall within the cutoff distance
	 * of some atom from the query set.
	 * 
	 * @param target Set of groups including the ligands
	 * @param query Atom selection
	 * @param cutoff Distance from query atoms to consider, in angstroms
	 * @return All groups from the target with at least one atom within cutoff of a query atom
	 * @see StructureTools#DEFAULT_LIGAND_PROXIMITY_CUTOFF
	 */
	public static List<Group> getLigandsByProximity(Collection<Group> target, Atom[] query, double cutoff) {
		// Geometric hashing of the reduced structure
		Grid grid = new Grid(cutoff);
		grid.addAtoms(query);

		List<Group> ligands = new ArrayList<>();
		for(Group g :target ) {
			// don't worry about waters
			if(g.isWater()) {
				continue;
			}

			if(g.isPolymeric() ) {
				// Polymers aren't ligands
				continue;
			}

			// It is a ligand!

			// Check that it's within cutoff of something in reduced
			List<Atom> groupAtoms = g.getAtoms();
			if( ! grid.hasAnyContact(Calc.atomsToPoints(groupAtoms))) {
				continue;
			}

			ligands.add(g);
		}
		return ligands;
	}
	
	/**
	 * Adds a particular group to a structure. A new chain will be created if necessary.
	 * 
	 * <p>When adding multiple groups, pass the return value of one call as the
	 * chainGuess parameter of the next call for efficiency.
	 * <pre>
	 * Chain guess = null;
	 * for(Group g : groups) {
	 *     guess = addGroupToStructure(s, g, guess );
	 * }
	 * </pre>
	 * @param s structure to receive the group
	 * @param g group to add
	 * @param chainGuess (optional) If not null, should be a chain from s. Used
	 *  to improve performance when adding many groups from the same chain
	 * @param clone Indicates whether the input group should be cloned before
	 *  being added to the new chain
	 * @return the chain g was added to
	 */
	public static Chain addGroupToStructure(Structure s, Group g, int model, Chain chainGuess, boolean clone ) {
		synchronized(s) {
			// Find or create the chain
			String chainId = g.getChainId();
			assert !chainId.isEmpty();
			Chain chain;
			if(chainGuess != null && chainGuess.getId() == chainId) {
				// previously guessed chain
				chain = chainGuess;
			} else {
				// Try to guess
				chain = s.getChain(chainId, model);
				if(chain == null) {
					// no chain found
					chain = new ChainImpl();
					chain.setId(chainId);

					Chain oldChain = g.getChain();
					chain.setName(oldChain.getName());

					EntityInfo oldEntityInfo = oldChain.getEntityInfo();

					EntityInfo newEntityInfo = s.getEntityById(oldEntityInfo.getMolId());
					if( newEntityInfo == null ) {
						newEntityInfo = new EntityInfo(oldEntityInfo);
						s.addEntityInfo(newEntityInfo);
					}
					newEntityInfo.addChain(chain);
					chain.setEntityInfo(newEntityInfo);
					
					// TODO Do the seqres need to be cloned too? -SB 2016-10-7
					chain.setSeqResGroups(oldChain.getSeqResGroups());
					chain.setSeqMisMatches(oldChain.getSeqMisMatches());
					
					s.addChain(chain,model);
				}
			}

			// Add cloned group
			if(clone) {
				g = (Group)g.clone();
			}
			chain.addGroup(g);

			return chain;
		}
	}

	/**
	 * Add a list of groups to a new structure. Chains will be automatically
	 * created in the new structure as needed.
	 * @param s structure to receive the group
	 * @param g group to add
	 * @param clone Indicates whether the input groups should be cloned before
	 *  being added to the new chain
	 */
	public static void addGroupsToStructure(Structure s, Collection<Group> groups, int model, boolean clone) {
		Chain chainGuess = null;
		for(Group g : groups) {
			chainGuess = addGroupToStructure(s, g, model, chainGuess, clone);
		}
	}
	
	/**
	 * Expand a set of atoms into all groups from the same structure.
	 * 
	 * If the structure is set, only the first atom is used (assuming all
	 * atoms come from the same original structure).
	 * If the atoms aren't linked to a structure (for instance, for cloned atoms),
	 * searches all chains of all atoms for groups.
	 * @param atoms Sample of atoms
	 * @return All groups from all chains accessible from the input atoms
	 */
	public static Set<Group> getAllGroupsFromSubset(Atom[] atoms) {
		return getAllGroupsFromSubset(atoms,null);
	}
	/**
	 * Expand a set of atoms into all groups from the same structure.
	 * 
	 * If the structure is set, only the first atom is used (assuming all
	 * atoms come from the same original structure).
	 * If the atoms aren't linked to a structure (for instance, for cloned atoms),
	 * searches all chains of all atoms for groups.
	 * @param atoms Sample of atoms
	 * @param types Type of groups to return (useful for getting only ligands, for instance).
	 *  Null gets all groups.
	 * @return All groups from all chains accessible from the input atoms
	 */
	public static Set<Group> getAllGroupsFromSubset(Atom[] atoms,GroupType types) {
		// Get the full structure
		Structure s = null;
		if (atoms.length > 0) {
			Group g = atoms[0].getGroup();
			if (g != null) {
				Chain c = g.getChain();
				if (c != null) {
					s = c.getStructure();
				}
			}
		}
		// Collect all groups from the structure
		Set<Chain> allChains = new HashSet<>();
		if( s != null ) {
			allChains.addAll(s.getChains());
		}
		// In case the structure wasn't set, need to use ca chains too
		for(Atom a : atoms) {
			Group g = a.getGroup();
			if(g != null) {
				Chain c = g.getChain();
				if( c != null ) {
					allChains.add(c);
				}
			}
		}

		if(allChains.isEmpty() ) {
			return Collections.emptySet();
		}
		
		// Extract all ligand groups
		Set<Group> full = new HashSet<>();
		for(Chain c : allChains) {
			if(types == null) {
				full.addAll(c.getAtomGroups());
			} else {
				full.addAll(c.getAtomGroups(types));
			}
		}

		return full;
	}


	/**
	 * Returns and array of all non-Hydrogen atoms in the given Structure,
	 * optionally including HET atoms or not. Waters are not included.
	 *
	 * @param s
	 * @param hetAtoms
	 *            if true HET atoms are included in array, if false they are not
	 * @return
	 */
	public static final Atom[] getAllNonHAtomArray(Structure s, boolean hetAtoms) {
		AtomIterator iter = new AtomIterator(s);
		return getAllNonHAtomArray(s, hetAtoms, iter);
	}
	/**
	 * Returns and array of all non-Hydrogen atoms in the given Structure,
	 * optionally including HET atoms or not. Waters are not included.
	 *
	 * @param s
	 * @param hetAtoms
	 *            if true HET atoms are included in array, if false they are not
	 * @param modelNr Model number to draw atoms from
	 * @return
	 */
	public static final Atom[] getAllNonHAtomArray(Structure s, boolean hetAtoms, int modelNr) {
		AtomIterator iter = new AtomIterator(s,modelNr);
		return getAllNonHAtomArray(s, hetAtoms, iter);
	}
	private static final Atom[] getAllNonHAtomArray(Structure s, boolean hetAtoms, AtomIterator iter) {
		List<Atom> atoms = new ArrayList<Atom>();

		while (iter.hasNext()) {
			Atom a = iter.next();
			if (a.getElement() == Element.H)
				continue;

			Group g = a.getGroup();

			if (g.isWater())
				continue;

			if (!hetAtoms && g.getType().equals(GroupType.HETATM))
				continue;

			atoms.add(a);
		}
		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * Returns and array of all non-Hydrogen atoms in the given Chain,
	 * optionally including HET atoms or not Waters are not included.
	 *
	 * @param c
	 * @param hetAtoms
	 *            if true HET atoms are included in array, if false they are not
	 * @return
	 */
	public static final Atom[] getAllNonHAtomArray(Chain c, boolean hetAtoms) {
		List<Atom> atoms = new ArrayList<Atom>();

		for (Group g : c.getAtomGroups()) {
			if (g.isWater())
				continue;
			for (Atom a : g.getAtoms()) {

				if (a.getElement() == Element.H)
					continue;

				if (!hetAtoms && g.getType().equals(GroupType.HETATM))
					continue;

				atoms.add(a);
			}
		}
		return atoms.toArray(new Atom[atoms.size()]);
	}
	
	/**
	 * Returns and array of all non-Hydrogen atoms coordinates in the given Chain,
	 * optionally including HET atoms or not Waters are not included.
	 *
	 * @param c
	 * @param hetAtoms
	 *            if true HET atoms are included in array, if false they are not
	 * @return
	 */
	public static final Point3d[] getAllNonHCoordsArray(Chain c, boolean hetAtoms) {
		List<Point3d> atoms = new ArrayList<Point3d>();

		for (Group g : c.getAtomGroups()) {
			if (g.isWater())
				continue;
			for (Atom a : g.getAtoms()) {

				if (a.getElement() == Element.H)
					continue;

				if (!hetAtoms && g.getType().equals(GroupType.HETATM))
					continue;

				atoms.add(a.getCoordsAsPoint3d());
			}
		}
		return atoms.toArray(new Point3d[atoms.size()]);
	}

	/**
	 * Adds to the given atoms list, all atoms of groups that contained all
	 * requested atomNames, i.e. if a group does not contain all of the
	 * requested atom names, its atoms won't be added.
	 *
	 * @param atomNames
	 * @param chains
	 * @param atoms
	 */
	private static void extractAtoms(String[] atomNames, List<Chain> chains,
			List<Atom> atoms) {

		for (Chain c : chains) {

			for (Group g : c.getAtomGroups()) {

				// a temp container for the atoms of this group
				List<Atom> thisGroupAtoms = new ArrayList<Atom>();
				// flag to check if this group contains all the requested atoms.
				boolean thisGroupAllAtoms = true;
				for (String atomName : atomNames) {
					Atom a = g.getAtom(atomName);

					if (a == null) {
						// this group does not have a required atom, skip it...
						thisGroupAllAtoms = false;
						break;
					}
					thisGroupAtoms.add(a);
				}
				if (thisGroupAllAtoms) {
					// add the atoms of this group to the array.
					for (Atom a : thisGroupAtoms) {
						atoms.add(a);
					}
				}

			}
		}
	}

	/**
	 * Returns an array of the requested Atoms from the Chain object. Iterates
	 * over all groups and checks if the requested atoms are in this group, no
	 * matter if this is a AminoAcid or Hetatom group. If the group does not
	 * contain all requested atoms then no atoms are added for that group.
	 *
	 * @param c
	 *            the Chain to get the atoms from
	 *
	 * @param atomNames
	 *            contains the atom names to be used.
	 * @return an Atom[] array
	 */
	public static final Atom[] getAtomArray(Chain c, String[] atomNames) {

		List<Atom> atoms = new ArrayList<Atom>();

		for (Group g : c.getAtomGroups()) {

			// a temp container for the atoms of this group
			List<Atom> thisGroupAtoms = new ArrayList<Atom>();
			// flag to check if this group contains all the requested atoms.
			boolean thisGroupAllAtoms = true;
			for (String atomName : atomNames) {
				Atom a = g.getAtom(atomName);
				if (a == null) {
					logger.debug("Group " + g.getResidueNumber() + " ("
							+ g.getPDBName()
							+ ") does not have the required atom '" + atomName
							+ "'");
					// this group does not have a required atom, skip it...
					thisGroupAllAtoms = false;
					break;
				}
				thisGroupAtoms.add(a);
			}

			if (thisGroupAllAtoms) {
				// add the atoms of this group to the array.
				for (Atom a : thisGroupAtoms) {
					atoms.add(a);
				}
			}

		}
		return atoms.toArray(new Atom[atoms.size()]);

	}

	/**
	 * Returns an Atom array of the C-alpha atoms. Any atom that is a carbon and
	 * has CA name will be returned.
	 *
	 * @param c
	 *            the structure object
	 * @return an Atom[] array
	 * @see #getRepresentativeAtomArray(Chain)
	 */
	public static final Atom[] getAtomCAArray(Chain c) {
		List<Atom> atoms = new ArrayList<Atom>();

		for (Group g : c.getAtomGroups()) {
			if (g.hasAtom(CA_ATOM_NAME)
					&& g.getAtom(CA_ATOM_NAME).getElement() == Element.C) {
				atoms.add(g.getAtom(CA_ATOM_NAME));
			}
		}

		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * Gets a representative atom for each group that is part of the chain
	 * backbone. Note that modified aminoacids won't be returned as part of the
	 * backbone if the {@link org.biojava.nbio.structure.io.mmcif.ReducedChemCompProvider} was used to load the
	 * structure.
	 *
	 * For amino acids, the representative is a CA carbon. For nucleotides, the
	 * representative is the {@value #NUCLEOTIDE_REPRESENTATIVE}. Other group
	 * types will be ignored.
	 *
	 * @param c
	 * @return representative Atoms of the chain backbone
	 * @since Biojava 4.1.0
	 */
	public static final Atom[] getRepresentativeAtomArray(Chain c) {
		List<Atom> atoms = new ArrayList<Atom>();

		for (Group g : c.getAtomGroups()) {

			switch (g.getType()) {
			case AMINOACID:
				if (g.hasAtom(CA_ATOM_NAME)
						&& g.getAtom(CA_ATOM_NAME).getElement() == Element.C) {
					atoms.add(g.getAtom(CA_ATOM_NAME));
				}
				break;
			case NUCLEOTIDE:
				if (g.hasAtom(NUCLEOTIDE_REPRESENTATIVE)) {
					atoms.add(g.getAtom(NUCLEOTIDE_REPRESENTATIVE));
				}
				break;
			default:
				// don't add
			}
		}

		return atoms.toArray(new Atom[atoms.size()]);

	}

	/**
	 * Provides an equivalent copy of Atoms in a new array. Clones everything,
	 * starting with parent groups and chains. The chain will only contain
	 * groups that are part of the input array.
	 *
	 * @param ca
	 *            array of representative atoms, e.g. CA atoms
	 * @return Atom array
	 * @since Biojava 4.1.0
	 */
	public static final Atom[] cloneAtomArray(Atom[] ca) {
		Atom[] newCA = new Atom[ca.length];

		List<Chain> model = new ArrayList<Chain>();
		int apos = -1;
		for (Atom a : ca) {
			apos++;
			Group parentG = a.getGroup();
			Chain parentC = parentG.getChain();

			Chain newChain = null;
			for (Chain c : model) {
				if (c.getName().equals(parentC.getName())) {
					newChain = c;
					break;
				}
			}
			if (newChain == null) {
				newChain = new ChainImpl();
				newChain.setId(parentC.getId());
				newChain.setName(parentC.getName());
				model.add(newChain);
			}

			Group parentN = (Group) parentG.clone();

			newCA[apos] = parentN.getAtom(a.getName());
			try {
				// if the group doesn't exist yet, this produces a StructureException
				newChain.getGroupByPDB(parentN.getResidueNumber());
			} catch (StructureException e) {
				// the group doesn't exist yet in the newChain, let's add it
				newChain.addGroup(parentN);
			}

		}
		return newCA;
	}

	/**
	 * Clone a set of representative Atoms, but returns the parent groups
	 *
	 * @param ca
	 *            Atom array
	 * @return Group array
	 */
	public static Group[] cloneGroups(Atom[] ca) {
		Group[] newGroup = new Group[ca.length];

		List<Chain> model = new ArrayList<Chain>();
		int apos = -1;
		for (Atom a : ca) {
			apos++;
			Group parentG = a.getGroup();
			Chain parentC = parentG.getChain();

			Chain newChain = null;
			for (Chain c : model) {
				if (c.getName().equals(parentC.getName())) {
					newChain = c;
					break;
				}
			}
			if (newChain == null) {
				newChain = new ChainImpl();
				newChain.setName(parentC.getName());
				model.add(newChain);
			}

			Group ng = (Group) parentG.clone();
			newGroup[apos] = ng;
			newChain.addGroup(ng);
		}
		return newGroup;
	}

	/**
	 * Utility method for working with circular permutations. Creates a
	 * duplicated and cloned set of Calpha atoms from the input array.
	 *
	 * @param ca2
	 *            atom array
	 * @return cloned and duplicated set of input array
	 */
	public static Atom[] duplicateCA2(Atom[] ca2) {
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length * 2];

		int pos = 0;

		Chain c = null;
		String prevChainId = "";
		for (Atom a : ca2) {
			Group g = (Group) a.getGroup().clone(); // works because each group
			// has only a single atom

			if (c == null) {
				c = new ChainImpl();
				Chain orig = a.getGroup().getChain();
				c.setId(orig.getId());
				c.setName(orig.getName());
			} else {
				Chain orig = a.getGroup().getChain();
				if (!orig.getId().equals(prevChainId)) {
					c = new ChainImpl();
					c.setId(orig.getId());
					c.setName(orig.getName());
				}
			}

			c.addGroup(g);
			ca2clone[pos] = g.getAtom(a.getName());

			pos++;
		}

		// Duplicate ca2!
		c = null;
		prevChainId = "";
		for (Atom a : ca2) {
			Group g = (Group) a.getGroup().clone();

			if (c == null) {
				c = new ChainImpl();
				Chain orig = a.getGroup().getChain();
				c.setId(orig.getId());
				c.setName(orig.getName());
			} else {
				Chain orig = a.getGroup().getChain();
				if (!orig.getId().equals(prevChainId)) {
					c = new ChainImpl();
					c.setId(orig.getId());
					c.setName(orig.getName());
				}
			}

			c.addGroup(g);
			ca2clone[pos] = g.getAtom(a.getName());

			pos++;
		}

		return ca2clone;

	}

	/**
	 * Return an Atom array of the C-alpha atoms. Any atom that is a carbon and
	 * has CA name will be returned.
	 *
	 * @param s
	 *            the structure object
	 * @return an Atom[] array
	 * @see #getRepresentativeAtomArray(Structure)
	 */
	public static Atom[] getAtomCAArray(Structure s) {

		List<Atom> atoms = new ArrayList<Atom>();

		for (Chain c : s.getChains()) {
			for (Group g : c.getAtomGroups()) {
				if (g.hasAtom(CA_ATOM_NAME)
						&& g.getAtom(CA_ATOM_NAME).getElement() == Element.C) {
					atoms.add(g.getAtom(CA_ATOM_NAME));
				}
			}
		}

		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * Gets a representative atom for each group that is part of the chain
	 * backbone. Note that modified aminoacids won't be returned as part of the
	 * backbone if the {@link org.biojava.nbio.structure.io.mmcif.ReducedChemCompProvider} was used to load the
	 * structure.
	 *
	 * For amino acids, the representative is a CA carbon. For nucleotides, the
	 * representative is the {@value #NUCLEOTIDE_REPRESENTATIVE}. Other group
	 * types will be ignored.
	 *
	 * @param s
	 *            Input structure
	 * @return representative Atoms of the structure backbone
	 * @since Biojava 4.1.0
	 */
	public static Atom[] getRepresentativeAtomArray(Structure s) {

		List<Atom> atoms = new ArrayList<Atom>();

		for (Chain c : s.getChains()) {
			Atom[] chainAtoms = getRepresentativeAtomArray(c);
			for (Atom a : chainAtoms) {
				atoms.add(a);
			}
		}

		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * Return an Atom array of the main chain atoms: CA, C, N, O Any group that
	 * contains those atoms will be included, be it a standard aminoacid or not
	 *
	 * @param s
	 *            the structure object
	 * @return an Atom[] array
	 */
	public static Atom[] getBackboneAtomArray(Structure s) {

		List<Atom> atoms = new ArrayList<Atom>();

		for (Chain c : s.getChains()) {
			for (Group g : c.getAtomGroups()) {
				if (g.hasAminoAtoms()) {
					// this means we will only take atoms grom groups that have
					// complete backbones
					for (Atom a : g.getAtoms()) {
						switch (g.getType()) {
						case NUCLEOTIDE:
							// Nucleotide backbone
							if (a.getName().equals(C1_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(C2_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(C3_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(C4_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(O2_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(O3_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(O4_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(O5_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(OP1_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(OP2_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(P_ATOM_NAME))
								atoms.add(a);
							// TODO Allow C4* names as well as C4'? -SB 3/2015
							break;
						case AMINOACID:
						default:
							// we do it this way instead of with g.getAtom() to
							// be sure we always use the same order as original
							if (a.getName().equals(CA_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(C_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(N_ATOM_NAME))
								atoms.add(a);
							if (a.getName().equals(O_ATOM_NAME))
								atoms.add(a);
							break;
						}
					}
				}
			}

		}

		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * Convert three character amino acid codes into single character e.g.
	 * convert CYS to C. Valid 3-letter codes will be those of the standard 20
	 * amino acids plus MSE, CSE, SEC, PYH, PYL (see the {@link #aminoAcids}
	 * map)
	 *
	 * @return the 1 letter code, or null if the given 3 letter code does not
	 *         correspond to an amino acid code
	 * @param groupCode3
	 *            a three character amino acid representation String
	 * @see {@link #get1LetterCode(String)}
	 */
	public static final Character get1LetterCodeAmino(String groupCode3) {
		return aminoAcids.get(groupCode3);
	}

	/**
	 * Convert a three letter amino acid or nucleotide code into a single
	 * character code. If the code does not correspond to an amino acid or
	 * nucleotide, returns {@link #UNKNOWN_GROUP_LABEL}.
	 *
	 * Returned null for nucleotides prior to version 4.0.1.
	 *
	 * @param groupCode3
	 *            three letter representation
	 * @return The 1-letter abbreviation
	 */
	public static final Character get1LetterCode(String groupCode3) {

		Character code1;

		// is it a standard amino acid ?
		code1 = get1LetterCodeAmino(groupCode3);

		if (code1 == null) {
			// hm groupCode3 is not standard
			// perhaps it is a nucleotide?
			groupCode3 = groupCode3.trim();
			if (isNucleotide(groupCode3)) {
				code1 = nucleotides30.get(groupCode3);
				if (code1 == null) {
					code1 = nucleotides23.get(groupCode3);
				}
				if (code1 == null) {
					code1 = UNKNOWN_GROUP_LABEL;
				}
			} else {
				// does not seem to be so let's assume it is
				// nonstandard aminoacid and label it "X"
				// logger.warning("unknown group name "+groupCode3 );
				code1 = UNKNOWN_GROUP_LABEL;
			}
		}

		return code1;

	}

	/**
	 * Test if the three-letter code of an ATOM entry corresponds to a
	 * nucleotide or to an aminoacid.
	 *
	 * @param groupCode3
	 *            3-character code for a group.
	 *
	 */
	public static final boolean isNucleotide(String groupCode3) {
		String code = groupCode3.trim();
		return nucleotides30.containsKey(code)
				|| nucleotides23.containsKey(code);
	}

	/**
	 * Reduce a structure to provide a smaller representation . Only takes the
	 * first model of the structure. If chainName is provided only return a
	 * structure containing that Chain ID. Converts lower case chain IDs to
	 * upper case if structure does not contain a chain with that ID.
	 *
	 * @param s
	 * @param chainId
	 * @return Structure
	 * @since 3.0
	 * @deprecated Use {@link StructureIdentifier#reduce(Structure)} instead (v. 4.2.0)
	 */
	@Deprecated
	public static final Structure getReducedStructure(Structure s,
			String chainId) throws StructureException {
		// since we deal here with structure alignments,
		// only use Model 1...

		Structure newS = new StructureImpl();
		newS.setPDBCode(s.getPDBCode());
		newS.setPDBHeader(s.getPDBHeader());
		newS.setName(s.getName());
		newS.setSSBonds(s.getSSBonds());
		newS.setDBRefs(s.getDBRefs());
		newS.setSites(s.getSites());
		newS.setBiologicalAssembly(s.isBiologicalAssembly());
		newS.setEntityInfos(s.getEntityInfos());
		newS.setSSBonds(s.getSSBonds());
		newS.setSites(s.getSites());

		if (chainId != null)
			chainId = chainId.trim();

		if (chainId == null || chainId.equals("")) {
			// only get model 0
			List<Chain> model0 = s.getModel(0);
			for (Chain c : model0) {
				newS.addChain(c);
			}
			return newS;

		}

		Chain c = null;
		try {
			c = s.getChainByPDB(chainId);
		} catch (StructureException e) {
			logger.warn(e.getMessage() + ". Chain id " + chainId
					+ " did not match, trying upper case Chain id.");
			c = s.getChainByPDB(chainId.toUpperCase());

		}
		if (c != null) {
			newS.addChain(c);
			for (EntityInfo comp : s.getEntityInfos()) {
				if (comp.getChainIds() != null
						&& comp.getChainIds().contains(c.getChainID())) {
					// found matching entity info. set description...
					newS.getPDBHeader().setDescription(
							"Chain " + c.getChainID() + " of " + s.getPDBCode()
							+ " " + comp.getDescription());
				}
			}
		}

		return newS;
	}

	public static final String convertAtomsToSeq(Atom[] atoms) {

		StringBuilder buf = new StringBuilder();
		Group prevGroup = null;
		for (Atom a : atoms) {
			Group g = a.getGroup();
			if (prevGroup != null) {
				if (prevGroup.equals(g)) {
					// we add each group only once.
					continue;
				}
			}
			String code3 = g.getPDBName();
			Character code1 = get1LetterCodeAmino(code3);
			if (code1 == null)
				code1 = UNKNOWN_GROUP_LABEL;

			buf.append(code1);

			prevGroup = g;

		}
		return buf.toString();
	}

	/**
	 * Get a group represented by a ResidueNumber.
	 *
	 * @param struc
	 *            a {@link Structure}
	 * @param pdbResNum
	 *            a {@link ResidueNumber}
	 * @return a group in the structure that is represented by the pdbResNum.
	 * @throws StructureException
	 *             if the group cannot be found.
	 */
	public static final Group getGroupByPDBResidueNumber(Structure struc,
			ResidueNumber pdbResNum) throws StructureException {
		if (struc == null || pdbResNum == null) {
			throw new IllegalArgumentException("Null argument(s).");
		}

		Chain chain = struc.getPolyChainByPDB(pdbResNum.getChainName());

		return chain.getGroupByPDB(pdbResNum);
	}

	/**
	 * Returns the set of intra-chain contacts for the given chain for given
	 * atom names, i.e. the contact map. Uses a geometric hashing algorithm that
	 * speeds up the calculation without need of full distance matrix. The
	 * parsing mode {@link FileParsingParameters#setAlignSeqRes(boolean)} needs
	 * to be set to true for this to work.
	 *
	 * @param chain
	 * @param atomNames
	 *            the array with atom names to be used. Beware: CA will do both
	 *            C-alphas an Calciums if null all non-H atoms of non-hetatoms
	 *            will be used
	 * @param cutoff
	 * @return
	 */
	public static AtomContactSet getAtomsInContact(Chain chain,
			String[] atomNames, double cutoff) {
		Grid grid = new Grid(cutoff);

		Atom[] atoms = null;
		if (atomNames == null) {
			atoms = getAllNonHAtomArray(chain, false);
		} else {
			atoms = getAtomArray(chain, atomNames);
		}
		// If tha
		if(atoms.length==0){ 
			logger.warn("No atoms found for buidling grid!");
			return new AtomContactSet(cutoff);
		}
		grid.addAtoms(atoms);

		return grid.getAtomContacts();
	}

	/**
	 * Returns the set of intra-chain contacts for the given chain for all non-H
	 * atoms of non-hetatoms, i.e. the contact map. Uses a geometric hashing
	 * algorithm that speeds up the calculation without need of full distance
	 * matrix. The parsing mode
	 * {@link FileParsingParameters#setAlignSeqRes(boolean)} needs to be set to
	 * true for this to work.
	 *
	 * @param chain
	 * @param cutoff
	 * @return
	 */
	public static AtomContactSet getAtomsInContact(Chain chain, double cutoff) {
		return getAtomsInContact(chain, (String[]) null, cutoff);
	}

	/**
	 * Returns the set of intra-chain contacts for the given chain for C-alpha
	 * atoms (including non-standard aminoacids appearing as HETATM groups),
	 * i.e. the contact map. Uses a geometric hashing algorithm that speeds up
	 * the calculation without need of full distance matrix. The parsing mode
	 * {@link FileParsingParameters#setAlignSeqRes(boolean)} needs to be set to
	 * true for this to work.
	 *
	 * @param chain
	 * @param cutoff
	 * @return
	 * @see {@link #getRepresentativeAtomsInContact(Chain, double)}
	 */
	public static AtomContactSet getAtomsCAInContact(Chain chain, double cutoff) {
		Grid grid = new Grid(cutoff);

		Atom[] atoms = getAtomCAArray(chain);

		grid.addAtoms(atoms);

		return grid.getAtomContacts();
	}

	/**
	 * Returns the set of intra-chain contacts for the given chain for C-alpha
	 * or C3' atoms (including non-standard aminoacids appearing as HETATM
	 * groups), i.e. the contact map. Uses a geometric hashing algorithm that
	 * speeds up the calculation without need of full distance matrix.
	 *
	 * @param chain
	 * @param cutoff
	 * @return
	 * @since Biojava 4.1.0
	 */
	public static AtomContactSet getRepresentativeAtomsInContact(Chain chain,
			double cutoff) {
		Grid grid = new Grid(cutoff);

		Atom[] atoms = getRepresentativeAtomArray(chain);

		grid.addAtoms(atoms);

		return grid.getAtomContacts();
	}

	/**
	 * Returns the set of inter-chain contacts between the two given chains for
	 * the given atom names. Uses a geometric hashing algorithm that speeds up
	 * the calculation without need of full distance matrix. The parsing mode
	 * {@link FileParsingParameters#setAlignSeqRes(boolean)} needs to be set to
	 * true for this to work.
	 *
	 * @param chain1
	 * @param chain2
	 * @param atomNames
	 *            the array with atom names to be used. For Calphas use {"CA"},
	 *            if null all non-H atoms will be used. Note HET atoms are
	 *            ignored unless this parameter is null.
	 * @param cutoff
	 * @param hetAtoms
	 *            if true HET atoms are included, if false they are not
	 * @return
	 */
	public static AtomContactSet getAtomsInContact(Chain chain1, Chain chain2,
			String[] atomNames, double cutoff, boolean hetAtoms) {
		Grid grid = new Grid(cutoff);
		Atom[] atoms1 = null;
		Atom[] atoms2 = null;
		if (atomNames == null) {
			atoms1 = getAllNonHAtomArray(chain1, hetAtoms);
			atoms2 = getAllNonHAtomArray(chain2, hetAtoms);
		} else {
			atoms1 = getAtomArray(chain1, atomNames);
			atoms2 = getAtomArray(chain2, atomNames);
		}
		grid.addAtoms(atoms1, atoms2);

		return grid.getAtomContacts();
	}

	/**
	 * Returns the set of inter-chain contacts between the two given chains for
	 * all non-H atoms. Uses a geometric hashing algorithm that speeds up the
	 * calculation without need of full distance matrix. The parsing mode
	 * {@link FileParsingParameters#setAlignSeqRes(boolean)} needs to be set to
	 * true for this to work.
	 *
	 * @param chain1
	 * @param chain2
	 * @param cutoff
	 * @param hetAtoms
	 *            if true HET atoms are included, if false they are not
	 * @return
	 */
	public static AtomContactSet getAtomsInContact(Chain chain1, Chain chain2,
			double cutoff, boolean hetAtoms) {
		return getAtomsInContact(chain1, chain2, null, cutoff, hetAtoms);
	}

	/**
	 * Finds Groups in {@code structure} that contain at least one Atom that is
	 * within {@code radius} Angstroms of {@code centroid}.
	 *
	 * @param structure
	 *            The structure from which to find Groups
	 * @param centroid
	 *            The centroid of the shell
	 * @param excludeResidues
	 *            A list of ResidueNumbers to exclude
	 * @param radius
	 *            The radius from {@code centroid}, in Angstroms
	 * @param includeWater
	 *            Whether to include Groups whose <em>only</em> atoms are water
	 * @param useAverageDistance
	 *            When set to true, distances are the arithmetic mean (1-norm)
	 *            of the distances of atoms that belong to the group and that
	 *            are within the shell; otherwise, distances are the minimum of
	 *            these values
	 * @return A map of Groups within (or partially within) the shell, to their
	 *         distances in Angstroms
	 */
	public static Map<Group, Double> getGroupDistancesWithinShell(
			Structure structure, Atom centroid,
			Set<ResidueNumber> excludeResidues, double radius,
			boolean includeWater, boolean useAverageDistance) {

		// for speed, we avoid calculating square roots
		radius = radius * radius;

		Map<Group, Double> distances = new HashMap<Group, Double>();

		// we only need this if we're averaging distances
		// note that we can't use group.getAtoms().size() because some the
		// group's atoms be outside the shell
		Map<Group, Integer> atomCounts = new HashMap<Group, Integer>();

		for (Chain chain : structure.getChains()) {
			groupLoop: for (Group chainGroup : chain.getAtomGroups()) {

				// exclude water
				if (!includeWater && chainGroup.isWater())
					continue;

				// check blacklist of residue numbers
				for (ResidueNumber rn : excludeResidues) {
					if (rn.equals(chainGroup.getResidueNumber()))
						continue groupLoop;
				}

				for (Atom testAtom : chainGroup.getAtoms()) {

					// use getDistanceFast as we are doing a lot of comparisons
					double dist = Calc.getDistanceFast(centroid, testAtom);

					// if we're the shell
					if (dist <= radius) {
						if (!distances.containsKey(chainGroup))
							distances.put(chainGroup, Double.POSITIVE_INFINITY);
						if (useAverageDistance) {
							// sum the distance; we'll divide by the total
							// number later
							// here, we CANNOT use fastDistance (distance
							// squared) because we want the arithmetic mean
							distances.put(chainGroup, distances.get(chainGroup)
									+ Math.sqrt(dist));
							if (!atomCounts.containsKey(chainGroup))
								atomCounts.put(chainGroup, 0);
							atomCounts.put(chainGroup,
									atomCounts.get(chainGroup) + 1);
						} else {
							// take the minimum distance among all atoms of
							// chainGroup
							// note that we can't break here because we might
							// find a smaller distance
							if (dist < distances.get(chainGroup)) {
								distances.put(chainGroup, dist);
							}
						}
					}

				}
			}
		}

		if (useAverageDistance) {
			for (Map.Entry<Group, Double> entry : distances.entrySet()) {
				int count = atomCounts.get(entry.getKey());
				distances.put(entry.getKey(), entry.getValue() / count);
			}
		} else {
			// in this case we used getDistanceFast
			for (Map.Entry<Group, Double> entry : distances.entrySet()) {
				distances.put(entry.getKey(), Math.sqrt(entry.getValue()));
			}
		}

		return distances;

	}

	public static Set<Group> getGroupsWithinShell(Structure structure,
			Atom atom, Set<ResidueNumber> excludeResidues, double distance,
			boolean includeWater) {

		// square the distance to use as a comparison against getDistanceFast
		// which returns the square of a distance.
		distance = distance * distance;

		Set<Group> returnSet = new LinkedHashSet<Group>();
		for (Chain chain : structure.getChains()) {
			groupLoop: for (Group chainGroup : chain.getAtomGroups()) {
				if (!includeWater && chainGroup.isWater())
					continue;
				for (ResidueNumber rn : excludeResidues) {
					if (rn.equals(chainGroup.getResidueNumber()))
						continue groupLoop;
				}
				for (Atom atomB : chainGroup.getAtoms()) {

					// use getDistanceFast as we are doing a lot of comparisons
					double dist = Calc.getDistanceFast(atom, atomB);
					if (dist <= distance) {
						returnSet.add(chainGroup);
						break;
					}

				}
			}
		}
		return returnSet;
	}

	/**
	 * <p>
	 * Returns a Set of Groups in a structure within the distance specified of a
	 * given group.
	 * </p>
	 * <p>
	 * Updated 18-Sep-2015 sroughley to return a Set so only a unique set of
	 * Groups returned
	 *
	 * @param structure
	 *            The structure to work with
	 * @param group
	 *            The 'query' group
	 * @param distance
	 *            The cutoff distance
	 * @param includeWater
	 *            Should water residues be included in the output?
	 * @return {@link LinkedHashSet} of {@link Group}s within at least one atom
	 *         with {@code distance} of at least one atom in {@code group}
	 */
	public static Set<Group> getGroupsWithinShell(Structure structure,
			Group group, double distance, boolean includeWater) {

		Set<Group> returnList = new LinkedHashSet<Group>();

		Set<ResidueNumber> excludeGroups = new HashSet<ResidueNumber>();
		excludeGroups.add(group.getResidueNumber());
		for (Atom atom : group.getAtoms()) {
			Set<Group> set = getGroupsWithinShell(structure, atom,
					excludeGroups, distance, includeWater);
			returnList.addAll(set);
		}

		return returnList;
	}

	/**
	 * Remove all models from a Structure and keep only the first
	 *
	 * @param s
	 *            original Structure
	 * @return a structure that contains only the first model
	 * @since 3.0.5
	 */
	public static Structure removeModels(Structure s) {
		if (s.nrModels() == 1)
			return s;

		Structure n = new StructureImpl();
		// go through whole substructure and clone ...

		// copy structure data

		n.setPDBCode(s.getPDBCode());
		n.setName(s.getName());

		// TODO: do deep copying of data!
		n.setPDBHeader(s.getPDBHeader());
		n.setDBRefs(s.getDBRefs());

		n.setSites(s.getSites());

		n.setChains(s.getModel(0));

		return n;

	}

	/**
	 * Removes all polymeric and solvent groups from a list of groups
	 *
	 */
	public static List<Group> filterLigands(List<Group> allGroups) {

		List<Group> groups = new ArrayList<Group>();
		for (Group g : allGroups) {

			if ( g.isPolymeric())
				continue;

			if (!g.isWater()) {
				groups.add(g);
			}
		}

		return groups;
	}

	/**
	 * Short version of {@link #getStructure(String, PDBFileParser, AtomCache)}
	 * which creates new parsers when needed
	 *
	 * @param name
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 */
	public static Structure getStructure(String name) throws IOException,
	StructureException {
		return StructureTools.getStructure(name, null, null);
	}

	/**
	 * Flexibly get a structure from an input String. The intent of this method
	 * is to allow any reasonable string which could refer to a structure to be
	 * correctly parsed. The following are currently supported:
	 * <ol>
	 * <li>Filename (if name refers to an existing file)
	 * <li>PDB ID
	 * <li>SCOP domains
	 * <li>PDP domains
	 * <li>Residue ranges
	 * <li>Other formats supported by AtomCache
	 * </ol>
	 *
	 * @param name
	 *            Some reference to the protein structure
	 * @param parser
	 *            A clean PDBFileParser to use if it is a file. If null, a
	 *            PDBFileParser will be instantiated if needed.
	 * @param cache
	 *            An AtomCache to use if the structure can be fetched from the
	 *            PDB. If null, a AtomCache will be instantiated if needed.
	 * @return A Structure object
	 * @throws IOException
	 *             if name is an existing file, but doesn't parse correctly
	 * @throws StructureException
	 *             if the format is unknown, or if AtomCache throws an
	 *             exception.
	 */
	public static Structure getStructure(String name, PDBFileParser parser,
			AtomCache cache) throws IOException, StructureException {
		File f = new File(FileDownloadUtils.expandUserHome(name));
		if (f.exists()) {
			if (parser == null) {
				parser = new PDBFileParser();
			}
			InputStream inStream = new FileInputStream(f);
			return parser.parsePDBFile(inStream);
		} else {
			if (cache == null) {
				cache = new AtomCache();
			}
			return cache.getStructure(name);
		}
	}

	/**
	 * @deprecated  use {@link Chain#isProtein()} instead.
	 */
	@Deprecated
	public static boolean isProtein(Chain c) {

		return c.isProtein();
	}

	/**
	 * @deprecated use {@link Chain#isNucleicAcid()} instead.
 	 */
	@Deprecated
	public static boolean isNucleicAcid(Chain c) {
		return c.isNucleicAcid();
	}

	/**
	 * @deprecated use {@link Chain#getPredominantGroupType()} instead.
	 */
	@Deprecated
	public static GroupType getPredominantGroupType(Chain c) {
		return c.getPredominantGroupType();
	}

	/**
	 * @deprecated use {@link Chain#isWaterOnly()} instead.
	 */
	@Deprecated
	public static boolean isChainWaterOnly(Chain c) {
		return c.isWaterOnly();
	}

	/**
     * @deprecated  use {@link Chain#isPureNonPolymer()} instead.
	 */
	@Deprecated
	public static boolean isChainPureNonPolymer(Chain c) {

		return c.isPureNonPolymer();
	}

	/**
	 * Cleans up the structure's alternate location (altloc) groups. All alternate location groups should have all atoms (except
	 * in the case of microheterogenity) or when a deuterium exists.
	 * Ensure that all the alt loc groups have all the atoms in the main group.
	 * @param structure The Structure to be cleaned up
	 */
	public static void cleanUpAltLocs(Structure structure) {
		for (int i =0; i< structure.nrModels() ; i++){
			for (Chain chain : structure.getModel(i)) {
				for (Group group : chain.getAtomGroups()) {
					for (Group altLocGroup : group.getAltLocs()) { 
						for ( Atom groupAtom : group.getAtoms()) {
							// If this alt loc doesn't have this atom
							if (! altLocGroup.hasAtom(groupAtom.getName())) {
								// Fix for microheterogenity
								if (altLocGroup.getPDBName().equals(group.getPDBName())) {
									// If it's a Hydrogen then we check for it's Deuterated brother
									if(!hasDeuteratedEquiv(groupAtom, altLocGroup)){
										altLocGroup.addAtom(groupAtom);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Check to see if an Deuterated atom has a non deuterated brother in the group.
	 * @param atom the input atom that is putatively deuterium
	 * @param currentGroup the group the atom is in
	 * @return true if the atom is deuterated and it's hydrogen equive exists.
	 */
	public static boolean hasNonDeuteratedEquiv(Atom atom, Group currentGroup) {
		if(atom.getElement()==Element.D && currentGroup.hasAtom(replaceFirstChar(atom.getName(),'D', 'H'))) {
			// If it's deuterated and has a non-deuterated brother
			return true;
		}
		return false;
	}
	
	/**
	 * Check to see if a Hydrogen has a  Deuterated brother in the group.
	 * @param atom the input atom that is putatively hydorgen
	 * @param currentGroup the group the atom is in
	 * @return true if the atom is hydrogen and it's Deuterium equiv exists.
	 */
	public static boolean hasDeuteratedEquiv(Atom atom, Group currentGroup) {
		if(atom.getElement()==Element.H && currentGroup.hasAtom(replaceFirstChar(atom.getName(),'H', 'D'))) {
			// If it's hydrogen and has a deuterated brother
			return true;
		}
		return false;
	}

	private static String replaceFirstChar(String name, char c, char d) {
		if(name.charAt(0)==c){
			return name.replaceFirst(String.valueOf(c), String.valueOf(d));
		}
		return name;
	}

}
