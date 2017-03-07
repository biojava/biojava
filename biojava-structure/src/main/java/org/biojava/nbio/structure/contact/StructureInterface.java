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
package org.biojava.nbio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.nbio.structure.asa.GroupAsa;
import org.biojava.nbio.structure.io.FileConvert;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.MMCIFFileTools;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.model.AtomSite;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.xtal.CrystalTransform;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * An interface between 2 molecules (2 sets of atoms).
 *
 * @author duarte_j
 *
 */
public class StructureInterface implements Serializable, Comparable<StructureInterface> {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(StructureInterface.class);

	/**
	 * Interfaces with larger inverse self contact overlap score will be considered isologous
	 */
	private static final double SELF_SCORE_FOR_ISOLOGOUS = 0.3;

	private int id;
	private double totalArea;
	private AtomContactSet contacts;
	private GroupContactSet groupContacts;

	private Pair<Atom[]> molecules;

	/**
	 * The identifier for each of the atom arrays (usually a chain identifier, i.e. a single capital letter)
	 * Serves to identify the molecules within the Asymmetric Unit of the crystal
	 */
	private Pair<String> moleculeIds;

	/**
	 * The transformations (crystal operators) applied to each molecule (if applicable)
	 */
	private Pair<CrystalTransform> transforms;

	private Map<ResidueNumber, GroupAsa> groupAsas1;
	private Map<ResidueNumber, GroupAsa> groupAsas2;

	private StructureInterfaceCluster cluster;

	/**
	 * Constructs a StructureInterface
	 * @param firstMolecule the atoms of the first molecule
	 * @param secondMolecule the atoms of the second molecule
	 * @param firstMoleculeId an identifier that identifies the first molecule within the Asymmetric Unit
	 * @param secondMoleculeId an identifier that identifies the second molecule within the Asymmetric Unit
	 * @param contacts the contacts between the 2 molecules
	 * @param firstTransf the transformation (crystal operator) applied to first molecule
	 * @param secondTransf the transformation (crystal operator) applied to second molecule
	 */
	public StructureInterface(
			Atom[] firstMolecule, Atom[] secondMolecule,
			String firstMoleculeId, String secondMoleculeId,
			AtomContactSet contacts,
			CrystalTransform firstTransf, CrystalTransform secondTransf) {

		this.molecules = new Pair<Atom[]>(firstMolecule, secondMolecule);
		this.moleculeIds = new Pair<String>(firstMoleculeId,secondMoleculeId);
		this.contacts = contacts;
		this.transforms = new Pair<CrystalTransform>(firstTransf, secondTransf);
	}

	/**
	 * Constructs an empty StructureInterface
	 */
	public StructureInterface() {
		this.groupAsas1 = new TreeMap<ResidueNumber, GroupAsa>();
		this.groupAsas2 = new TreeMap<ResidueNumber, GroupAsa>();
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	/**
	 * Returns a pair of identifiers for each of the 2 member molecules that
	 * identify them uniquely in the crystal:
	 *   &lt;molecule id (asym unit id)&gt;+&lt;operator id&gt;+&lt;crystal translation&gt;
	 * @return
	 */
	public Pair<String> getCrystalIds() {
		return new Pair<String>(
			moleculeIds.getFirst()+transforms.getFirst().getTransformId()+transforms.getFirst().getCrystalTranslation(),
			moleculeIds.getSecond()+transforms.getSecond().getTransformId()+transforms.getSecond().getCrystalTranslation());
	}

	/**
	 * Returns the total area buried upon formation of this interface,
	 * defined as: 1/2[ (ASA1u-ASA1c) + (ASA2u-ASA2u) ] , with:
	 *  <p>ASAxu = ASA of first/second unbound chain</p>
	 *  <p>ASAxc = ASA of first/second complexed chain</p>
	 * In the area calculation HETATOM groups not part of the main protein/nucleotide chain
	 * are not included.
	 * @return
	 */
	public double getTotalArea() {
		return totalArea;
	}

	public void setTotalArea(double totalArea) {
		this.totalArea = totalArea;
	}

	public AtomContactSet getContacts() {
		return contacts;
	}

	public void setContacts(AtomContactSet contacts) {
		this.contacts = contacts;
	}

	public Pair<Atom[]> getMolecules() {
		return molecules;
	}

	public void setMolecules(Pair<Atom[]> molecules) {
		this.molecules = molecules;
	}

	/**
	 * Return the pair of identifiers identifying each of the 2 molecules of this interface
	 * in the asymmetry unit (usually the chain identifier if this interface is between 2 chains)
	 * @return
	 */
	public Pair<String> getMoleculeIds() {
		return moleculeIds;
	}

	public void setMoleculeIds(Pair<String> moleculeIds) {
		this.moleculeIds = moleculeIds;
	}

	/**
	 * Return the 2 crystal transform operations performed on each of the
	 * molecules of this interface.
	 * @return
	 */
	public Pair<CrystalTransform> getTransforms() {
		return transforms;
	}

	public void setTransforms(Pair<CrystalTransform> transforms) {
		this.transforms = transforms;
	}

	protected void setAsas(double[] asas1, double[] asas2, int nSpherePoints, int nThreads, int cofactorSizeToUse) {

		Atom[] atoms = getAtomsForAsa(cofactorSizeToUse);
		AsaCalculator asaCalc = new AsaCalculator(atoms,
				AsaCalculator.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);

		double[] complexAsas = asaCalc.calculateAsas();

		if (complexAsas.length!=asas1.length+asas2.length)
			throw new IllegalArgumentException("The size of ASAs of complex doesn't match that of ASAs 1 + ASAs 2");


		groupAsas1 = new TreeMap<ResidueNumber, GroupAsa>();
		groupAsas2 = new TreeMap<ResidueNumber, GroupAsa>();

		this.totalArea = 0;

		for (int i=0;i<asas1.length;i++) {
			Group g = atoms[i].getGroup();

			if (!g.getType().equals(GroupType.HETATM) ||
				isInChain(g)) {
				// interface area should be only for protein/nucleotide but not hetatoms that are not part of the chain
				this.totalArea += (asas1[i] - complexAsas[i]);
			}

			if (!groupAsas1.containsKey(g.getResidueNumber())) {
				GroupAsa groupAsa = new GroupAsa(g);
				groupAsa.addAtomAsaU(asas1[i]);
				groupAsa.addAtomAsaC(complexAsas[i]);
				groupAsas1.put(g.getResidueNumber(), groupAsa);
			} else {
				GroupAsa groupAsa = groupAsas1.get(g.getResidueNumber());
				groupAsa.addAtomAsaU(asas1[i]);
				groupAsa.addAtomAsaC(complexAsas[i]);
			}
		}

		for (int i=0;i<asas2.length;i++) {
			Group g = atoms[i+asas1.length].getGroup();

			if (!g.getType().equals(GroupType.HETATM) ||
				isInChain(g)) {
				// interface area should be only for protein/nucleotide but not hetatoms that are not part of the chain
				this.totalArea += (asas2[i] - complexAsas[i+asas1.length]);
			}

			if (!groupAsas2.containsKey(g.getResidueNumber())) {
				GroupAsa groupAsa = new GroupAsa(g);
				groupAsa.addAtomAsaU(asas2[i]);
				groupAsa.addAtomAsaC(complexAsas[i+asas1.length]);
				groupAsas2.put(g.getResidueNumber(), groupAsa);
			} else {
				GroupAsa groupAsa = groupAsas2.get(g.getResidueNumber());
				groupAsa.addAtomAsaU(asas2[i]);
				groupAsa.addAtomAsaC(complexAsas[i+asas1.length]);
			}
		}

		// our interface area definition: average of bsa of both molecules
		this.totalArea = this.totalArea/2.0;

	}

	protected Atom[] getFirstAtomsForAsa(int cofactorSizeToUse) {

		return getAllNonHAtomArray(molecules.getFirst(), cofactorSizeToUse);
	}

	protected Atom[] getSecondAtomsForAsa(int cofactorSizeToUse) {

		return getAllNonHAtomArray(molecules.getSecond(), cofactorSizeToUse);
	}

	protected Atom[] getAtomsForAsa(int cofactorSizeToUse) {
		Atom[] atoms1 = getFirstAtomsForAsa(cofactorSizeToUse);
		Atom[] atoms2 = getSecondAtomsForAsa(cofactorSizeToUse);

		Atom[] atoms = new Atom[atoms1.length+atoms2.length];
		for (int i=0;i<atoms1.length;i++) {
			atoms[i] = atoms1[i];
		}
		for (int i=0;i<atoms2.length;i++) {
			atoms[i+atoms1.length] = atoms2[i];
		}

		return atoms;
	}

	/**
	 * Returns and array of all non-Hydrogen atoms in the given molecule, including all
	 * main chain HETATOM groups. Non main-chain HETATOM groups with fewer than minSizeHetAtomToInclude
	 * non-Hydrogen atoms are not included.
	 * @param m
	 * @param minSizeHetAtomToInclude HETATOM groups (non main-chain) with fewer number of
	 * non-Hydrogen atoms are not included
	 * @return
	 */
	private static final Atom[] getAllNonHAtomArray(Atom[] m, int minSizeHetAtomToInclude) {
		List<Atom> atoms = new ArrayList<Atom>();

		for (Atom a:m){

			if (a.getElement()==Element.H) continue;

			Group g = a.getGroup();
			if (g.getType().equals(GroupType.HETATM) &&
				!isInChain(g) &&
				getSizeNoH(g)<minSizeHetAtomToInclude) {
				continue;
			}

			atoms.add(a);

		}
		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * Calculates the number of non-Hydrogen atoms in the given group
	 * @param g
	 * @return
	 */
	private static int getSizeNoH(Group g) {
		int size = 0;
		for (Atom a:g.getAtoms()) {
			if (a.getElement()!=Element.H)
				size++;
		}
		return size;
	}

	/**
	 * Returns true if the given group is part of the main chain, i.e. if it is
	 * a peptide-linked group or a nucleotide
	 * @param g
	 * @return
	 */
	private static boolean isInChain(Group g) {
		ChemComp chemComp = g.getChemComp();

		if (chemComp==null) {
			logger.warn("Warning: can't determine PolymerType for group "+g.getResidueNumber()+" ("+g.getPDBName()+"). Will consider it as non-nucleotide/non-protein type.");
			return false;
		}

		PolymerType polyType = chemComp.getPolymerType();
		for (PolymerType protOnlyType: PolymerType.PROTEIN_ONLY) {
			if (polyType==protOnlyType) return true;
		}
		for (PolymerType protOnlyType: PolymerType.POLYNUCLEOTIDE_ONLY) {
			if (polyType==protOnlyType) return true;
		}

		return false;
	}

	/**
	 * Tells whether the interface corresponds to one mediated by crystallographic symmetry,
	 * i.e. it is between symmetry-related molecules (with same chain identifier)
	 * @return
	 */
	public boolean isSymRelated() {
		return moleculeIds.getFirst().equals(moleculeIds.getSecond());
	}

	/**
	 * Returns true if the transformation applied to the second molecule of this interface
	 * has an infinite character (pure translation or screw rotation)
	 * and both molecules of the interface have the same asymmetric unit identifier (chain id): in such cases the
	 * interface would lead to infinite fiber-like (linear or helical) assemblies
	 * @return
	 */
	public boolean isInfinite() {
		return ((isSymRelated() && transforms.getSecond().getTransformType().isInfinite()));
	}

	/**
	 * Returns true if the 2 molecules of this interface are the same entity (i.e. homomeric interface), false
	 * otherwise (i.e. heteromeric interface)
	 * @return true if homomeric or if either of the entities is unknonw (null Compounds), false otherwise
	 */
	public boolean isHomomeric() {
		EntityInfo first = getParentChains().getFirst().getEntityInfo();
		EntityInfo second = getParentChains().getSecond().getEntityInfo();
		if (first==null || second==null) {
			logger.warn("Some compound of interface {} is null, can't determine whether it is homo/heteromeric. Consider it homomeric", getId());
			return true;
		}
		return
			first.getRepresentative().getId().equals(second.getRepresentative().getId());
	}

	/**
	 * Gets a map of ResidueNumbers to GroupAsas for all groups of first chain.
	 * @return
	 */
	public Map<ResidueNumber, GroupAsa> getFirstGroupAsas() {
		return groupAsas1;
	}

	/**
	 * Gets the GroupAsa for the corresponding residue number of first chain
	 * @param resNum
	 * @return
	 */
	public GroupAsa getFirstGroupAsa(ResidueNumber resNum) {
		return groupAsas1.get(resNum);
	}

	public void setFirstGroupAsa(GroupAsa groupAsa) {
		groupAsas1.put(groupAsa.getGroup().getResidueNumber(), groupAsa);
	}

	/**
	 * Gets a map of ResidueNumbers to GroupAsas for all groups of second chain.
	 * @return
	 */
	public Map<ResidueNumber, GroupAsa> getSecondGroupAsas() {
		return groupAsas2;
	}

	public void setSecondGroupAsa(GroupAsa groupAsa) {
		groupAsas2.put(groupAsa.getGroup().getResidueNumber(), groupAsa);
	}

	/**
	 * Gets the GroupAsa for the corresponding residue number of second chain
	 * @param resNum
	 * @return
	 */
	public GroupAsa getSecondGroupAsa(ResidueNumber resNum) {
		return groupAsas2.get(resNum);
	}

	/**
	 * Returns the residues belonging to the interface core, defined as those residues at
	 * the interface (BSA>0) and for which the BSA/ASA ratio is above the given bsaToAsaCutoff
	 * @param bsaToAsaCutoff
	 * @param minAsaForSurface the minimum ASA to consider a residue on the surface
	 * @return
	 */
	public Pair<List<Group>> getCoreResidues(double bsaToAsaCutoff, double minAsaForSurface) {

		List<Group> core1 = new ArrayList<Group>();
		List<Group> core2 = new ArrayList<Group>();

		for (GroupAsa groupAsa:groupAsas1.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface && groupAsa.getBsa()>0) {
				if (groupAsa.getBsaToAsaRatio()<bsaToAsaCutoff) {
					//rim1.add(groupAsa.getGroup());
				} else {
					core1.add(groupAsa.getGroup());
				}
			}
		}
		for (GroupAsa groupAsa:groupAsas2.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface && groupAsa.getBsa()>0) {
				if (groupAsa.getBsaToAsaRatio()<bsaToAsaCutoff) {
					//rim2.add(groupAsa.getGroup());
				} else {
					core2.add(groupAsa.getGroup());
				}
			}
		}

		return new Pair<List<Group>>(core1, core2);
	}

	/**
	 * Returns the residues belonging to the interface rim, defined as those residues at
	 * the interface (BSA>0) and for which the BSA/ASA ratio is below the given bsaToAsaCutoff
	 * @param bsaToAsaCutoff
	 * @param minAsaForSurface the minimum ASA to consider a residue on the surface
	 * @return
	 */
	public Pair<List<Group>> getRimResidues(double bsaToAsaCutoff, double minAsaForSurface) {

		List<Group> rim1 = new ArrayList<Group>();
		List<Group> rim2 = new ArrayList<Group>();

		for (GroupAsa groupAsa:groupAsas1.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface && groupAsa.getBsa()>0) {
				if (groupAsa.getBsaToAsaRatio()<bsaToAsaCutoff) {
					rim1.add(groupAsa.getGroup());
				} else {
					//core1.add(groupAsa.getGroup());
				}
			}
		}
		for (GroupAsa groupAsa:groupAsas2.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface && groupAsa.getBsa()>0) {
				if (groupAsa.getBsaToAsaRatio()<bsaToAsaCutoff) {
					rim2.add(groupAsa.getGroup());
				} else {
					//core2.add(groupAsa.getGroup());
				}
			}
		}

		return new Pair<List<Group>>(rim1, rim2);
	}

	/**
	 * Returns the residues belonging to the interface, i.e. the residues
	 * at the surface with BSA>0
	 * @param minAsaForSurface the minimum ASA to consider a residue on the surface
	 * @return
	 */
	public Pair<List<Group>> getInterfacingResidues(double minAsaForSurface) {

		List<Group> interf1 = new ArrayList<Group>();
		List<Group> interf2 = new ArrayList<Group>();

		for (GroupAsa groupAsa:groupAsas1.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface && groupAsa.getBsa()>0) {
				interf1.add(groupAsa.getGroup());
			}
		}
		for (GroupAsa groupAsa:groupAsas2.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface && groupAsa.getBsa()>0) {
				interf2.add(groupAsa.getGroup());
			}
		}

		return new Pair<List<Group>>(interf1, interf2);
	}

	/**
	 * Returns the residues belonging to the surface
	 * @param minAsaForSurface the minimum ASA to consider a residue on the surface
	 * @return
	 */
	public Pair<List<Group>> getSurfaceResidues(double minAsaForSurface) {
		List<Group> surf1 = new ArrayList<Group>();
		List<Group> surf2 = new ArrayList<Group>();

		for (GroupAsa groupAsa:groupAsas1.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface) {
				surf1.add(groupAsa.getGroup());
			}
		}
		for (GroupAsa groupAsa:groupAsas2.values()) {

			if (groupAsa.getAsaU()>minAsaForSurface) {
				surf2.add(groupAsa.getGroup());
			}
		}

		return new Pair<List<Group>>(surf1, surf2);
	}

	public StructureInterfaceCluster getCluster() {
		return cluster;
	}

	public void setCluster(StructureInterfaceCluster cluster) {
		this.cluster = cluster;
	}

	/**
	 * Calculates the contact overlap score between this StructureInterface and
	 * the given one.
	 * The two sides of the given StructureInterface need to match this StructureInterface
	 * in the sense that they must come from the same Compound (Entity), i.e.
	 * their residue numbers need to align with 100% identity, except for unobserved
	 * density residues. The SEQRES indices obtained through {@link EntityInfo#getAlignedResIndex(Group, Chain)} are
	 * used to match residues, thus if no SEQRES is present or if {@link FileParsingParameters#setAlignSeqRes(boolean)}
	 * is not used, this calculation is not guaranteed to work properly.
	 * @param other
	 * @param invert if false the comparison will be done first-to-first and second-to-second,
	 * if true the match will be first-to-second and second-to-first
	 * @return the contact overlap score, range [0.0,1.0]
	 */
	public double getContactOverlapScore(StructureInterface other, boolean invert) {

		Structure thisStruct = getParentStructure();
		Structure otherStruct = other.getParentStructure();

		if (thisStruct!=otherStruct) {
			// in the current implementation, comparison between different structure doesn't make much sense
			// and won't even work since the compounds of both will never match. We warn because it
			// really is not what this is intended for at the moment
			logger.warn("Comparing interfaces from different structures, contact overlap score will be 0");
			return 0;
		}

		Pair<Chain> thisChains = getParentChains();
		Pair<Chain> otherChains = other.getParentChains();

		if (thisChains.getFirst().getEntityInfo() == null || thisChains.getSecond().getEntityInfo() == null ||
			otherChains.getFirst().getEntityInfo() == null || otherChains.getSecond().getEntityInfo() == null ) {
			// this happens in cases like 2uub
			logger.warn("Found chains with null compounds while comparing interfaces {} and {}. Contact overlap score for them will be 0.",
					this.getId(), other.getId());
			return 0;
		}

		Pair<EntityInfo> thisCompounds = new Pair<EntityInfo>(thisChains.getFirst().getEntityInfo(), thisChains.getSecond().getEntityInfo());
		Pair<EntityInfo> otherCompounds = new Pair<EntityInfo>(otherChains.getFirst().getEntityInfo(), otherChains.getSecond().getEntityInfo());

		if ( (  (thisCompounds.getFirst() == otherCompounds.getFirst()) &&
				(thisCompounds.getSecond() == otherCompounds.getSecond())   )  ||
			 (  (thisCompounds.getFirst() == otherCompounds.getSecond()) &&
				(thisCompounds.getSecond() == otherCompounds.getFirst())   )	) {

			int common = 0;
			GroupContactSet thisContacts = getGroupContacts();
			GroupContactSet otherContacts = other.getGroupContacts();

			for (GroupContact thisContact:thisContacts) {

				ResidueIdentifier first = null;
				ResidueIdentifier second = null;

				if (!invert) {
					first = new ResidueIdentifier(thisContact.getPair().getFirst());

					second = new ResidueIdentifier(thisContact.getPair().getSecond());
				} else {
					first = new ResidueIdentifier(thisContact.getPair().getSecond());

					second = new ResidueIdentifier(thisContact.getPair().getFirst());
				}

				if (otherContacts.hasContact(first,second)) {
					common++;
				}
			}
			return (2.0*common)/(thisContacts.size()+otherContacts.size());
		} else {
			logger.debug("Chain pairs {},{} and {},{} belong to different compound pairs, contact overlap score will be 0 ",
					thisChains.getFirst().getId(),thisChains.getSecond().getId(),
					otherChains.getFirst().getId(),otherChains.getSecond().getId());
			return 0.0;
		}
	}

	public GroupContactSet getGroupContacts() {
		if (groupContacts==null) {
			this.groupContacts  = new GroupContactSet(contacts);
		}
		return this.groupContacts;
	}

	/**
	 * Tell whether the interface is isologous, i.e. it is formed
	 * by the same patches of same Compound on both sides.
	 *
	 * @return true if isologous, false if heterologous
	 */
	public boolean isIsologous() {
		double scoreInverse = this.getContactOverlapScore(this, true);
		logger.debug("Interface {} contact overlap score with itself inverted: {}",
				getId(), scoreInverse);
		return (scoreInverse>SELF_SCORE_FOR_ISOLOGOUS);
	}

	/**
	 * Finds the parent chains by looking up the references of first atom of each side of this interface
	 * @return
	 */
	public Pair<Chain> getParentChains() {
		Atom[] firstMol = this.molecules.getFirst();
		Atom[] secondMol = this.molecules.getSecond();
		if (firstMol.length==0 || secondMol.length==0) {
			logger.warn("No atoms found in first or second molecule, can't get parent Chains");
			return null;
		}

		return new Pair<Chain>(firstMol[0].getGroup().getChain(), secondMol[0].getGroup().getChain());
	}

	/**
	 * Finds the parent compounds by looking up the references of first atom of each side of this interface
	 * @return
	 */
	public Pair<EntityInfo> getParentCompounds() {
		Pair<Chain> chains = getParentChains();
		if (chains == null) {
			logger.warn("Could not find parents chains, compounds will be null");
			return null;
		}
		return new Pair<EntityInfo>(chains.getFirst().getEntityInfo(), chains.getSecond().getEntityInfo());
	}

	private Structure getParentStructure() {
		Atom[] firstMol = this.molecules.getFirst();
		if (firstMol.length==0) {
			logger.warn("No atoms found in first molecule, can't get parent Structure");
			return null;
		}
		return firstMol[0].getGroup().getChain().getStructure();
	}

	/**
	 * Return a String representing the 2 molecules of this interface in PDB format.
	 * If the molecule ids (i.e. chain ids) are the same for both molecules, then the second
	 * one will be replaced by the next letter in alphabet (or A for Z)
	 * @return
	 */
	public String toPDB() {

		String molecId1 = getMoleculeIds().getFirst();
		String molecId2 = getMoleculeIds().getSecond();

		if (molecId2.equals(molecId1)) {
			// if both chains are named equally we want to still named them differently in the output pdb file
			// so that molecular viewers can handle properly the 2 chains as separate entities
			char letter = molecId1.charAt(0);
			if (letter!='Z' && letter!='z') {
				molecId2 = Character.toString((char)(letter+1)); // i.e. next letter in alphabet
			} else {
				molecId2 = Character.toString((char)(letter-25)); //i.e. 'A' or 'a'
			}
		}

		StringBuilder sb = new StringBuilder();
		for (Atom atom:this.molecules.getFirst()) {
			sb.append(FileConvert.toPDB(atom, molecId1));
		}
		sb.append("TER");
		sb.append(System.getProperty("line.separator"));
		for (Atom atom:this.molecules.getSecond()) {
			sb.append(FileConvert.toPDB(atom,molecId2));
		}
		sb.append("TER");
		sb.append(System.getProperty("line.separator"));
		sb.append("END");
		sb.append(System.getProperty("line.separator"));
		return sb.toString();
	}

	/**
	 * Return a String representing the 2 molecules of this interface in mmCIF format.
	 * If the molecule ids (i.e. chain ids) are the same for both molecules, then the second
	 * one will be written as chainId_operatorId (with operatorId taken from {@link #getTransforms()}
	 * @return
	 */
	public String toMMCIF() {
		StringBuilder sb = new StringBuilder();

		String molecId1 = getMoleculeIds().getFirst();
		String molecId2 = getMoleculeIds().getSecond();

		if (isSymRelated()) {
			// if both chains are named equally we want to still named them differently in the output mmcif file
			// so that molecular viewers can handle properly the 2 chains as separate entities
			molecId2 = molecId2 + "_" +getTransforms().getSecond().getTransformId();
		}

		sb.append(SimpleMMcifParser.MMCIF_TOP_HEADER).append("BioJava_interface_").append(getId()).append(System.getProperty("line.separator"));

		sb.append(FileConvert.getAtomSiteHeader());

		// we reassign atom ids if sym related (otherwise atom ids would be duplicated and some molecular viewers can't cope with that)
		int atomId = 1;
		List<AtomSite> atomSites = new ArrayList<>();
		for (Atom atom:this.molecules.getFirst()) {
			if (isSymRelated()) {
				atomSites.add(MMCIFFileTools.convertAtomToAtomSite(atom, 1, molecId1, molecId1, atomId));
			} else {
				atomSites.add(MMCIFFileTools.convertAtomToAtomSite(atom, 1, molecId1, molecId1));
			}
			atomId++;
		}
		for (Atom atom:this.molecules.getSecond()) {
			if (isSymRelated()) {
				atomSites.add(MMCIFFileTools.convertAtomToAtomSite(atom, 1, molecId2, molecId2, atomId));
			} else {
				atomSites.add(MMCIFFileTools.convertAtomToAtomSite(atom, 1, molecId2, molecId2));
			}
			atomId++;
		}

		sb.append(MMCIFFileTools.toMMCIF(atomSites,AtomSite.class));

		return sb.toString();
	}

	@Override
	public int compareTo(StructureInterface o) {
		// this will sort descending on interface areas
		return (Double.compare(o.totalArea,this.totalArea));
	}

	@Override
	public String toString() {
		return String.format("StructureInterface %d (%s, %.0f A, <%s; %s>)", id, moleculeIds,totalArea,transforms.getFirst().toXYZString(),transforms.getSecond().toXYZString());
	}

}
