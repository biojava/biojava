package org.biojava.bio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.asa.AsaCalculator;
import org.biojava.bio.structure.asa.GroupAsa;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava.bio.structure.xtal.CrystalTransform;

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
	
	private int id;
	private double totalArea;
	private AtomContactSet contacts;
	
	private Pair<Atom[]> molecules;

	/**
	 * The identifier for each of the atom arrays (usually a chain identifier, i.e. a single capital letter)
	 * Servers to identify the molecules within the Asymmetric Unit of the crystal
	 */
	private Pair<String> moleculeIds; 
	
	/**
	 * The transformations (crystal operators) applied to each molecule (if applicable)
	 */
	private Pair<CrystalTransform> transforms;  
		
	private Map<ResidueNumber, GroupAsa> groupAsas1;
	private Map<ResidueNumber, GroupAsa> groupAsas2;
	
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

	public Pair<String> getMoleculeIds() {
		return moleculeIds;
	}
	
	/**
	 * Return the 2 crystal transform operations performed on each of the 
	 * chains of this interface.
	 * @return
	 */
	public Pair<CrystalTransform> getTransforms() {
		return transforms;
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
	 * and both chains of the interface have the same PDB chain code: in such cases the 
	 * interface would lead to infinite fiber-like (linear or helical) assemblies
	 * @return
	 */
	public boolean isInfinite() {
		return ((isSymRelated() && transforms.getSecond().getTransformType().isInfinite()));
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
	
	/**
	 * Gets a map of ResidueNumbers to GroupAsas for all groups of second chain.
	 * @return
	 */
	public Map<ResidueNumber, GroupAsa> getSecondGroupAsas() {
		return groupAsas2;
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
