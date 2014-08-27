package org.biojava.bio.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.structure.asa.AsaCalculator;
import org.biojava.bio.structure.asa.GroupAsa;
import org.biojava.bio.structure.contact.AtomContactSet;
import org.biojava.bio.structure.contact.Pair;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava.bio.structure.xtal.CrystalTransform;


/**
 * An interface between 2 polymer (protein/nucleotide) chains.
 * 
 * @author duarte_j
 *
 */
public class ChainInterface implements Serializable, Comparable<ChainInterface> { 

	private static final long serialVersionUID = 1L;
	
	private int id;
	private double totalArea;
	private AtomContactSet contacts;
	
	private Pair<Chain> chains;
	
	private Pair<CrystalTransform> transforms;  // the transformations (crystal operators) applied to each molecule
		
	private Map<ResidueNumber, GroupAsa> groupAsas1;
	private Map<ResidueNumber, GroupAsa> groupAsas2;
	
	
	public ChainInterface(Chain firstMolecule, Chain secondMolecule, AtomContactSet contacts, CrystalTransform firstTransf, CrystalTransform secondTransf) {
		this.chains = new Pair<Chain>(firstMolecule, secondMolecule);
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
	 * Returns a pair of identifiers for each of the 2 member chains that
	 * identify them uniquely in the crystal:
	 *   &lt;chain id&gt;+&lt;operator id&gt;+&lt;crystal translation&gt;
	 * @return
	 */
	public Pair<String> getCrystalIds() {
		return new Pair<String>(
			chains.getFirst().getChainID()+transforms.getFirst().getTransformId()+transforms.getFirst().getCrystalTranslation(),
			chains.getSecond().getChainID()+transforms.getSecond().getTransformId()+transforms.getSecond().getCrystalTranslation());
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

	public Pair<Chain> getChains() {
		return chains;
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
		
		Atom[] atoms = getAtoms(cofactorSizeToUse);
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
	
	protected Atom[] getFirstAtoms(int cofactorSizeToUse) {
		
		Atom[] atoms1 = getAllNonHAtomArray(chains.getFirst(), cofactorSizeToUse);		
		
		return atoms1;		
	}
	
	protected Atom[] getSecondAtoms(int cofactorSizeToUse) {
		
		Atom[] atoms2 = getAllNonHAtomArray(chains.getSecond(), cofactorSizeToUse);
		
		return atoms2;
	}
	
	protected Atom[] getAtoms(int cofactorSizeToUse) {
		Atom[] atoms1 = getFirstAtoms(cofactorSizeToUse);
		Atom[] atoms2 = getSecondAtoms(cofactorSizeToUse);
		
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
	 * Returns and array of all non-Hydrogen atoms in the given Chain, including all
	 * main chain HETATOM groups. Non main-chain HETATOM groups with fewer than minSizeHetAtomToInclude
	 * non-Hydrogen atoms are not included.
	 * @param c
	 * @param minSizeHetAtomToInclude HETATOMs (non main-chain) with fewer number of 
	 * non-Hydrogen atoms are not included
	 * @return
	 */
	private static final Atom[] getAllNonHAtomArray(Chain c, int minSizeHetAtomToInclude) {
		List<Atom> atoms = new ArrayList<Atom>();
		
		for (Group g:c.getAtomGroups()){
			
			if (g.getType().equals(GroupType.HETATM) && 
				!isInChain(g) &&
				getSizeNoH(g)<minSizeHetAtomToInclude) {
				continue;	
			}
			
			for (Atom a:g.getAtoms()) {

				if (a.getElement()==Element.H) continue;

				atoms.add(a);
			}
		}
		return (Atom[]) atoms.toArray(new Atom[atoms.size()]);			
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
			// TODO is there a better solution? we should at least do this through a logger once we have slf4j in BJ
			System.err.println("Warning: can't determine PolymerType for group "+g.getResidueNumber()+" ("+g.getPDBName()+"). Will consider it as non-nucleotide/non-protein type.");
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
	 * i.e. it is between symmetry-related chains (with same chain identifier)
	 * @return
	 */
	public boolean isSymRelated() {
		return chains.getFirst().getChainID().equals(chains.getSecond().getChainID());
	}

	/**
	 * Returns true if the transformation applied to the second chain of this interface
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
	public int compareTo(ChainInterface o) {
		// this will sort descending on interface areas
		return (Double.compare(o.totalArea,this.totalArea));
	}
	
}
