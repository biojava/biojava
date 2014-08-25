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
import org.biojava.bio.structure.xtal.CrystalTransform;


/**
 * An interface between 2 protein chains.
 * 
 * @author duarte_j
 *
 */
public class ChainInterface implements Serializable, Comparable<ChainInterface> { 

	private static final long serialVersionUID = 1L;
	
	private int id;
	private double interfaceArea;
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
	 *   <chain id>+<operator id>+<crystal translation>
	 * @return
	 */
	public Pair<String> getCrystalIds() {
		return new Pair<String>(
			chains.getFirst().getChainID()+transforms.getFirst().getTransformId()+transforms.getFirst().getCrystalTranslation(),
			chains.getSecond().getChainID()+transforms.getSecond().getTransformId()+transforms.getSecond().getCrystalTranslation());
	}

	public double getInterfaceArea() {
		return interfaceArea;
	}

	public void setInterfaceArea(double interfaceArea) {
		this.interfaceArea = interfaceArea;
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

	public Pair<CrystalTransform> getTransforms() {
		return transforms;
	}

	protected void setAsas(double[] asas1, double[] asas2, int nSpherePoints, int nThreads, boolean hetAtoms, int cofactorSizeToUse) {
		
		Atom[] atoms = getAtoms(hetAtoms, cofactorSizeToUse);
		AsaCalculator asaCalc = new AsaCalculator(atoms, 
				AsaCalculator.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);
		
		double[] complexAsas = asaCalc.calculateAsas();
		
		if (complexAsas.length!=asas1.length+asas2.length) 
			throw new IllegalArgumentException("The size of ASAs of complex doesn't match that of ASAs 1 + ASAs 2");
		
		
		groupAsas1 = new TreeMap<ResidueNumber, GroupAsa>();
		groupAsas2 = new TreeMap<ResidueNumber, GroupAsa>();
		
		
		this.interfaceArea = 0;
		
		for (int i=0;i<asas1.length;i++) {			
			this.interfaceArea += (asas1[i] - complexAsas[i]);
			
			Group g = atoms[i].getGroup();
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
			this.interfaceArea += (asas2[i] - complexAsas[i+asas1.length]);
			
			Group g = atoms[i+asas1.length].getGroup();
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
		this.interfaceArea = this.interfaceArea/2.0;
		 
	}
	
	protected Atom[] getFirstAtoms(boolean hetAtoms, int cofactorSizeToUse) {
		// TODO we need to take hetAtoms into account (those within peptide chain or nucleotide chain)
		// TODO we should do this based on the cofactor size cutoff given
		boolean  USE_HETATOMS_FOR_ASA = false;
				
		Atom[] atoms1 = StructureTools.getAllNonHAtomArray(chains.getFirst(), USE_HETATOMS_FOR_ASA);		
		
		return atoms1;		
	}
	
	protected Atom[] getSecondAtoms(boolean hetAtoms, int cofactorSizeToUse) {
		// TODO we need to take hetAtoms into account (those within peptide chain or nucleotide chain)
		// TODO we should do this based on the cofactor size cutoff given
		boolean  USE_HETATOMS_FOR_ASA = false;
						
		Atom[] atoms2 = StructureTools.getAllNonHAtomArray(chains.getSecond(), USE_HETATOMS_FOR_ASA);
		
		return atoms2;
	}
	
	protected Atom[] getAtoms(boolean hetAtoms, int cofactorSizeToUse) {
		Atom[] atoms1 = getFirstAtoms(hetAtoms, cofactorSizeToUse);
		Atom[] atoms2 = getSecondAtoms(hetAtoms, cofactorSizeToUse);
		
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

	public Map<ResidueNumber, GroupAsa> getFirstGroupAsas() {
		return groupAsas1;
	}
	
	public GroupAsa getFirstGroupAsa(ResidueNumber resNum) {
		return groupAsas1.get(resNum);
	}
	
	public Map<ResidueNumber, GroupAsa> getSecondGroupAsas() {
		return groupAsas2;
	}

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
		return (Double.compare(o.interfaceArea,this.interfaceArea));
	}
	
}
