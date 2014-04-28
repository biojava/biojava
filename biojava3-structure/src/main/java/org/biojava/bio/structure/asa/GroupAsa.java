package org.biojava.bio.structure.asa;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;

/**
 * A class to store the results of ASA calculations, it can
 * hold ASA values per atom present in {@link org.biojava.bio.structure.Group}
 * 
 * @author duarte_j
 *
 */
public class GroupAsa {

	
	// ASA in extended tripeptide conformation (GLY-X-GLY) 
	private static final HashMap<Character,Double> tripeptAsa = initTriPeptAsas();
	
	private static HashMap<Character,Double>  initTriPeptAsas() {
		// ASA in extended tripeptide conformation (GLY-X-GLY) from Miller et al JMB 1987 (for calculation of relative ASAs)
		HashMap<Character,Double> map = new HashMap<Character,Double>();
		map.put('A', 113.0);
		map.put('R', 241.0);
		map.put('N', 158.0);
		map.put('D', 151.0);
		map.put('C', 140.0);
		map.put('Q', 189.0);
		map.put('E', 183.0);
		map.put('G',  85.0);
		map.put('H', 194.0);
		map.put('I', 182.0);
		map.put('L', 180.0);
		map.put('K', 211.0);
		map.put('M', 204.0);
		map.put('F', 218.0);
		map.put('P', 143.0);
		map.put('S', 122.0);
		map.put('T', 146.0);
		map.put('W', 259.0);
		map.put('Y', 229.0);
		map.put('V', 160.0);
		return map;
	}
	
	
	
	
	private Group g;
	
	/**
	 * ASA of uncomplexed residue
	 */
	private double asaU;
	
	/**
	 * ASA of complexed residue
	 */
	private double asaC;
	
	/**
	 * The individual atoms uncomplexed ASAs
	 */
	private List<Double> atomAsaUs;
	
	/**
	 * The individual atoms complexed ASAs
	 */
	private List<Double> atomAsaCs;
	
	public GroupAsa(Group g) {
		this.g = g;
		
		int groupNoHSize = getGroupNoHSize();
		atomAsaUs = new ArrayList<Double>(groupNoHSize);
		atomAsaCs = new ArrayList<Double>(groupNoHSize);
	}
	
	private int getGroupNoHSize() {
		int count = 0;
		for (Atom atom:g.getAtoms()) {
			if (atom.getElement()!=Element.H) count++;
		}
		return count;
	}

	public Group getGroup() {
		return g;
	}
	
	/**
	 * Returns the ASA of the residue in the uncomplexed state 
	 * @return
	 */
	public double getAsaU() {
		return asaU;
	}

	public void setAsaU(double asaU) {
		this.asaU = asaU;
	}

	/**
	 * Returns the ASA of the residue in the complexed state
	 * @return
	 */
	public double getAsaC() {
		return asaC;
	}

	public void setAsaC(double asaC) {
		this.asaC = asaC;
	}
	
	public void addAtomAsaU(double asa) {
		this.asaU += asa;
		this.atomAsaUs.add(asa);
	}

	public void addAtomAsaC(double asa) {
		this.asaC += asa;
		this.atomAsaCs.add(asa);
	}
	
	/**
	 * Returns the BSA value for this group, i.e. the difference between ASA uncomplexed and ASA complexed
	 * @return
	 */
	public double getBsa() {
		return (asaU-asaC);
	}

	/**
	 * Returns the bsa/asa(uncomplexed) ratio, i.e. the ratio of burial of a residue upon complexation
	 * @return
	 */
	public double getBsaToAsaRatio() {
		return getBsa()/asaU;
	}
	
	/**
	 * Returns the relative (uncomplexed) ASA, i.e. the ASA of the residue 
	 * with respect to its ASA in an extended tri-peptide conformation (GLY-x-GLY)
	 * @return
	 */
	public double getRelativeAsaU() {
		if (!g.getType().equals(GroupType.AMINOACID)) 
			throw new IllegalArgumentException("Can not calculate relative ASA for non amino-acid");

		char aa = ((AminoAcid)g).getAminoType();

		return (asaU/tripeptAsa.get(aa));

	}
	
	/**
	 * Returns the relative (complexed) ASA, i.e. the ASA of the residue 
	 * with respect to its ASA in an extended tri-peptide conformation (GLY-x-GLY) 
	 * @return
	 */
	public double getRelativeAsaC() {
		if (!g.getType().equals(GroupType.AMINOACID)) 
			throw new IllegalArgumentException("Can not calculate relative ASA for non amino-acid");

		char aa = ((AminoAcid)g).getAminoType();

		return (asaC/tripeptAsa.get(aa));

	}
	
}
