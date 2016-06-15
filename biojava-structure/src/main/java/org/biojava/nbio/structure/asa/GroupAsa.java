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
package org.biojava.nbio.structure.asa;

import org.biojava.nbio.structure.*;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * A class to store the results of ASA calculations, it can
 * hold ASA values per atom present in {@link org.biojava.nbio.structure.Group}
 *
 * @author duarte_j
 *
 */
public class GroupAsa implements Serializable {


	private static final long serialVersionUID = 1L;

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

	public List<Double> getAtomAsaUs() {
		return atomAsaUs;
	}

	public void setAtomAsaUs(List<Double> atomAsaUs) {
		this.atomAsaUs = atomAsaUs;
		this.asaU = 0;
		for (Double atomAsaU : atomAsaUs) {
			this.asaU += atomAsaU;
		}
	}

	public List<Double> getAtomAsaCs() {
		return atomAsaCs;
	}

	public void setAtomAsaCs(List<Double> atomAsaCs) {
		this.atomAsaCs = atomAsaCs;
		this.asaC = 0;
		for (Double atomAsaC : atomAsaCs) {
			this.asaC += atomAsaC;
		}
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

	@Override
	public Object clone() {
		GroupAsa n = new GroupAsa(this.g);
		n.setAsaC(this.getAsaC());
		n.setAsaU(this.getAsaU());
		n.atomAsaUs = new ArrayList<Double>(this.atomAsaUs.size());
		n.atomAsaCs = new ArrayList<Double>(this.atomAsaCs.size());
		for (int i=0;i<this.atomAsaUs.size();i++) {
			n.atomAsaUs.add(this.atomAsaUs.get(i));
		}
		for (int i=0;i<this.atomAsaCs.size();i++) {
			n.atomAsaCs.add(this.atomAsaCs.get(i));
		}

		return n;
	}
}
