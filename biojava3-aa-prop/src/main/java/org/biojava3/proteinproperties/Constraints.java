package org.biojava3.proteinproperties;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

public class Constraints {
	private static AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
	public static Map<AminoAcidCompound, Double> aa2MolecularWeight = new HashMap<AminoAcidCompound, Double>();
	public static Map<AminoAcidCompound, Double> aa2Hydrophathicity = new HashMap<AminoAcidCompound, Double>();
	public static Map<AminoAcidCompound, Double> aa2PKa = new HashMap<AminoAcidCompound, Double>();
	
	static{
		initMolecularWeight();
		initHydropathicity();
		initPKa();
	}
	
	private static void initMolecularWeight(){
		/*
		 * Initialize molecular weights - http://au.expasy.org/tools/findmod/findmod_masses.html#AA
		 */
//		Alanine (A)	71.03711	71.0788
		aa2MolecularWeight.put(aaSet.getCompoundForString("A"), 71.0788);
//		Arginine (R)	156.10111	156.1875
		aa2MolecularWeight.put(aaSet.getCompoundForString("R"), 156.1875);
//		Asparagine (N)	114.04293	114.1038
		aa2MolecularWeight.put(aaSet.getCompoundForString("N"), 114.1038);
//		Aspartic acid (D)	115.02694	115.0886
		aa2MolecularWeight.put(aaSet.getCompoundForString("D"), 115.0886);
//		Cysteine (C)	103.00919	103.1388
		aa2MolecularWeight.put(aaSet.getCompoundForString("C"), 103.1388);
//		Glutamic acid (E)	129.04259	129.1155
		aa2MolecularWeight.put(aaSet.getCompoundForString("E"), 129.1155);
//		Glutamine (Q)	128.05858	128.1307
		aa2MolecularWeight.put(aaSet.getCompoundForString("Q"), 128.1307);
//		Glycine (G)	57.02146	57.0519
		aa2MolecularWeight.put(aaSet.getCompoundForString("G"), 57.0519);
//		Histidine (H)	137.05891	137.1411
		aa2MolecularWeight.put(aaSet.getCompoundForString("H"), 137.1411);
//		Isoleucine (I)	113.08406	113.1594
		aa2MolecularWeight.put(aaSet.getCompoundForString("I"), 113.1594);
//		Leucine (L)	113.08406	113.1594
		aa2MolecularWeight.put(aaSet.getCompoundForString("L"), 113.1594);
//		Lysine (K)	128.09496	128.1741
		aa2MolecularWeight.put(aaSet.getCompoundForString("K"), 128.1741);
//		Methionine (M)	131.04049	131.1926
		aa2MolecularWeight.put(aaSet.getCompoundForString("M"), 131.1926);
//		Phenylalanine (F)	147.06841	147.1766
		aa2MolecularWeight.put(aaSet.getCompoundForString("F"), 147.1766);
//		Proline (P)	97.05276	97.1167
		aa2MolecularWeight.put(aaSet.getCompoundForString("P"), 97.1167);
//		Serine (S)	87.03203	87.0782
		aa2MolecularWeight.put(aaSet.getCompoundForString("S"), 87.0782);
//		Threonine (T)	101.04768	101.1051
		aa2MolecularWeight.put(aaSet.getCompoundForString("T"), 101.1051);
//		Tryptophan (W)	186.07931	186.2132
		aa2MolecularWeight.put(aaSet.getCompoundForString("W"), 186.2132);
//		Tyrosine (Y)	163.06333	163.1760
		aa2MolecularWeight.put(aaSet.getCompoundForString("Y"), 163.1760);
//		Valine (V)	99.06841	99.1326
		aa2MolecularWeight.put(aaSet.getCompoundForString("V"), 99.1326);
	}
	
	private static void initHydropathicity(){
		/*
		 * Initialize hydropathicity - http://au.expasy.org/tools/pscale/Hphob.Doolittle.html
		 */
//		Ala(A):  1.800  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("A"), 1.800);
//		Arg(R): -4.500  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("R"), -4.500);
//		Asn(N): -3.500  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("N"), -3.500);
//		Asp(D): -3.500  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("D"), -3.500);
//		Cys(C):  2.500  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("C"), 2.500);
//		Gln(E): -3.500  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("E"), -3.500);
//		Glu(Q): -3.500  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("Q"), -3.500);
//		Gly(G): -0.400  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("G"), -0.400);
//		His(H): -3.200  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("H"), -3.200);
//		Ile(I):  4.500  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("I"), 4.500);
//		Leu(L):  3.800  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("L"), 3.800);
//		Lys(K): -3.900  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("K"), -3.900);
//		Met(M):  1.900  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("M"), 1.900);
//		Phe(F):  2.800  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("F"), 2.800);
//		Pro(P): -1.600  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("P"), -1.600);
//		Ser(S): -0.800  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("S"), -0.800);
//		Thr(T): -0.700  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("T"), -0.700);
//		Trp(W): -0.900  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("W"), -0.900);
//		Tyr(Y): -1.300  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("Y"), -1.300);
//		Val(V):  4.200  
		aa2Hydrophathicity.put(aaSet.getCompoundForString("V"), 4.200);
	}

	private static void initPKa(){
//		K, Lys	10.5	
		aa2PKa.put(aaSet.getCompoundForString("K"), 10.5);
//		D, Asp	3.86
		aa2PKa.put(aaSet.getCompoundForString("D"), 3.86);
//		R, Arg	12.4	
		aa2PKa.put(aaSet.getCompoundForString("R"), 12.4);
//		E, Glu	4.25
		aa2PKa.put(aaSet.getCompoundForString("E"), 4.25);
//		H, His	6.00	
		aa2PKa.put(aaSet.getCompoundForString("H"), 6.00);
//		C, Cys	8.33
		aa2PKa.put(aaSet.getCompoundForString("C"), 8.33);
//		Y, Tyr	10.0
		aa2PKa.put(aaSet.getCompoundForString("Y"), 10.0);
	}
}
