package org.biojava3.aaproperties;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.aaproperties.PeptideProperties.SingleLetterAACode;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

/**
 * This class is used to support the implementation of properties stated in IPeptideProperties. 
 * It initializes several values that would be needed for the computation of properties such as
 * <p/>
 * Molecular weight<br/>
 * Instability index<br/>
 * Hydropathy value<br/>
 * pKa<br/>
 * 
 * @author kohchuanhock
 * @version 2011.05.21
 * @see IPeptideProperties
 */
public class Constraints {
	private static AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
	//A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V
	public static AminoAcidCompound A = aaSet.getCompoundForString("A");
	public static AminoAcidCompound R = aaSet.getCompoundForString("R");
	public static AminoAcidCompound N = aaSet.getCompoundForString("N");
	public static AminoAcidCompound D = aaSet.getCompoundForString("D");
	public static AminoAcidCompound C = aaSet.getCompoundForString("C");
	public static AminoAcidCompound E = aaSet.getCompoundForString("E");
	public static AminoAcidCompound Q = aaSet.getCompoundForString("Q");
	public static AminoAcidCompound G = aaSet.getCompoundForString("G");
	public static AminoAcidCompound H = aaSet.getCompoundForString("H");
	public static AminoAcidCompound I = aaSet.getCompoundForString("I");
	public static AminoAcidCompound L = aaSet.getCompoundForString("L");
	public static AminoAcidCompound K = aaSet.getCompoundForString("K");
	public static AminoAcidCompound M = aaSet.getCompoundForString("M");
	public static AminoAcidCompound F = aaSet.getCompoundForString("F");
	public static AminoAcidCompound P = aaSet.getCompoundForString("P");
	public static AminoAcidCompound S = aaSet.getCompoundForString("S");
	public static AminoAcidCompound T = aaSet.getCompoundForString("T");
	public static AminoAcidCompound W = aaSet.getCompoundForString("W");
	public static AminoAcidCompound Y = aaSet.getCompoundForString("Y");
	public static AminoAcidCompound V = aaSet.getCompoundForString("V");

	public static Map<AminoAcidCompound, Double> aa2ExtinctionCoefficient = new HashMap<AminoAcidCompound, Double>();
	public static Map<AminoAcidCompound, Double> aa2MolecularWeight = new HashMap<AminoAcidCompound, Double>();
	public static Map<AminoAcidCompound, Double> aa2Hydrophathicity = new HashMap<AminoAcidCompound, Double>();
	public static Map<AminoAcidCompound, Double> aa2PKa = new HashMap<AminoAcidCompound, Double>();
	public static Map<String, Double> diAA2Instability = new HashMap<String, Double>();

	public static Map<AminoAcidCompound, Double> aa2NTerminalPka = new HashMap<AminoAcidCompound, Double>();
	public static Map<AminoAcidCompound, Double> aa2CTerminalPka = new HashMap<AminoAcidCompound, Double>();
	
	static{
		initMolecularWeight();
		initHydropathicity();
		initPKa();
		initInstability();
		initExtinctionCoefficient();
	}

	/** 
	 * Does the initialization of molecular weights based on http://au.expasy.org/tools/findmod/findmod_masses.html#AA
	 */
	public static void initMolecularWeight(){
		//		Alanine (A)	71.03711	71.0788
		aa2MolecularWeight.put(A, 71.0788);
		//		Arginine (R)	156.10111	156.1875
		aa2MolecularWeight.put(R, 156.1875);
		//		Asparagine (N)	114.04293	114.1038
		aa2MolecularWeight.put(N, 114.1038);
		//		Aspartic acid (D)	115.02694	115.0886
		aa2MolecularWeight.put(D, 115.0886);
		//		Cysteine (C)	103.00919	103.1388
		aa2MolecularWeight.put(C, 103.1388);
		//		Glutamic acid (E)	129.04259	129.1155
		aa2MolecularWeight.put(E, 129.1155);
		//		Glutamine (Q)	128.05858	128.1307
		aa2MolecularWeight.put(Q, 128.1307);
		//		Glycine (G)	57.02146	57.0519
		aa2MolecularWeight.put(G, 57.0519);
		//		Histidine (H)	137.05891	137.1411
		aa2MolecularWeight.put(H, 137.1411);
		//		Isoleucine (I)	113.08406	113.1594
		aa2MolecularWeight.put(I, 113.1594);
		//		Leucine (L)	113.08406	113.1594
		aa2MolecularWeight.put(L, 113.1594);
		//		Lysine (K)	128.09496	128.1741
		aa2MolecularWeight.put(K, 128.1741);
		//		Methionine (M)	131.04049	131.1926
		aa2MolecularWeight.put(M, 131.1926);
		//		Phenylalanine (F)	147.06841	147.1766
		aa2MolecularWeight.put(F, 147.1766);
		//		Proline (P)	97.05276	97.1167
		aa2MolecularWeight.put(P, 97.1167);
		//		Serine (S)	87.03203	87.0782
		aa2MolecularWeight.put(S, 87.0782);
		//		Threonine (T)	101.04768	101.1051
		aa2MolecularWeight.put(T, 101.1051);
		//		Tryptophan (W)	186.07931	186.2132
		aa2MolecularWeight.put(W, 186.2132);
		//		Tyrosine (Y)	163.06333	163.1760
		aa2MolecularWeight.put(Y, 163.1760);
		//		Valine (V)	99.06841	99.1326
		aa2MolecularWeight.put(V, 99.1326);
	}

	/**
	 * Does the initialization of hydropathicity based on http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html
	 */
	private static void initHydropathicity(){
		//		Ala(A):  1.800  
		aa2Hydrophathicity.put(A, 1.800);
		//		Arg(R): -4.500  
		aa2Hydrophathicity.put(R, -4.500);
		//		Asn(N): -3.500  
		aa2Hydrophathicity.put(N, -3.500);
		//		Asp(D): -3.500  
		aa2Hydrophathicity.put(D, -3.500);
		//		Cys(C):  2.500  
		aa2Hydrophathicity.put(C, 2.500);
		//		Gln(E): -3.500  
		aa2Hydrophathicity.put(E, -3.500);
		//		Glu(Q): -3.500  
		aa2Hydrophathicity.put(Q, -3.500);
		//		Gly(G): -0.400  
		aa2Hydrophathicity.put(G, -0.400);
		//		His(H): -3.200  
		aa2Hydrophathicity.put(H, -3.200);
		//		Ile(I):  4.500  
		aa2Hydrophathicity.put(I, 4.500);
		//		Leu(L):  3.800  
		aa2Hydrophathicity.put(L, 3.800);
		//		Lys(K): -3.900  
		aa2Hydrophathicity.put(K, -3.900);
		//		Met(M):  1.900  
		aa2Hydrophathicity.put(M, 1.900);
		//		Phe(F):  2.800  
		aa2Hydrophathicity.put(F, 2.800);
		//		Pro(P): -1.600  
		aa2Hydrophathicity.put(P, -1.600);
		//		Ser(S): -0.800  
		aa2Hydrophathicity.put(S, -0.800);
		//		Thr(T): -0.700  
		aa2Hydrophathicity.put(T, -0.700);
		//		Trp(W): -0.900  
		aa2Hydrophathicity.put(W, -0.900);
		//		Tyr(Y): -1.300  
		aa2Hydrophathicity.put(Y, -1.300);
		//		Val(V):  4.200  
		aa2Hydrophathicity.put(V, 4.200);
	}

	/**
	 * Does the initialization of PKa based on
	 * http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#NetCharge
	 */
	private static void initPKaInnovagen(){
		/*
		 * A.Lehninger, Principles of Biochemistry, 4th Edition (2005), Chapter 3, page78, Table 3-1.
		 */
		//(NH2-)	9.69	(-COOH)	2.34
		aa2CTerminalPka.put(G, 2.34);
		aa2CTerminalPka.put(A, 2.34);
		aa2CTerminalPka.put(P, 1.99);
		aa2CTerminalPka.put(V, 2.32);
		aa2CTerminalPka.put(L, 2.36);
		aa2CTerminalPka.put(I, 2.36);
		aa2CTerminalPka.put(M, 2.28);
		
		aa2CTerminalPka.put(F, 1.83);
		aa2CTerminalPka.put(Y, 2.20);
		aa2CTerminalPka.put(W, 2.38);
		
		aa2CTerminalPka.put(S, 2.21);
		aa2CTerminalPka.put(T, 2.11);
		aa2CTerminalPka.put(C, 1.96);
		aa2CTerminalPka.put(N, 2.02);
		aa2CTerminalPka.put(Q, 2.17);
		
		aa2CTerminalPka.put(K, 2.18);
		aa2CTerminalPka.put(H, 1.82);
		aa2CTerminalPka.put(R, 2.17);
		
		aa2CTerminalPka.put(D, 1.88);
		aa2CTerminalPka.put(E, 2.19);
		
		aa2NTerminalPka.put(G, 9.60);
		aa2NTerminalPka.put(A, 9.69);
		aa2NTerminalPka.put(P, 10.96);
		aa2NTerminalPka.put(V, 9.62);
		aa2NTerminalPka.put(L, 9.60);
		aa2NTerminalPka.put(I, 9.68);
		aa2NTerminalPka.put(M, 9.21);
		
		aa2NTerminalPka.put(F, 9.13);
		aa2NTerminalPka.put(Y, 9.11);
		aa2NTerminalPka.put(W, 9.39);
		
		aa2NTerminalPka.put(S, 9.15);
		aa2NTerminalPka.put(T, 9.62);
		aa2NTerminalPka.put(C, 10.28);
		aa2NTerminalPka.put(N, 8.80);
		aa2NTerminalPka.put(Q, 9.13);
		
		aa2NTerminalPka.put(K, 8.95);
		aa2NTerminalPka.put(H, 9.17);
		aa2NTerminalPka.put(R, 9.04);
		
		aa2NTerminalPka.put(D, 9.60);
		aa2NTerminalPka.put(E, 9.67);
		
		//		K, Lys	10.53	
		aa2PKa.put(K, 10.53);
		//		D, Asp	3.65
		aa2PKa.put(D, 3.65);
		//		R, Arg	12.48	
		aa2PKa.put(R, 12.48);
		//		E, Glu	4.25
		aa2PKa.put(E, 4.25);
		//		H, His	6.00	
		aa2PKa.put(H, 6.00);
		//		C, Cys	8.18
		aa2PKa.put(C, 8.18);
		//		Y, Tyr	10.07
		aa2PKa.put(Y, 10.07);
	}
	
	private static void initPKa(){
		initPKaInnovagen();
	}

	/** 
	 * Does the initialization of dipeptide instability index based on the following paper
	 * 
	 * Guruprasad, K., Reddy, B.V.B. and Pandit, M.W. (1990) 
	 * Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. 
	 * Protein Eng. 4,155-161. Table III.
	 */
	private static void initInstability(){
		double[][] instability = {
				//W		C		M		H		Y		F		Q		N		I		R		D		P		T		K		E		V		S		G		A		L
				{1.0,	1.0,	24.68,	24.68,	1.0,	1.0,	1.0,	13.34,	1.0,	1.0,	1.0,	1.0,	-14.03, 1.0, 	1.0, 	-7.49, 	1.0, 	-9.37, 	-14.03, 13.34},
				{24.68, 1.0, 	33.6, 	33.6, 	1.0, 	1.0, 	-6.54, 	1.0, 	1.0, 	1.0, 	20.26,	20.26, 	33.6, 	1.0, 	1.0, 	-6.54, 	1.0, 	1.0, 	1.0, 	20.26},
				{1.0, 	1.0, 	-1.88, 	58.28, 	24.68, 	1.0, 	-6.54, 	1.0, 	1.0, 	-6.54, 	1.0, 	44.94, 	-1.88, 	1.0, 	1.0, 	1.0, 	44.94, 	1.0, 	13.34, 	1.0},
				{-1.88, 1.0, 	1.0,	1.0,	44.94,	-9.37,	1.0,	24.68,	44.94,	1.0,	1.0,	-1.88,	-6.54,	24.68,	1.0,	1.0,	1.0,	-9.37,	1.0,	1.0},
				{-9.37,	1.0,	44.94,	13.34,	13.34,	1.0,	1.0,	1.0,	1.0,	-15.91,	24.68,	13.34,	-7.49,	1.0,	-6.54,	1.0,	1.0,	-7.49,	24.68,	1.0},
				{1.0,	1.0,	1.0,	1.0,	33.6,	1.0,	1.0,	1.0,	1.0,	1.0,	13.34,	20.26,	1.0,	-14.03,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0},
				{1.0,	-6.54,	1.0,	1.0,	-6.54,	-6.54,	20.26,	1.0,	1.0,	1.0,	20.26,	20.26,	1.0,	1.0,	20.26,	-6.54,	44.94,	1.0,	1.0,	1.0},
				{-9.37,	-1.88,	1.0,	1.0,	1.0,	-14.03,	-6.54,	1.0,	44.94,	1.0,	1.0,	-1.88,	-7.49,	24.68,	1.0,	1.0,	1.0,	-14.03,	1.0,	1.0},
				{1.0,	1.0,	1.0,	13.34,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	-1.88,	1.0,	-7.49,	44.94,	-7.49,	1.0,	1.0,	1.0,	20.26},
				{58.28,	1.0,	1.0,	20.26,	-6.54,	1.0,	20.26,	13.34,	1.0,	58.28,	1.0,	20.26,	1.0,	1.0,	1.0,	1.0,	44.94,	-7.49,	1.0,	1.0},
				{1.0,	1.0,	1.0,	1.0,	1.0,	-6.54,	1.0,	1.0,	1.0,	-6.54,	1.0,	1.0,	-14.03,	-7.49,	1.0,	1.0,	20.26,	1.0,	1.0,	1.0},
				{-1.88,	-6.54,	-6.54,	1.0,	1.0,	20.26,	20.26,	1.0,	1.0,	-6.54,	-6.54,	20.26,	1.0,	1.0,	18.38,	20.26,	20.26,	1.0,	20.26,	1.0},
				{-14.03,1.0,	1.0,	1.0,	1.0,	13.34,	-6.54,	-14.03,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	20.26,	1.0,	1.0,	-7.49,	1.0,	1.0},
				{1.0,	1.0,	33.6,	1.0,	1.0,	1.0,	24.68,	1.0,	-7.49,	33.6,	1.0,	-6.54,	1.0,	1.0,	1.0,	-7.49,	1.0,	-7.49,	1.0,	-7.49},
				{-14.03,44.94,	1.0,	-6.54,	1.0,	1.0,	20.26,	1.0,	20.26,	1.0,	20.26,	20.26,	1.0,	1.0,	33.6,	1.0,	20.26,	1.0,	1.0,	1.0},
				{1.0,	1.0,	1.0,	1.0,	-6.54,	1.0,	1.0,	1.0,	1.0,	1.0,	-14.03,	20.26,	-7.49,	-1.88,	1.0,	1.0,	1.0,	-7.49,	1.0,	1.0},
				{1.0,	33.6,	1.0,	1.0,	1.0,	1.0,	20.26,	1.0,	1.0,	20.26,	1.0,	44.94,	1.0,	1.0,	20.26,	1.0,	20.26,	1.0,	1.0,	1.0},
				{13.34,	1.0,	1.0,	1.0,	-7.49,	1.0,	1.0,	-7.49,	-7.49,	1.0,	1.0,	1.0,	-7.49,	-7.49,	-6.54,	1.0,	1.0,	13.34,	-7.49,	1.0},
				{1.0,	44.94,	1.0,	-7.49,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	-7.49,	20.26,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0},
				{24.68,	1.0,	1.0,	1.0,	1.0,	1.0,	33.6,	1.0,	1.0,	20.26,	1.0,	20.26,	1.0,	-7.49,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0}
		};

		SingleLetterAACode[] aa = SingleLetterAACode.values();
		for(int i = 0; i < aa.length; i++){
			for(int j = 0; j < aa.length; j++){
				diAA2Instability.put("" + aa[i] + aa[j], instability[i][j]);
			}
		}
	}
	
	/**
	 * Does the initialization of extinction coefficient based on
	 * http://au.expasy.org/tools/protparam-doc.html
	 */
	public static void initExtinctionCoefficient(){
		aa2ExtinctionCoefficient.put(Y, 1490.0);
		aa2ExtinctionCoefficient.put(W, 5500.0);
		aa2ExtinctionCoefficient.put(C, 125.0);
	}
}
