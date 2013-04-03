/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * http://www.ncbi.nlm.nih.gov/pubmed/14534198
 * </pre>
 * 
 * Thanks to Yuzhen Ye and A. Godzik for granting permission to freely use and redistribute this code.
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
 *
 * Created on Jun 17, 2009
 * Created by Andreas Prlic - RCSB PDB 
 * 
 */

package org.biojava.bio.structure.align.fatcat.calc;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

/** a class that performs calculations on AFPCHains
 * 
 * @author Andreas Prlic
 *
 */
public class AFPCalculator
{
	public static final boolean debug = FatCatAligner.debug;


	public static final  void extractAFPChains(FatCatParameters params, AFPChain afpChain,Atom[] ca1,Atom[] ca2){



		List<AFP> afpSet = new ArrayList<AFP>();
		afpChain.setAfpSet(afpSet);

		if ( debug )
			System.err.println("nr of atoms ca1: " + ca1.length + " ca2: " +  ca2.length);



		int     p1, p2;
		int n0, n, n1, n2;
		double  filter1;
		double rmsd = 0;
		//double[] r = new double[9]; // rotation matrix
		//double[] t = new double[3]; // shift vector

		Matrix r = new Matrix(3,3);
		Atom   t = new AtomImpl();


		int sparse = params.getSparse();
		int maxTra = params.getMaxTra();
		int fragLen = params.getFragLen();
		double disFilter = params.getDisFilter();
		double rmsdCut = params.getRmsdCut();
		double badRmsd = params.getBadRmsd();
		double fragScore = params.getFragScore();

		int     add = sparse + 1; //if add > 1, use sparse sampling
		n0 = n = n1 = n2 = 0;

		int minLen = 0;

		int prot1Length = ca1.length;
		int prot2Length = ca2.length;

		if(prot1Length < prot2Length)
			minLen = prot1Length;
		else
			minLen = prot2Length;
		afpChain.setMinLen(minLen);

		afpChain.setBlockResList(new int[maxTra+1][2][minLen]);
		afpChain.setFocusRes1(new int[minLen]);
		afpChain.setFocusRes2(new int[minLen]);

		for(p1 = 0; p1 < prot1Length - fragLen; p1 += add )    {
			for(p2 = 0; p2 < prot2Length - fragLen; p2 += add)     {
				n0 ++;
				filter1 = getEnd2EndDistance(ca1, ca2, p1, p1 + fragLen - 1, p2, p2 + fragLen - 1);
				//difference bewteen end-to-end distances
				if(filter1 > disFilter) { n1 ++; continue; }
				boolean filter2 = filterTerminal(ca1,ca2, p1, p1 + fragLen - 1, p2, p2 + fragLen - 1, fragLen, minLen);
				if(filter2)     {
					n2 ++;
					continue;

				} //be cautious to use this filter !!

				// here FATCAT does a a jacobi transformation
				//rmsd = kearsay(fragLen, ca1[p1], ca2[p2], r, t);
				// we use the BioJava SVD instead...

				//
				rmsd = getRmsd(ca1,ca2,fragLen, p1,p2,r,t);

				//printf("afp %d: p1 %d p2 %d rmsd %f end-to-end dis %f\n", afpSet.size(), p1, p2, rmsd, filter1);

				if(rmsd < rmsdCut)      {
					AFP     afptmp = new AFP();
					afptmp.setP1(p1);
					afptmp.setP2(p2);
					afptmp.setFragLen(fragLen);
					afptmp.setRmsd(rmsd);
					afptmp.setM(r);
					afptmp.setT(t.getCoords());
					afptmp.setScore(scoreAfp(afptmp,badRmsd,fragScore));
					afpSet.add(afptmp);
					n ++;
				}
			}
		}

		int afpNum = afpSet.size();

		if(debug) {
			String msg = String.format("possible AFP-pairs %d, remain %d after filter 1 remove %d; filter 2 remove %d\n",
					n0, afpNum, n1, n2);
			System.err.println(msg);
		}


	}

	/**
	 * filter 1 for AFP extration: the distance of end-to-end
	 * @param p1b
	 * @param p1e
	 * @param p2b
	 * @param p2e
	 * @return
	 */
	private static final double getEnd2EndDistance(Atom[] ca1, Atom[] ca2, int p1b, int p1e, int p2b, int p2e)
	{

		double min = 99;
		try {
			double dist1 = Calc.getDistance(ca1[p1b], ca1[p1e]);
			double dist2 = Calc.getDistance(ca2[p2b], ca2[p2e]);
			min = dist1 - dist2;

		} catch (Exception e){
			e.printStackTrace();
		}
		return Math.abs(min);
	}

	/**
	 * filter 2 for AFP extration: the context
	 * @param p1b
	 * @param p1e
	 * @param p2b
	 * @param p2e
	 * @return
	 */

	private static final  boolean filterTerminal(Atom[] ca1, Atom[] ca2, int p1b, int p1e, int p2b, int p2e, int fragLen, int minLen)
	{
		int     d1 = (p1b < p2b)?p1b:p2b;
		int     d2 = (ca1.length - p1e) < (ca2.length - p2e)?(ca1.length - p1e):(ca2.length - p2e);
		int     d3 = d1 + d2 + fragLen; //maximum alignment length from current AFP


		/// DO NOT DO Math.round() this will give different results to FATCAT....
		int     d4 = (int)(0.3 * minLen);

		if(d3 < d4)     {
			//if (debug){
			//   System.out.println("filterTerminal: " + d3 + " " + d4 +" " + p1b + " " + p1e + " " + p2b + " " +  p2e + " " + fragLen + " " + minLen );
			//}
			return true;
		}

		return false;
	}

	private static final double getRmsd(Atom[] ca1, Atom[] ca2, int fragLen, int p1, int p2, Matrix m, Atom t) {


		double rmsd = 99.9;
		try {
			Atom[] catmp1 = getFragment(ca1, p1, fragLen,false);
			Atom[] catmp2 = getFragment(ca2, p2, fragLen,true); // clone the atoms for fragment 2 so we can manipulate them...

			if ( catmp1 == null) {
				System.err.println("could not get fragment for ca1 " + p1 + " " + fragLen );
				return rmsd;

			}

			if ( catmp2 == null) {
				System.err.println("could not get fragment for ca2 " + p2 + " " + fragLen );
				return rmsd;

			}

			SVDSuperimposer svd = new SVDSuperimposer(catmp1, catmp2);

			m = svd.getRotation();
			t = svd.getTranslation();

			for (Atom a : catmp2){
				Calc.rotate(a,m);
				Calc.shift(a,t);

			}

			rmsd = SVDSuperimposer.getRMS(catmp1,catmp2);

		} catch (Exception e){
			e.printStackTrace();
		}

		return rmsd;
	}

	/** get a continue subset of Atoms based by the starting position and the length
	 *
	 * @param caall
	 * @param pos ... the start position
	 * @param fragmentLength .. the length of the subset to extract.
	 * @param clone: returns a copy of the atom (in case the coordinate get manipulated...)
	 * @return an Atom[] array
	 */
	private static final Atom[] getFragment(Atom[] caall, int pos, int fragmentLength , boolean clone){

		if ( pos+fragmentLength > caall.length)
			return null;

		Atom[] tmp = new Atom[fragmentLength];

		for (int i=0;i< fragmentLength;i++){
			if (clone){
				tmp[i] = (Atom)caall[i+pos].clone();
			} else {
				tmp[i] = (Atom)caall[i+pos];
			}
		}
		return tmp;

	}


	/**
	 * Assign score to each AFP
	 */

	private static final double scoreAfp(AFP afp, double badRmsd, double fragScore)
	{
		//longer AFP with low rmsd is better
		double  s, w;
		//s = (rmsdCut - afptmp.rmsd) * afptmp.len; //the same scroing strategy as that in the post-processing
		w = afp.getRmsd() / badRmsd;
		w = w * w;
		s = fragScore * (1.0 - w);
		return s;
	}

	//------------------------------------------------------------------
	//Sort the AFPs in increase of their diagonals(i,j)
	//------------------------------------------------------------------
	public static final  void sortAfps(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
	{


		List<AFP> afpSet = afpChain.getAfpSet();

		if ( debug)
			System.err.println("entering sortAfps");

		int pro1Len = ca1.length;
		int pro2Len = ca2.length;

		afpChain.setAfpIndex(      new int[pro1Len][pro2Len]); //the index of (i,j) pair in AFP list, otherwise -1
		afpChain.setAfpAftIndex(   new int[pro1Len][pro2Len]);  //the index of AFP (i,j*) nearest to (i,j), j*<j. if a AFP exits for (i,j), it equals to afpIndex
		afpChain.setAfpBefIndex(   new int[pro1Len][pro2Len]); //the index of AFP (i,j*) nearest to (i,j), j*>j. if a AFP exits for (i,j), it equals to afpIndex

		int[][] afpIndex       = afpChain.getAfpIndex();
		int[][] afpAftIndex    = afpChain.getAfpAftIndex();
		int[][] afpBefIndex    = afpChain.getAfpBefIndex();

		for(int i = 0; i < pro1Len; i ++)   {
			for(int j = 0; j < pro2Len; j ++)   {

				afpIndex[i][j] = afpAftIndex[i][j] = afpBefIndex[i][j] = -1;
			}
		}

		//index the AFP for easy extraction of compatible AFPs
		int afpNum = afpSet.size();

		int     b0 = 0;
		for(int a = 0; a < afpNum; a ++)    {
			if(a == afpNum - 1 || afpSet.get(a).getP1() != afpSet.get(a+1).getP1())   {
				int i = afpSet.get(a).getP1();
				for(int b = b0; b <= a; b ++)       {
					int j = afpSet.get(b).getP2();
					afpIndex[i][j]=b ;
					afpBefIndex[i][j]=b;
					afpAftIndex[i][j]=b;
					if(afpSet.get(b).getP1() != i)    {
						System.err.println(String.format("Warning: wrong afp index %d %d\n", i, afpSet.get(b).getP1()));
						return;
					}
				}
				for(int k = 1; k < pro2Len; k ++)   {
					if( afpBefIndex[i][k] == -1){
						afpBefIndex[i][k] = afpBefIndex[i][k-1];
					}
				}
				for(int k = pro2Len - 2; k >= 0; k --)      {
					if(afpAftIndex[i][k] == -1) {
						afpAftIndex[i][k] =  afpAftIndex[i][k+1];
					}
				}
				b0 = a + 1;
			}
		}

		if ( debug)
			System.err.println("done sortAfps");


	}
}
