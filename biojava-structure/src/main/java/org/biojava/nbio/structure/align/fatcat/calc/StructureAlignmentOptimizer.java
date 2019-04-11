/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * https://www.ncbi.nlm.nih.gov/pubmed/14534198
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

package org.biojava.nbio.structure.align.fatcat.calc;


import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.geometry.SuperPositions;



public class StructureAlignmentOptimizer
{

	//private static final boolean showAlig = false;

	int pro1Len;
	int pro2Len;
	int maxLen;
	Atom[] cod1;
	Atom[] cod2;

	int[][] equSet;
	int equLen;
	int equLen0;
	double[][]sij;

	int maxKeepStep;
	int keepStep;

	double  Dc;   //the criteria for structural equivalent residues, eg. 3.0 (CE), 6.0(ProSup)
	double  rmsdCut;//the criteria for stoping optimization
	double  increase;
	double  stopLenPer;
	double  stopRmsdPer;
	double  stopRmsd;

	double  gapIni;
	double  gapExt;

	double rmsd;

	private static final boolean debug = FatCatAligner.debug;

	/**
	 * optimize the structural alignment by update the equivalent residues
	 * and then run dynamic programming
	 * input: len1 the length of structure 1; c1: the structure information of 1
	 *        len2 the length of structure 2; c2: the structure information of 2
	 *        iniLen and iniSet is the length and list of initial equivalent residues
	 */

	public StructureAlignmentOptimizer(int b1, int end1, Atom[] c1, int b2, int end2, Atom[] c2,
			int iniLen, int[][] iniSet) throws StructureException{

		//System.err.println("optimizing range:" + b1 + "-" + end1 + "("+ (end1-b1) + ") b2:  " + b2 + "-" + end2+ "("+ (end2-b2) + ") iniLen " + iniLen);
		//System.out.println("ca1: " + c1.length + " ca2: " + c2.length);

		int len1 = end1-b1;
		int len2 = end2-b2;

		//System.err.println("len1: " + len1 + " len2:" + len2);

		pro1Len = len1;
		pro2Len = len2;

		cod1 = new Atom[len1];
		cod2 = new Atom[len2];

		for(int i = 0; i < len1; i ++)      {
			Atom a = c1[i+b1];
			//cod1[i] = (Atom)a.clone();
			Group parent = (Group)a.getGroup().clone();
			//cod1[i].setParent(parent);
			cod1[i] = parent.getAtom(a.getName());
			//cod1[i] = c1[i];
		}
		for(int i = 0; i < len2; i ++)      {
			Atom a = c2[i+b2];
			//cod2[i]= (Atom)a.clone();
			Group parent = (Group)a.getGroup().clone();
			//cod2[i].setParent(parent);
			cod2[i] = parent.getAtom(a.getName());
			//cod2[i] = c2[i];
		}


		//initial equivalent sets
		maxLen = (len1 < len2)?len1:len2;

		equSet = new int[2][maxLen];
		for(int i = 0; i < iniLen; i ++)    {

			equSet[0][i] = iniSet[0][i];
			equSet[1][i] = iniSet[1][i];
			if(iniSet[0][i] > len1 || iniSet[1][i] > len2)  {
				throw new RuntimeException(String.format("StructureAlignmentOptimizer: focus exceeds the protein 1 or 2 length!"));
			}
		}

		equLen    = iniLen;
		equLen0   = equLen;

		setParameters();

		sij =  new double[pro1Len][pro2Len];

//      if (showAlig)
//         showCurrentAlignment(iniLen, equSet, "initial alignment");

	}

//   private void showCurrentAlignment(int len, int[][]set, String title){
//      BiojavaJmol jmol = new BiojavaJmol();
//      jmol.setTitle(title);
//
//      Chain c1 = new ChainImpl();
//      c1.setName("A");
//
//      Chain c2 = new ChainImpl();
//      c2.setName("B");
//      for(int i = 0; i < len; i ++)    {
//         Atom a1 = cod1[set[0][i]];
//         Atom a2 = cod2[set[1][i]];
//
//         Group g1 = a1.getParent();
//         Group g2 = a2.getParent();
//
//         try {
//            Group n1 = new AminoAcidImpl();
//            n1.setPDBCode(g1.getPDBCode());
//            n1.setPDBName(g1.getPDBName());
//            n1.addAtom(a1);
//
//            Group n2 = new AminoAcidImpl();
//            n2.setPDBCode(g2.getPDBCode());
//            n2.setPDBName(g2.getPDBName());
//            n2.addAtom(a2);
//
//
//            c1.addGroup(n1);
//            c2.addGroup(n2);
//         } catch (Exception e){
//            //
//         }
//
//      }
//
//      Structure s = new StructureImpl();
//      s.setPDBCode(title);
//      List<Chain> model1 = new ArrayList<Chain>();
//      model1.add(c1);
//      List<Chain> model2 = new ArrayList<Chain>();
//      model2.add(c2);
//      s.addModel(model1);
//      s.addModel(model2);
//      s.setNmr(true);
//      jmol.setStructure(s);
//      jmol.evalString("select *; backbone 0.4; wireframe off; spacefill off; " +
//      "select not protein and not solvent; spacefill on;");
//      jmol.evalString("select */1 ; color red; model 1; ");
//
//      // now show both models again.
//      jmol.evalString("model 0;");
//
//   }


	/** run the optimization
	 *
	 * @param maxi maximum nr. of iterations
	 * @throws StructureException
	 */
	public void runOptimization(int maxi) throws StructureException{
		superimposeBySet();
		if ( debug)
			System.err.println("   initial rmsd " + rmsd);

//      if (showAlig)
//         showCurrentAlignment(equLen, equSet, "after initial superimposeBySet Len:" +equLen + " rmsd:" +rmsd);

		maxKeepStep = 4;
		keepStep = 0;

		optimize(maxi);
	}

	/**
	 * refer CE, similarity = Dc - dij, Dc is increased by 0.5 each cycle,
	 * optimization continues until either
	 * i)alignment length is less than 95% of alignment length before optimization
	 * ii)rmsd is less than 110% of rmsd at the cycle when condition i) was first satisfied
	 */
	private void setParameters()
	{
		Dc = 3.0; //Dc = 2.0
		increase = 0.5;
		stopLenPer = 0.95;
		stopRmsdPer = 1.1;
		stopRmsd = -1.0;
		rmsdCut = 3.0;
		gapIni = 5.0;
		gapExt = 0.5;
	}


	/**
	 * superimpose two structures according to the equivalent residues
	 */
	private void superimposeBySet ()
	throws StructureException
	{

		//extract the coordinations of equivalent residues
		Atom[] tmp1 = new Atom[equLen];
		Atom[] tmp2 = new Atom[equLen];
		int     i,  r1, r2;


		for(i = 0; i < equLen; i ++)    {
			r1 = equSet[0][i];
			r2 = equSet[1][i];

			tmp1[i] =       cod1[ r1 ];
			tmp2[i] = (Atom)cod2[ r2 ].clone(); // have to be cloned!
			//tmp2[i] = cod2[ r2 ];


			/*try {
				System.out.println("before superimpos: " + equSet[0][i]+"-"+ equSet[1][i]+ " dist:" + Calc.getDistance(tmp1[i], cod2[equSet[1][i]]));
			} catch (Exception e){
				e.printStackTrace();
			}*/
		}

		//superimpose the equivalent residues
		Matrix4d trans = SuperPositions.superpose(Calc.atomsToPoints(tmp1),
				Calc.atomsToPoints(tmp2));

		Calc.transform(tmp2, trans);

		// weird, why does it take the RMSD before the rotation?
		// the rmsd is only for the subset contained in the tmp arrays.
		rmsd = Calc.rmsd(tmp1,tmp2);

		//System.err.println("rmsd after superimpose by set: " + rmsd);

		//transform structure 2 according to the superimposition of the equivalent residues
		Calc.transform(cod2, trans);

//      for(i = 0; i < equLen; i ++)    {
//         try {
//            System.err.println("after superimpos: " + equSet[0][i]+"-"+ equSet[1][i]+ " dist:" + Calc.getDistance(tmp1[i], cod2[equSet[1][i]]));
//         } catch (Exception e){
//            e.printStackTrace();
//         }
//      }

	}


	private void optimize(int maxi) throws StructureException
	{
		long optStart = System.currentTimeMillis();
		if ( debug)
			System.out.println("Optimizing up to " + maxi + " iterations.. ");
		boolean ifstop = true;;
		int     i, alnLen;
		alnLen = 0;

		int[][]     alnList =  new int[2][maxLen];
		for(i = 0; i < maxi; i ++)      {

			//if ( debug){
			//   System.out.println("optimize iteration: " + i);
			//}

			calMatrix();

			FCAlignHelper aln = new FCAlignHelper(sij,pro1Len,pro2Len,gapIni, gapExt);

			//ALIGN0 *aln = new ALIGN0(sij, pro1Len, pro2Len, gapIni, gapExt);
			alnLen = aln.getAlignPos(alnList);
			if(alnLen < 3)  ifstop = true; //very rare, mark by Y.Y on 5/1/03
			else    ifstop = defineEquPos(alnLen, alnList);

			if(ifstop)      break;
			Dc += increase;

//         if (showAlig)
//            if ( i == 0 )
//               showCurrentAlignment(alnLen, alnList,  "optimizing alignment - after " + i + " iterations alnLen:" + alnLen + " rmsd " + rmsd);
		}

		if  (debug){
			if(i < maxi)    {
				System.out.println(String.format("   optimize converged at %d iterations\n", i));
			}
			else    System.out.println("   optimize stop without convergence\n");
			System.out.println("optimization time: " + (System.currentTimeMillis() - optStart) + " ms.");
		}

//      if (showAlig)
//         showCurrentAlignment(alnLen, alnList,  "optimizing alignment - after " + i + " iterations alnLen:" + alnLen + " rmsd " + rmsd);
	}

	//--------------------------------------------------------------------------------------------------------
	//the definition of matrix between residues:
	//              Sij = Dc^2 - Dij^2 if Dij <= Dc
	//                    0            else
	//--------------------------------------------------------------------------------------------------------
	private void calMatrix() throws StructureException
	{
		int     i, j;
		double  dis;
		for(i = 0; i < pro1Len; i ++)   {
			for(j = 0; j < pro2Len; j ++)   {
				dis = Calc.getDistance(cod1[i],cod2[j]);

				if(dis < Dc) {
					sij[i][j] = Dc - dis;
				}
				else  {
					sij[i][j] = 0;
				}
			}
		}
	}

	/**
	 * the equivalent residues: residues where Dij &lt;= Dc and i,j is an aligned pair
	 * use the previous superimposing
	 */
	private boolean defineEquPos(int alnLen, int[][] alnList)
	throws StructureException
	{
		int     i, r1, r2;
		int     equLenOld = equLen;
		int[][]    equSetOld = new int[2][equLenOld];
		for(i = 0; i < equLen; i ++)    {
			equSetOld[0][i] = equSet[0][i];
			equSetOld[1][i] = equSet[1][i];
		}
		double  rmsdOld = rmsd;
		double  dis;
		equLen = 0;
		//if (debug)
		//   System.out.println(String.format(" OPT: Dc %f, equLenOld %d, rmsdOld %f, alnLen %d", Dc, equLenOld, rmsdOld, alnLen));
		for(i = 0; i < alnLen; i ++)    {
			r1 = alnList[0][i];
			r2 = alnList[1][i];
			dis = Calc.getDistance(cod1[r1],cod2[r2]);
			if(dis <= Dc)   {
				//System.out.println(r1 + "-"  + r2 + " d:" + dis);
				equSet[0][equLen] = r1;
				equSet[1][equLen] = r2;
				equLen ++;
			}
		}

		superimposeBySet();

		//if (debug)
		//   System.out.println(String.format(" OPT: new equLen %d rmsd %f", equLen, rmsd));

		boolean     ifstop = false;

//      if (debug) {
//         System.out.print(" OPT: rmsd diff: " + Math.abs(rmsd - rmsdOld) + " equLens: " + equLenOld + ":"+ equLen);
//         if ( Math.abs(rmsd - rmsdOld) < 1e-10)
//            System.out.println(" NO DIFF!");
//         else
//            System.out.println(" DIFF!");
//      }

		if((Math.abs(rmsd - rmsdOld) < 1e-10 ) && (equLenOld == equLen)) keepStep ++;
		else    keepStep = 0;

		if(keepStep > maxKeepStep)      {
			ifstop = true; //converge
		} //allowing up to maxKeepStep instead of 1 is essential for some special cases
		else if(stopRmsd < 0) {
			ifstop = false; //condition 1, continue
		}
		else if((rmsd <= stopRmsd * stopRmsdPer) || (rmsd < rmsdCut)) {
			ifstop = false; //condition 2, continue
		} //rmsdCut is adopted or not? to be tuned
		else    {
			ifstop = true; //get worse
		}


		if((stopRmsd <0) && (equLen >= stopLenPer * equLen0))    {
			//System.err.println("stopRmsd: " + stopRmsd + " setting to rmsd:" + rmsd);
			stopRmsd = rmsd; //condition 1
		}

		return ifstop;
	}


	public double optimizeResult(int[] optLen, int optLenpos, int[][] list)
	{
		optLen[optLenpos] = equLen;

		for(int i = 0; i < equLen; i ++)        {
			list[0][i] = equSet[0][i];
			list[1][i] = equSet[1][i];
		}
		return rmsd;
	}

}
