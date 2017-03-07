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

package org.biojava.nbio.structure.align;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.model.AFP;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.geometry.Matrices;
import org.biojava.nbio.structure.geometry.SuperPositions;
import org.biojava.nbio.structure.jama.Matrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

//import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;

public class AFPTwister {
	private final static Logger logger = LoggerFactory
			.getLogger(AFPTwister.class);

	// private static final boolean showAlignmentSteps = false;

	/**
	 * calculate the total rmsd of the blocks output a merged pdb file for both
	 * proteins protein 1, in chain A protein 2 is twisted according to the
	 * twists detected, in chain B
	 * 
	 * @return twisted Groups
	 */
	public static Group[] twistPDB(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
			throws StructureException {
		// --------------------------------------------------------

		if (afpChain.isShortAlign())
			return new Group[0];

		List<AFP> afpSet = afpChain.getAfpSet();

		int blockNum = afpChain.getBlockNum();

		int i, b2, e2;

		// superimposing according to the initial AFP-chaining
		Atom[] origCA = StructureTools.cloneAtomArray(ca2);
		Atom[] iniTwistPdb = StructureTools.cloneAtomArray(ca2);

		int[] blockResSize = afpChain.getBlockResSize();

		int[][][] blockResList = afpChain.getBlockResList();

		int[] afpChainList = afpChain.getAfpChainList();
		int[] block2Afp = afpChain.getBlock2Afp();
		int[] blockSize = afpChain.getBlockSize();
		int[] focusAfpList = afpChain.getFocusAfpList();

		if (focusAfpList == null) {
			focusAfpList = new int[afpChain.getMinLen()];
			afpChain.setFocusAfpList(focusAfpList);
		}

		int focusAfpn = 0;
		e2 = 0;
		b2 = 0;

		logger.debug("blockNUm at twister: ", blockNum);

		for (int bk = 0; bk < blockNum; bk++) {

			// THIS IS TRANSFORMING THE ORIGINAL ca2 COORDINATES, NO CLONING...
			// copies the atoms over to iniTwistPdb later on in modifyCod

			transformOrigPDB(blockResSize[bk], blockResList[bk][0],
					blockResList[bk][1], ca1, ca2, null, -1);

			// transform pro2 according to comparison of pro1 and pro2 at give
			// residues
			if (bk > 0) {
				b2 = e2;
			}

			logger.debug("b2 is {} before  modifyCon", b2);

			if (bk < (blockNum - 1)) {

				// bend at the middle of two consecutive AFPs
				int afpPos = afpChainList[block2Afp[bk] + blockSize[bk] - 1];
				AFP a1 = afpSet.get(afpPos);
				e2 = a1.getP2();

				int afpPos2 = afpChainList[block2Afp[bk + 1]];
				AFP a2 = afpSet.get(afpPos2);
				e2 = (a2.getP2() - e2) / 2 + e2;
				logger.debug("e2 : {}", e2);
			} else {
				// last one is until the end...
				e2 = ca2.length;
			}

			// this copies the coordinates over into iniTwistPdb
			cloneAtomRange(iniTwistPdb, ca2, b2, e2);

			// bound[bk] = e2;
			for (i = 0; i < blockSize[bk]; i++) {
				focusAfpList[focusAfpn] = afpChainList[block2Afp[bk] + i];
				focusAfpn++;
			}
		}

		int focusResn = afp2Res(afpChain, focusAfpn, focusAfpList, 0);

		afpChain.setTotalLenIni(focusResn);

		logger.debug(String.format("calrmsdini for %d residues", focusResn));

		double totalRmsdIni = calCaRmsd(ca1, iniTwistPdb, focusResn,
				afpChain.getFocusRes1(), afpChain.getFocusRes2());

		afpChain.setTotalRmsdIni(totalRmsdIni);
		logger.debug("got iniRMSD: {}", totalRmsdIni);
		if (totalRmsdIni == 5.76611141613097) {
			logger.debug("{}", afpChain.getAlnseq1());
			logger.debug("{}", afpChain.getAlnsymb());
			logger.debug("{}", afpChain.getAlnseq2());
		}

		afpChain.setFocusAfpList(focusAfpList);
		afpChain.setBlock2Afp(block2Afp);
		afpChain.setAfpChainList(afpChainList);

		return twistOptimized(afpChain, ca1, origCA);
	}

	/**
	 * superimposing according to the optimized alignment
	 *
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return Group array twisted.
	 * @throws StructureException
	 */
	public static Group[] twistOptimized(AFPChain afpChain, Atom[] ca1,
			Atom[] ca2) throws StructureException {

		Atom[] optTwistPdb = new Atom[ca2.length];

		int gPos = -1;
		for (Atom a : ca2) {
			gPos++;
			optTwistPdb[gPos] = a;
		}

		int blockNum = afpChain.getBlockNum();

		int b2 = 0;
		int e2 = 0;
		int focusResn = 0;
		int[] focusRes1 = afpChain.getFocusRes1();
		int[] focusRes2 = afpChain.getFocusRes2();

		if (focusRes1 == null) {
			focusRes1 = new int[afpChain.getCa1Length()];
			afpChain.setFocusRes1(focusRes1);
		}
		if (focusRes2 == null) {
			focusRes2 = new int[afpChain.getCa2Length()];
			afpChain.setFocusRes2(focusRes2);
		}

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();

		for (int bk = 0; bk < blockNum; bk++) {
			// THIS IS TRANSFORMING THE ORIGINAL ca2 COORDINATES, NO CLONING...
			// copies the atoms over to iniTwistPdb later on in modifyCod
			transformOrigPDB(optLen[bk], optAln[bk][0], optAln[bk][1], ca1,
					ca2, afpChain, bk);

			// transform pro2 according to comparison of pro1 and pro2 at give
			// residues
			if (bk > 0) {
				b2 = e2;
			}
			if (bk < blockNum - 1) { // bend at the middle of two consecutive
										// blocks
				e2 = optAln[bk][1][optLen[bk] - 1];
				e2 = (optAln[bk + 1][1][0] - e2) / 2 + e2;
			} else {
				e2 = ca2.length;
			}
			cloneAtomRange(optTwistPdb, ca2, b2, e2);
			for (int i = 0; i < optLen[bk]; i++) {
				focusRes1[focusResn] = optAln[bk][0][i];
				focusRes2[focusResn] = optAln[bk][1][i];
				focusResn++;
			}
		}
		int totalLenOpt = focusResn;
		logger.debug("calrmsdopt for {} residues", focusResn);
		double totalRmsdOpt = calCaRmsd(ca1, optTwistPdb, focusResn, focusRes1,
				focusRes2);
		logger.debug("got opt RMSD: {}", totalRmsdOpt);
		int optLength = afpChain.getOptLength();

		if (totalLenOpt != optLength) {
			logger.warn("Final alignment length is different {} {}",
					totalLenOpt, optLength);
		}
		logger.debug("final alignment length {}, rmsd {}", focusResn,
				totalRmsdOpt);

		afpChain.setTotalLenOpt(totalLenOpt);
		afpChain.setTotalRmsdOpt(totalRmsdOpt);

		return StructureTools.cloneGroups(optTwistPdb);

	}

	/**
	 * transform the coordinates in the ca2 according to the superimposing of
	 * the given position pairs. No Cloning, transforms input atoms.
	 */
	// orig name: transPdb
	private static void transformOrigPDB(int n, int[] res1, int[] res2,
			Atom[] ca1, Atom[] ca2, AFPChain afpChain, int blockNr)
			throws StructureException {
		logger.debug(
				"transforming original coordinates {} len1: {} res1: {} len2: {} res2: {}",
				n, ca1.length, res1.length, ca2.length, res2.length);

		Atom[] cod1 = getAtoms(ca1, res1, n, false);
		Atom[] cod2 = getAtoms(ca2, res2, n, false);

		// double *cod1 = pro1->Cod4Res(n, res1);
		// double *cod2 = pro2->Cod4Res(n, res2);

		Matrix4d transform = SuperPositions.superpose(Calc.atomsToPoints(cod1),
				Calc.atomsToPoints(cod2));

		Matrix r = Matrices.getRotationJAMA(transform);
		Atom t = Calc.getTranslationVector(transform);

		logger.debug("transPdb: transforming orig coordinates with matrix: {}",
				r);

		if (afpChain != null) {
			Matrix[] ms = afpChain.getBlockRotationMatrix();
			if (ms == null)
				ms = new Matrix[afpChain.getBlockNum()];

			ms[blockNr] = r;

			Atom[] shifts = afpChain.getBlockShiftVector();
			if (shifts == null)
				shifts = new Atom[afpChain.getBlockNum()];
			shifts[blockNr] = t;

			afpChain.setBlockRotationMatrix(ms);
			afpChain.setBlockShiftVector(shifts);
		}

		for (Atom a : ca2)
			Calc.transform(a.getGroup(), transform);

	}

	// like Cod4Res
	// most likely the clone flag is not needed
	private static Atom[] getAtoms(Atom[] ca, int[] positions, int length,
			boolean clone) {

		List<Atom> atoms = new ArrayList<Atom>();
		for (int i = 0; i < length; i++) {
			int p = positions[i];
			Atom a;
			if (clone) {
				a = (Atom) ca[p].clone();
				a.setGroup((Group) ca[p].getGroup().clone());
			} else {
				a = ca[p];
			}
			atoms.add(a);
		}
		return atoms.toArray(new Atom[atoms.size()]);
	}

	/**
	 * Clones and copies the Atoms from p2 into p1 range is between r1 and r2
	 *
	 * @param p1
	 * @param p2
	 * @param r1
	 * @param r2
	 */
	// orig name: modifyCod
	private static void cloneAtomRange(Atom[] p1, Atom[] p2, int r1, int r2)
			throws StructureException {

		logger.debug("modifyCod from: {} to: {}", r1, r2);

		// special clone method, can;t use StructureTools.cloneCAArray, since we
		// access the data
		// slightly differently here.
		List<Chain> model = new ArrayList<Chain>();
		for (int i = r1; i < r2; i++) {

			Group g = p2[i].getGroup();
			Group newG = (Group) g.clone();

			p1[i] = newG.getAtom(p2[i].getName());
			Chain parentC = g.getChain();

			Chain newChain = null;

			for (Chain c : model) {
				if (c.getName().equals(parentC.getName())) {
					newChain = c;
					break;
				}
			}

			if (newChain == null) {
				newChain = new ChainImpl();
				newChain.setName(parentC.getName());
				model.add(newChain);
			}

			newChain.addGroup(newG);

		} // modify caCod

	}

	/**
	 * Return the rmsd of the CAs between the input pro and this protein, at
	 * given positions. quite similar to transPdb but while that one transforms
	 * the whole ca2, this one only works on the res1 and res2 positions.
	 *
	 * Modifies the coordinates in the second set of Atoms (pro).
	 *
	 * @return rmsd of CAs
	 */
	private static double calCaRmsd(Atom[] ca1, Atom[] pro, int resn,
			int[] res1, int[] res2) throws StructureException {

		Atom[] cod1 = getAtoms(ca1, res1, resn, false);
		Atom[] cod2 = getAtoms(pro, res2, resn, false);

		if (cod1.length == 0 || cod2.length == 0) {
			logger.info("length of atoms  == 0!");
			return 99;
		}

		Matrix4d transform = SuperPositions.superpose(Calc.atomsToPoints(cod1),
				Calc.atomsToPoints(cod2));

		for (Atom a : cod2)
			Calc.transform(a.getGroup(), transform);

		return Calc.rmsd(cod1, cod2);
	}

	/**
	 * Set the list of equivalent residues in the two proteins given a list of
	 * AFPs
	 *
	 * WARNING: changes the values for FocusRes1, focusRes2 and FocusResn in
	 * afpChain!
	 *
	 * @param afpChain
	 *            the AFPChain to store resuts
	 * @param afpn
	 *            nr of afp
	 * @param afpPositions
	 * @param listStart
	 * @return nr of eq residues
	 */

	public static int afp2Res(AFPChain afpChain, int afpn, int[] afpPositions,
			int listStart) {
		int[] res1 = afpChain.getFocusRes1();
		int[] res2 = afpChain.getFocusRes2();
		int minLen = afpChain.getMinLen();

		int n = 0;

		List<AFP> afpSet = afpChain.getAfpSet();

		for (int i = listStart; i < listStart + afpn; i++) {
			int a = afpPositions[i];
			for (int j = 0; j < afpSet.get(a).getFragLen(); j++) {
				if (n >= minLen) {
					throw new RuntimeException(
							"Error: too many residues in AFPChainer.afp2Res!");
				}
				res1[n] = afpSet.get(a).getP1() + j;
				res2[n] = afpSet.get(a).getP2() + j;
				n++;
			}
		}

		afpChain.setFocusRes1(res1);
		afpChain.setFocusRes2(res2);
		afpChain.setFocusResn(n);

		if (n == 0) {
			logger.warn("n=0!!! + " + afpn + " " + listStart + " "
					+ afpPositions.length);
		}
		return n;
	}

}
