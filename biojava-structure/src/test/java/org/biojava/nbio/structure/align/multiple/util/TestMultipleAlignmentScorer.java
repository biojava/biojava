package org.biojava.nbio.structure.align.multiple.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.ReferenceSuperimposer;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the correctness of various Score calculations for 
 * {@link MultipleAlignment}s.
 * <p>
 * Currently tested:
 * <ul><li>Reference and total RMSD
 * <li>Reference and average TM-Score
 * <li>MC-Score
 * </ul>
 * 
 * @author Aleix Lafita
 *
 */
public class TestMultipleAlignmentScorer {

	@Test
	public void testRefRMSD() throws Exception{
		
		//Identity Test: RefRMSD has to be equal to 0.0
		MultipleAlignment identMSA = identityMSTA();
		double refRMSD1 = MultipleAlignmentScorer.getRefRMSD(identMSA, 0);
		assertEquals(0.0, refRMSD1, 0.000001);
		
		//Simple Test: RefRMSD has to be equal to sqrt(2.5)
		MultipleAlignment simpleMSA = simpleMSTA();
		double refRMSD2 = MultipleAlignmentScorer.getRefRMSD(simpleMSA, 0);
		assertEquals(1.5811388, refRMSD2, 0.000001);
		
		//Gapped Test: RefRMSD has to be equal to 1.0
		MultipleAlignment gappedMSA = gappedMSTA();
		double refRMSD3 = MultipleAlignmentScorer.getRefRMSD(gappedMSA, 0);
		assertEquals(1.0, refRMSD3, 0.000001);
	}
	
	@Test
	public void testRMSD() throws Exception{
		
		//Identity Test: RMSD has to be equal to 0.0
		MultipleAlignment identMSA = identityMSTA();
		double RMSD1 = MultipleAlignmentScorer.getRMSD(identMSA);
		assertEquals(0.0, RMSD1, 0.000001);
		
		//Simple Test: RMSD has to be equal to 1.3228756
		MultipleAlignment simpleMSA = simpleMSTA();
		double RMSD2 = MultipleAlignmentScorer.getRMSD(simpleMSA);
		assertEquals(1.3228756, RMSD2, 0.000001);
		
		//Gapped Test: RMSD has to be equal to 1.0
		MultipleAlignment gappedMSA = gappedMSTA();
		double RMSD3 = MultipleAlignmentScorer.getRMSD(gappedMSA);
		assertEquals(0.8340576, RMSD3, 0.000001);
	}
	
	@Test
	public void testRefTMScore() throws Exception{
		
		//Identity Test: RefTM-Score has to be equal to 1.0
		MultipleAlignment identMSA = identityMSTA();
		double refTM1 = MultipleAlignmentScorer.getRefTMScore(identMSA, 0);
		assertEquals(1.0, refTM1, 0.000001);
		
		//Simple Test: RefTM-Score has to be equal to 0.6831032
		MultipleAlignment simpleMSA = simpleMSTA();
		double refTM2 = MultipleAlignmentScorer.getRefTMScore(simpleMSA, 0);
		assertEquals(0.6831032, refTM2, 0.000001);
		
		//Simple Test: RefTM-Score has to be equal to 0.6831032
		MultipleAlignment gappedMSA = gappedMSTA();
		double refTM3 = MultipleAlignmentScorer.getRefTMScore(gappedMSA, 0);
		assertEquals(0.2672780, refTM3, 0.000001);
	}
	
	@Test
	public void testAvgTMScore() throws Exception{
		
		//Identity Test: AvgTM-Score has to be equal to 1.0
		MultipleAlignment identMSA = identityMSTA();
		double avgTM1 = MultipleAlignmentScorer.getAvgTMScore(identMSA);
		assertEquals(1.0, avgTM1, 0.000001);
		
		//Simple Test: AvgTM-Score has to be equal to 0.7261305
		MultipleAlignment simpleMSA = simpleMSTA();
		double avgTM2 = MultipleAlignmentScorer.getAvgTMScore(simpleMSA);
		assertEquals(0.7261305, avgTM2, 0.000001);
		
		//Simple Test: AvgTM-Score has to be equal to 0.2831998
		MultipleAlignment gappedMSA = gappedMSTA();
		double avgTM3 = MultipleAlignmentScorer.getAvgTMScore(gappedMSA);
		assertEquals(0.2831998, avgTM3, 0.000001);
	}
	
	@Test
	public void testMCScore() throws Exception {
		
		//Identity Test: MultipleMC-Score has to be equal to 576.21
		MultipleAlignment identMSA = identityMSTA();
		double MC1 = MultipleAlignmentScorer.getMCScore(identMSA, 10, 5, 10);
		assertEquals(574.98, MC1, 0.01);
		
		//Simple Test: MultipleMC-Score has to be equal to 2351.93
		MultipleAlignment simpleMSA = simpleMSTA();
		double MC2 = MultipleAlignmentScorer.getMCScore(simpleMSA, 10, 5, 10);
		assertEquals(2197.20, MC2, 0.01);
		
		//Simple Test: MultipleMC-Score has to be equal to 489.74
		MultipleAlignment gappedMSA = gappedMSTA();
		double MC3 = MultipleAlignmentScorer.getMCScore(gappedMSA, 10, 5, 10);
		assertEquals(478.52, MC3, 0.01);
	}
	
	/**
	 * Generates an identity MultipleAlignment: 3 structures with 
	 * the same Atoms and perfectly aligned, so that TM-score = 1 
	 * and RMSD = 0.
	 * @return MultipleAlignment identity
	 * @throws StructureException 
	 */
	private MultipleAlignment identityMSTA() throws StructureException {
		
		//Generate the identical Atom arrays
		List<Atom[]> atomArrays = new ArrayList<Atom[]>(20);
		for (int i=0; i<3; i++) atomArrays.add(makeDummyCA(20));
		
		//Generate the identity alignment (1-1-1,2-2-2,etc)
		List<List<Integer>> alnRes = new ArrayList<List<Integer>>(3);
		for (int str=0; str<3; str++){
			List<Integer> chain = new ArrayList<Integer>(20);
			for (int res=0; res<20; res++) chain.add(res);
			alnRes.add(chain);
		}
		
		//MultipleAlignment generation
		MultipleAlignment msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(atomArrays);
		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(alnRes);
		
		//Superimpose the alignment (which should give the identity matrices)
		ReferenceSuperimposer imposer = new ReferenceSuperimposer();
		imposer.superimpose(msa);
		
		return msa;
	}
	
	/**
	 * Generates a simple MultipleAlignment: 3 structures with the same 
	 * Atoms but incorreclty aligned (offset of 1 position) without gaps.
	 * 
	 * @return MultipleAlignment simple MSTA
	 * @throws StructureException 
	 */
	private MultipleAlignment simpleMSTA() throws StructureException{
		
		//Generate three identical Atom arrays
		List<Atom[]> atomArrays = new ArrayList<Atom[]>(52);
		for (int i=0; i<3; i++) atomArrays.add(makeDummyCA(52));
		
		//Generate the incorrect alignment (0-1-2,1-2-3,etc)
		List<List<Integer>> alnRes = new ArrayList<List<Integer>>(3);
		for (int str=0; str<3; str++){
			List<Integer> chain = new ArrayList<Integer>(50);
			for (int res=0; res<50; res++) chain.add(res+str);
			alnRes.add(chain);
		}
		
		//MultipleAlignment generation
		MultipleAlignment msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(atomArrays);
		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(alnRes);
		
		//We want the identity transfromations to maintain the missalignment
		Matrix4d ident = new Matrix4d();
		ident.setIdentity();
		bs.setTransformations(Arrays.asList(ident,ident,ident));
		
		return msa;
	}
	
	/**
	 * Generates a simple MultipleAlignment: 3 structures with the 
	 * same Atoms but incorreclty aligned with gaps.
	 * 
	 * @return MultipleAlignment gapped MSTA
	 * @throws StructureException 
	 */
	private MultipleAlignment gappedMSTA() throws StructureException{
		
		//Generate three identical Atom arrays
		List<Atom[]> atomArrays = new ArrayList<Atom[]>(30);
		for (int i=0; i<3; i++) atomArrays.add(makeDummyCA(30));
		
		//Generate alignment with nulls and some missalignments
		List<List<Integer>> alnRes = new ArrayList<List<Integer>>(3);
		List<Integer> chain1 = Arrays.asList(
				1, 2, 	 3, 5, 8, 10, 12, 	 15, 17,  19, 22, null, 24, 27);
		List<Integer> chain2 = Arrays.asList(
				1, null, 3, 6, 9, 11, 12,   15, null, 18, 22, 24,   26, 28);
		List<Integer> chain3 = Arrays.asList(
				1, 2,    4, 7, 9, 10, null, 15, null, 17, 22, 24,   26, 28);
		
		alnRes.add(chain1);
		alnRes.add(chain2);
		alnRes.add(chain3);
		
		//MultipleAlignment generation
		MultipleAlignment msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(atomArrays);
		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(alnRes);
		
		//We want the identity transfromations to mantain the missalignments
		Matrix4d ident = new Matrix4d();
		ident.setIdentity();
		bs.setTransformations(Arrays.asList(ident,ident,ident));
		
		return msa;
	}
	
	/**
	 * Makes dummy CA atoms at 1A intervals. Only the x coordinate increments
	 * by one at each consecutive Atom.
	 */
	private Atom[] makeDummyCA(int len) {
		Atom[] ca1;
		Chain chain1 = new ChainImpl();
		ca1 = new Atom[len];
		for(int i=0;i<len;i++) {
			ca1[i] = new AtomImpl();
			ca1[i].setName("CA");
			ca1[i].setCoords(new double[] { i, 0, 0 });
			Group aa = new AminoAcidImpl();
			aa.setPDBName("GLY");
			aa.setResidueNumber( ResidueNumber.fromString(i+""));
			aa.addAtom(ca1[i]);
			chain1.addGroup(aa);
		}
		return ca1;
	}
}
