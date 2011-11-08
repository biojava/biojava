/**
 * 
 */
package org.biojava.bio.structure.align.ce;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;

/**
 * @author Spencer Bliven
 *
 */
public class OptimalCECPMainTest extends TestCase {

	AtomCache cache = new AtomCache();
	
	/* (non-Javadoc)
	 * @see junit.framework.TestCase#setUp()
	 */
	protected void setUp() throws Exception {
		super.setUp();
	}
	
	/**
	 * Basic test that alignPermuted(..., 0) is equivalent to a normal CE alignment.
	 * 
	 * Also checks that {@link AFPChain#equals(Object)} is working the way we expect.
	 * @throws IOException
	 * @throws StructureException
	 */
	public void testUnpermuted() throws IOException, StructureException {
		String name1, name2;
		
		//small case
		name1 = "d1qdmA1";
		name2 = "d1nklA_";
		
		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		// Calculate all alignments initially
		OptimalCECPMain cecp = new OptimalCECPMain();
		Atom[] ca2clone = cache.getAtoms(name2);
		AFPChain cp0 = cecp.alignPermuted(ca1, ca2clone, (OptimalCECPParameters) cecp.getParameters(), 0);

		CeMain ce = new CeMain();
		AFPChain nocp = ce.align(ca1,ca2);
		
		assertEquals(nocp,cp0);
	}
	
	/**
	 * Very basic test of {@link OptimalCECPMain#permuteOptAln(AFPChain, int)}
	 * 
	 * It should do nothing on unpermuted alignments.
	 * @throws NoSuchMethodException 
	 * @throws SecurityException 
	 * @throws StructureException 
	 * @throws IOException 
	 * @throws InvocationTargetException 
	 * @throws IllegalAccessException 
	 * @throws IllegalArgumentException 
	 */
	public void testPermuteOptAlnUnpermuted() throws SecurityException, NoSuchMethodException, StructureException, IOException, IllegalArgumentException, IllegalAccessException, InvocationTargetException {
		//test private member using reflection
		Method permuteOptAln = OptimalCECPMain.class.getDeclaredMethod(
				"permuteOptAln", AFPChain.class, int.class);
		permuteOptAln.setAccessible(true);

		String name1, name2;
		name1 = "d1qdmA1";
		name2 = "d1nklA_";
		
		CeMain ce = (CeCPMain) StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);

		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);
		
		AFPChain afpChain = ce.align(ca1, ca2);
		AFPChain afpChain2 = (AFPChain) afpChain.clone();
		
		permuteOptAln.invoke(null, afpChain2, 0);
		
		assertEquals("Permuting by 0 changed the alignment!",afpChain, afpChain2);
	}
	
	/**
	 * Checks that individual alignments performed by alignOptimal are consistent
	 * with the alignments returned by individual calls to alignPermuted.
	 * 
	 * This addresses a bug involving multiple calls to align() on the same
	 * CE instance.
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	public void testOptimalAlignmentConsistency() throws IOException, StructureException {
		String name1, name2;
		OptimalCECPMain ce;
		AFPChain afpChain;
		int[] cps; //CP points to check for consistency
		
		//small case
		name1 = "d1qdmA1";
		name2 = "d1nklA_";
		cps = new int[] {0,1,2,41,5,38};
		
		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		// Calculate all alignments initially
		ce = new OptimalCECPMain();
		AFPChain[] alignments = new AFPChain[ca2.length];
		ce.alignOptimal(ca1, ca2, (OptimalCECPParameters) ce.getParameters(), alignments);
		
		for(int cp : cps) {
			// fresh instance to avoid contamination
			ce = new OptimalCECPMain();
			
			// new copy of ca2, since alignPermuted has side effects
			Atom[] ca2clone = cache.getAtoms(name2);
			afpChain = ce.alignPermuted(ca1, ca2clone, (OptimalCECPParameters) ce.getParameters(), cp);

			assertEquals("Alignment "+cp+" differs.",afpChain, alignments[cp]);
		}
		
	}
	
	/**
	 * Tests private {@link OptimalCECPMain#permuteArray(Object[], int)}
	 * @throws Exception
	 */
	public void testPermuteArray() throws Exception {
		//test private member using reflection
		Method permuteArray = OptimalCECPMain.class.getDeclaredMethod(
				"permuteArray", Object[].class, int.class);
		permuteArray.setAccessible(true);

		String[] arr0 = new String[] {"A","B","C","D","E","F"};
		String[] arr1 = new String[] {"B","C","D","E","F","A"};
		String[] arr5 = new String[] {"F","A","B","C","D","E"};

		String[] arrP;

		arrP = Arrays.copyOf(arr0, arr0.length);
		assertTrue("Shallow equals!",Arrays.deepEquals(arr0, arrP));
		permuteArray.invoke(null, arrP, 0);
		assertTrue(String.format("Permuting by 0 gave %s%s%s%s%s%s",(Object[])arrP),
				Arrays.deepEquals(arr0, arrP));
		
		arrP = Arrays.copyOf(arr0, arr0.length);
		permuteArray.invoke(null, arrP, 1);
		assertTrue(String.format("Permuting by 1 gave %s%s%s%s%s%s",(Object[])arrP),
				Arrays.deepEquals(arr1, arrP));

		arrP = Arrays.copyOf(arr0, arr0.length);
		permuteArray.invoke(null, arrP, 5);
		assertTrue(String.format("Permuting by 7 gave %s%s%s%s%s%s",(Object[])arrP),
				Arrays.deepEquals(arr5, arrP));

		arrP = Arrays.copyOf(arr0, arr0.length);
		permuteArray.invoke(null, arrP, -1);
		assertTrue(String.format("Permuting by -1 gave %s%s%s%s%s%s",(Object[])arrP),
				Arrays.deepEquals(arr5, arrP));

		try {
			arrP = Arrays.copyOf(arr0, arr0.length);
			permuteArray.invoke(null, arrP, 6);
			fail("Illegal index. Should throw exception.");
		} catch( InvocationTargetException e) {
			if( ! (e.getCause() instanceof ArrayIndexOutOfBoundsException)) {
				throw e;
			}
		}
	}

	/**
	 * Tests private {@link OptimalCECPMain#permuteOptAln(AFPChain, int)}
	 */
	public void testPermuteOptAln() throws Exception {
		//test private member using reflection
		Method permuteOptAln = OptimalCECPMain.class.getDeclaredMethod(
				"permuteOptAln", AFPChain.class, int.class);
		permuteOptAln.setAccessible(true);
		
		// Two structures with nearly 100% sequence identity
		/* 
		 * Aligned (0-based index):
		 * 	3LB9.A	1HV1
		 * 	------	----
		 * 	0-62	122-184
		 * 	63		0
		 * 	65-181	1-117
		 * 
		 * unaligned:
		 * 	64		-
		 * 	-		118-121
		 * 
		 * PDB numbering:
		 * 	+2		+1
		 *
		 */
		String name1, name2;
		name1 = "3LB9.A";
		name2 = "1HV1";
		
		CeCPMain ce = (CeCPMain) StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);

		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);
		
		// Create permuted CA chain
		Method permuteArray = OptimalCECPMain.class.getDeclaredMethod(
				"permuteArray", Object[].class, int.class);
		permuteArray.setAccessible(true);
		Atom[] ca2p = StructureTools.cloneCAArray(ca2);
		permuteArray.invoke(null, ca2p, 63);
		
		AFPChain cpAlignment = ce.align(ca1, ca2);
		//System.out.println(cpAlignment.toCE(ca1, ca2));
		//printOptAln(cpAlignment);
		
		assertNotNull(cpAlignment);
		
		
		int[] optLen = cpAlignment.getOptLen();
		int[][][] optAln = cpAlignment.getOptAln();
		
		
		assertEquals("Wrong total length",181,cpAlignment.getOptLength());
		assertEquals("Wrong number of blocks",2, cpAlignment.getBlockNum());
		assertEquals("Wrong block 0 length",63,optLen[0]);
		assertEquals("Wrong block 1 length",118,optLen[1]);

		//just test some key positions in each block
		assertEquals("Wrong residue at start of block 0, protein 0",0,optAln[0][0][0]);
		assertEquals("Wrong residue at start of block 0, protein 1",122,optAln[0][1][0]);
		assertEquals("Wrong residue at end of block 0, protein 0",62,optAln[0][0][62]);
		assertEquals("Wrong residue at end of block 0, protein 1",184,optAln[0][1][62]);

		assertEquals("Wrong residue at start of block 1, protein 0",63,optAln[1][0][0]);
		assertEquals("Wrong residue at start of block 1, protein 1",0,optAln[1][1][0]);
		assertEquals("Wrong residue at pos 1 of block 1, protein 0",65,optAln[1][0][1]);
		assertEquals("Wrong residue at pos 1 of block 1, protein 1",1,optAln[1][1][1]);
		assertEquals("Wrong residue at pos 54 of block 1, protein 0",118,optAln[1][0][54]);
		assertEquals("Wrong residue at pos 54 of block 1, protein 1",54,optAln[1][1][54]);
		assertEquals("Wrong residue at pos 55 of block 1, protein 0",119,optAln[1][0][55]);
		assertEquals("Wrong residue at pos 55 of block 1, protein 1",55,optAln[1][1][55]);
		assertEquals("Wrong residue at end of block 1, protein 0",181,optAln[1][0][117]);
		assertEquals("Wrong residue at end of block 1, protein 1",117,optAln[1][1][117]);

		
		// permute! should align at 0,0
		//System.out.println("Permuting by 63 residues...");
		permuteOptAln.invoke(null, cpAlignment, 63);
		//System.out.println(cpAlignment.toCE(ca1, ca2p));
		//printOptAln(cpAlignment);
		
		optLen = cpAlignment.getOptLen();
		optAln = cpAlignment.getOptAln();
		
		assertEquals("Wrong total length",181,cpAlignment.getOptLength());
		assertEquals("Wrong number of blocks",2, cpAlignment.getBlockNum());
		assertEquals("Wrong block 0 length",63,optLen[0]);
		assertEquals("Wrong block 1 length",118,optLen[1]);		

		//just test some key positions in each block
		assertEquals("Wrong residue at start of block 0, protein 0",0,optAln[0][0][0]);
		assertEquals("Wrong residue at start of block 0, protein 1",0,optAln[0][1][0]);
		assertEquals("Wrong residue at end of block 0, protein 0",62,optAln[0][0][62]);
		assertEquals("Wrong residue at end of block 0, protein 1",62,optAln[0][1][62]);

		assertEquals("Wrong residue at start of block 1, protein 0",63,optAln[1][0][0]);
		assertEquals("Wrong residue at start of block 1, protein 1",63,optAln[1][1][0]);
		assertEquals("Wrong residue at pos 1 of block 1, protein 0",65,optAln[1][0][1]);
		assertEquals("Wrong residue at pos 1 of block 1, protein 1",64,optAln[1][1][1]);
		assertEquals("Wrong residue at end of block 1, protein 0",181,optAln[1][0][117]);
		assertEquals("Wrong residue at end of block 1, protein 1",180,optAln[1][1][117]);

		
		// undo permutation
		//System.out.println("Permuting by -63 residues...");
		permuteOptAln.invoke(null, cpAlignment, -63);
		//System.out.println(cpAlignment.toCE(ca1, ca2));
		//printOptAln(cpAlignment);
		
		optLen = cpAlignment.getOptLen();
		optAln = cpAlignment.getOptAln();
		
		assertEquals("Wrong total length",181,cpAlignment.getOptLength());
		assertEquals("Wrong number of blocks",2, cpAlignment.getBlockNum());
		assertEquals("Wrong block 0 length",63,optLen[0]);
		assertEquals("Wrong block 1 length",118,optLen[1]);

		//just test some key positions in each block
		assertEquals("Wrong residue at start of block 0, protein 0",0,optAln[0][0][0]);
		assertEquals("Wrong residue at start of block 0, protein 1",122,optAln[0][1][0]);
		assertEquals("Wrong residue at end of block 0, protein 0",62,optAln[0][0][62]);
		assertEquals("Wrong residue at end of block 0, protein 1",184,optAln[0][1][62]);

		assertEquals("Wrong residue at start of block 1, protein 0",63,optAln[1][0][0]);
		assertEquals("Wrong residue at start of block 1, protein 1",0,optAln[1][1][0]);
		assertEquals("Wrong residue at pos 1 of block 1, protein 0",65,optAln[1][0][1]);
		assertEquals("Wrong residue at pos 1 of block 1, protein 1",1,optAln[1][1][1]);
		assertEquals("Wrong residue at pos 54 of block 1, protein 0",118,optAln[1][0][54]);
		assertEquals("Wrong residue at pos 54 of block 1, protein 1",54,optAln[1][1][54]);
		assertEquals("Wrong residue at pos 55 of block 1, protein 0",119,optAln[1][0][55]);
		assertEquals("Wrong residue at pos 55 of block 1, protein 1",55,optAln[1][1][55]);
		assertEquals("Wrong residue at end of block 1, protein 0",181,optAln[1][0][117]);
		assertEquals("Wrong residue at end of block 1, protein 1",117,optAln[1][1][117]);

	}

	/**
	 * Print an AFPChain manually for debugging
	 * @param cpAlignment
	 */
	@SuppressWarnings("unused")
	private static void printOptAln(AFPChain cpAlignment) {
		int[] optLen = cpAlignment.getOptLen();
		int[][][] optAln = cpAlignment.getOptAln();
		
		for(int block=0;block<cpAlignment.getBlockNum();block++) {
			for(int pos=0;pos<optLen[block]; pos++) {
				System.out.format("%s\t%s\n", optAln[block][0][pos],
						optAln[block][1][pos]);
			}
			System.out.println();
		}
	}
}
