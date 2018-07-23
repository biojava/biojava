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
package org.biojava.nbio.structure.io.mmtf;


import org.junit.Test;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.BondImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.ExperimentalTechnique;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.io.mmtf.MmtfUtils;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.xtal.BravaisLattice;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
/**
 * Test the MMTF utils class
 * @author Anthony Bradley
 *
 */
public class TestMmtfUtils {

	/**
	 * Integration test to see that the microheterogenity is being dealt with correctly.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void microHeterogenity() throws IOException, StructureException {
		MmtfUtils.setUpBioJava();
		Structure inputStructure = StructureIO.getStructure("4ck4");
		// Count the number of groups
		Group before = inputStructure.getChains().get(0).getAtomGroup(17);
		assertTrue(inputStructure.getChains().get(0).getAtomGroup(17).hasAltLoc());
		List<Atom> totalAtoms = new ArrayList<>(getAllAtoms(inputStructure));
		int totGroups = 0;
		int totAtomsCounter = 0;
		Set<Atom> totAtoms = new HashSet<>();
		for (Chain c : inputStructure.getChains()) {
			totGroups += c.getAtomGroups().size();
			for (Group g: c.getAtomGroups() ){
				totAtomsCounter+=g.getAtoms().size();
				totAtoms.addAll(g.getAtoms());
				for (Group alt : g.getAltLocs()) {
					totAtomsCounter+=alt.getAtoms().size();
					totAtoms.addAll(alt.getAtoms());
				}
			}
		}
		// Now "fix" the microheterogenity
		MmtfUtils.fixMicroheterogenity(inputStructure);
		assertEquals(before, inputStructure.getChains().get(0).getAtomGroup(17));
		assertFalse(inputStructure.getChains().get(0).getAtomGroup(17).hasAltLoc());
		assertFalse(inputStructure.getChains().get(0).getAtomGroup(18).hasAltLoc());
		int totGroupsAfter = 0;
		int totAtomsCounterAfter = 0;
		Set<Atom> totAtomsAfter = new HashSet<>();
		for (Chain c : inputStructure.getChains()) {
			totGroupsAfter += c.getAtomGroups().size();
			for (Group g: c.getAtomGroups() ){
				totAtomsCounterAfter+=g.getAtoms().size();
				totAtomsAfter.addAll(g.getAtoms());
				for (Group alt : g.getAltLocs()) {
					totAtomsAfter.addAll(alt.getAtoms());
					totAtomsCounterAfter+=alt.getAtoms().size();
				}
			}
		}
		// Find the atoms after the fix.
		List<Atom> totalAtomsAfter = new ArrayList<>(getAllAtoms(inputStructure));
		// Get all of the duplicate atoms
		Set<Atom> duplicates = findDuplicates(totalAtomsAfter);
		for (Atom a : duplicates) {
			System.out.println(a);
		}
		// There should be no duplicates
		assertEquals(duplicates.size(), 0);
		assertEquals(totalAtoms.size(), totalAtomsAfter.size());
		// Check there are two more groups afterwards
		assertEquals(totGroupsAfter-2, totGroups);
		// Check there are no more atoms afterwards
		assertEquals(totAtomsAfter.size(), totAtoms.size());
		// Check the counter are the same too
		assertEquals(totAtomsCounterAfter, totAtomsCounter);

	}

	/**
	 * Function to get all the atoms in the strucutre as a list.
	 * @param bioJavaStruct the biojava structure
	 * @return a list of all the unique atoms in the structure
	 */
	private List<Atom> getAllAtoms(Structure bioJavaStruct) {
		// Get all the atoms
		List<Atom> theseAtoms = new ArrayList<Atom>();
		for (int i=0; i<bioJavaStruct.nrModels(); i++){
			List<Chain> chains = bioJavaStruct.getModel(i);
			for (Chain c : chains) {
				for (Group g : c.getAtomGroups()) {
					for(Atom a: MmtfUtils.getAtomsForGroup(g)){
						theseAtoms.add(a);					
					}
				}
			}
		}
		return theseAtoms;
	}

	/**
	 * Test that getting the space group info as a string works.
	 */
	@Test
	public void testGetSpaceGroupAsString() {
		assertEquals("NA", MmtfUtils.getSpaceGroupAsString(null));
		SpaceGroup spaceGroup = new SpaceGroup(21, 1, 1, "P212121", "P 21 21 21", BravaisLattice.TRICLINIC);
		assertEquals("P212121", MmtfUtils.getSpaceGroupAsString(spaceGroup));
	}

	/**
	 * Test that getting the unit cell as an array of doubles works.
	 */
	@Test
	public void testGetUnitCellAsArray() {
		PDBCrystallographicInfo xtalInfo = new PDBCrystallographicInfo();
		CrystalCell cell = new CrystalCell();
		cell.setA(1.0);
		cell.setB(2.0);
		cell.setC(3.0);
		cell.setAlpha(4.0);
		cell.setBeta(5.0);
		cell.setGamma(6.0);
		float[] testArray = new float[] {1.0f,2.0f,3.0f,4.0f,5.0f,6.0f};
		xtalInfo.setCrystalCell(cell);
		float[] outputArray = MmtfUtils.getUnitCellAsArray(xtalInfo);
		assertArrayEquals(testArray, outputArray, 0.0f);
	}

	/**
	 * Test getting the list of experimental methods as string array.
	 */
	@Test
	public void testGetExperimentalMethods() {
		Set<ExperimentalTechnique> experimentalTechniques = new HashSet<>();
		experimentalTechniques.add(ExperimentalTechnique.XRAY_DIFFRACTION);
		experimentalTechniques.add(ExperimentalTechnique.ELECTRON_MICROSCOPY);
		String[] techniques = MmtfUtils.techniquesToStringArray(experimentalTechniques);
		String[] testTechniques = {"X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"};
		Arrays.sort(techniques);
		Arrays.sort(testTechniques);
		assertArrayEquals(testTechniques, techniques);
	}

	/**
	 * Test the conversion of a matrix to an array of doubles.
	 */
	@Test 
	public void testConvertToDoubleArray() {
		Matrix4d matrix4d = new Matrix4d();
		matrix4d.m00 = 0.0;
		matrix4d.m01 = 0.1;
		matrix4d.m02 = 0.2;
		matrix4d.m03 = 0.3;
		matrix4d.m10 = 1.0;
		matrix4d.m11 = 1.1;
		matrix4d.m12 = 1.2;
		matrix4d.m13 = 1.3;
		matrix4d.m20 = 2.0;
		matrix4d.m21 = 2.1;
		matrix4d.m22 = 2.2;
		matrix4d.m23 = 2.3;
		matrix4d.m30 = 3.0;
		matrix4d.m31 = 3.1;
		matrix4d.m32 = 3.2;
		matrix4d.m33 = 3.3;
		double[] testData = new double[] {0.0, 0.1, 0.2, 0.3,
				1.0, 1.1, 1.2, 1.3,
				2.0, 2.1, 2.2, 2.3,
				3.0, 3.1, 3.2, 3.3};
		assertArrayEquals(testData,MmtfUtils.convertToDoubleArray(matrix4d), 0.0);
	}

	/**
	 * Test to check the conversion of BioassemblyInfo to a primitive map.
	 */
	@Test
	public void testMakePrimitiveBioasembly() {
		double[] testData = new double[] {0.0, 0.1, 0.2, 0.3,
				1.0, 1.1, 1.2, 1.3,
				2.0, 2.1, 2.2, 2.3,
				3.0, 3.1, 3.2, 3.3};
		BioAssemblyInfo bioAssemblyInfo = new BioAssemblyInfo();
		List<BiologicalAssemblyTransformation> transforms = new ArrayList<>();
		BiologicalAssemblyTransformation biologicalAssemblyTransformation = new BiologicalAssemblyTransformation();
		biologicalAssemblyTransformation.setChainId("C");
		biologicalAssemblyTransformation.setTransformationMatrix(new Matrix4d(testData));
		transforms.add(biologicalAssemblyTransformation);
		bioAssemblyInfo.setTransforms(transforms);
		// Map the chain to the second index
		Map<String, Integer> chainIdToIndexMap = new HashMap<>();
		chainIdToIndexMap.put("C", 2);
		// Now do the conversion and test they are the same
		Map<double[], int[]> transMap = MmtfUtils.getTransformMap(bioAssemblyInfo, chainIdToIndexMap);
		assertArrayEquals(testData, (double[]) transMap.keySet().toArray()[0], 0.0);
		assertArrayEquals(new int[] {2} , (int[]) transMap.values().toArray()[0]);
	}


	/**
	 * Test getting the data as an appropriately formatted string.
	 */
	public void testGetIsoDateString() {
		Date inputDate = new Date();
		inputDate.setTime(86500);
		// One day after 
		assertEquals("1970-01-02",MmtfUtils.dateToIsoString(inputDate));
	}

	/**
	 * Test getting the number of groups from a structure.
	 */
	@Test
	public void testGetNumGroups() {
		Structure structure = new StructureImpl();
		Chain chain = new ChainImpl();
		Group groupOne = new  AminoAcidImpl();
		Group groupTwo = new HetatomImpl();
		Group groupThree = new NucleotideImpl();
		structure.addChain(chain);
		chain.addGroup(groupOne);
		chain.addGroup(groupTwo);
		chain.addGroup(groupThree);
		assertEquals(3,MmtfUtils.getNumGroups(structure));
	}


	/**
	 * Test getting the correct atoms from a group
	 */
	@Test
	public void testGetAtomsForGroup() {
		Group group = new AminoAcidImpl();
		Group altLoc = new AminoAcidImpl();
		Atom atomOne = new AtomImpl();
		atomOne.setX(1.00);
		Atom atomTwo = new AtomImpl();
		atomTwo.setX(2.00);
		Atom atomThree = new AtomImpl();
		atomThree.setX(3.00);
		atomThree.setAltLoc('B');
		Atom atomFour = new AtomImpl();
		atomFour.setX(4.00);
		List<Atom> inputList = new ArrayList<>();
		inputList.add(atomOne);
		inputList.add(atomTwo);
		inputList.add(atomFour);
		inputList.add(atomThree);
		group.addAtom(atomOne);
		group.addAtom(atomTwo);
		group.addAtom(atomFour);
		altLoc.addAtom(atomOne);
		altLoc.addAtom(atomTwo);
		altLoc.addAtom(atomThree);
		group.addAltLoc(altLoc);
		List<Atom> atomList = MmtfUtils.getAtomsForGroup(group);
		assertEquals(inputList, atomList);
	}


	/**
	 * Test getting the number of bonds from a list of atoms.
	 */
	@Test
	public void testGetNumBondsFromGroup() {
		List<Atom> atoms = new ArrayList<>();
		Atom atomOne = new AtomImpl();
		Atom atomTwo = new AtomImpl();
		Atom atomThree = new AtomImpl();
		atoms.add(atomOne);
		atoms.add(atomTwo);
		atoms.add(atomThree);
		// Make the same bond twice iwth different atom orders
		new BondImpl(atomOne, atomTwo, 2);
		new BondImpl(atomTwo, atomOne, 2);
		// Make the same bond twice
		new BondImpl(atomOne, atomThree, 2);
		new BondImpl(atomOne, atomThree, 2);
		// Make this bond twice with different orders
		new BondImpl(atomTwo, atomThree, 2);		
		new BondImpl(atomTwo, atomThree, 1);
		assertEquals(3, MmtfUtils.getNumBondsInGroup(atoms));
	}
	/**
	 * Test that getting the secondary structure type works.
	 */
	@Test
	public void testGetSetSecStructType() {
		Group group = new AminoAcidImpl();
		MmtfUtils.setSecStructType(group, 0);
		assertEquals(MmtfUtils.getSecStructType(group), 0);
		MmtfUtils.setSecStructType(group, 1);
		assertEquals(MmtfUtils.getSecStructType(group), 1);
		MmtfUtils.setSecStructType(group, 2);
		assertEquals(MmtfUtils.getSecStructType(group), 2);
		MmtfUtils.setSecStructType(group, 3);
		assertEquals(MmtfUtils.getSecStructType(group), 3);
		MmtfUtils.setSecStructType(group, 4);
		assertEquals(MmtfUtils.getSecStructType(group), 4);
		MmtfUtils.setSecStructType(group, 5);
		assertEquals(MmtfUtils.getSecStructType(group), 5);
		MmtfUtils.setSecStructType(group, 6);
		assertEquals(MmtfUtils.getSecStructType(group), 6);
		MmtfUtils.setSecStructType(group, 7);
		assertEquals(MmtfUtils.getSecStructType(group), 7);
		// Now test two null possibilities
		Group newGroup = new AminoAcidImpl();
		MmtfUtils.setSecStructType(newGroup, -1);
		assertEquals(MmtfUtils.getSecStructType(newGroup), -1);	
		// Now test two null possibilities
		Group newerGroup = new AminoAcidImpl();
		MmtfUtils.setSecStructType(newerGroup, 10);
		assertEquals(MmtfUtils.getSecStructType(newerGroup), -1);	
	}

	/**
	 * Test that setting the secondary structure types behaves as expected.
	 */
	@Test
	public void testGetSecStructTypeFromDsspIndex(){
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(0).name,"pi Helix");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(1).name,"Bend");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(2).name,"alpha Helix");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(3).name,"Extended");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(4).name,"3-10 Helix");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(5).name,"Bridge");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(6).name,"Turn");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(7).name,"Coil");
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(-1), null);
		assertEquals(MmtfUtils.getSecStructTypeFromDsspIndex(10), null);

	}

	/**
	 * Test that getting the structure data info works.
	 */
	@Test
	public void testGetStructureInfo() {
		Structure structure = new StructureImpl();
		Chain chain = new ChainImpl();
		chain.setId("A");
		Map<String,Integer> testMap = new HashMap<>();
		testMap.put("A", 0);
		List<Chain> chainList = new ArrayList<>();
		chainList.add(chain);
		Group group = new AminoAcidImpl();
		chain.addGroup(group);
		Atom atomOne = new AtomImpl();
		Atom atomTwo = new AtomImpl();
		List<Atom> atomList = new ArrayList<>();
		atomList.add(atomOne);
		atomList.add(atomTwo);
		new BondImpl(atomOne, atomTwo, 1);
		structure.addChain(chain);
		group.addAtom(atomOne);
		group.addAtom(atomTwo);
		// Get the structure
		MmtfSummaryDataBean mmtfSummaryDataBean = MmtfUtils.getStructureInfo(structure);
		assertEquals(mmtfSummaryDataBean.getAllAtoms(), atomList);
		assertEquals(testMap, mmtfSummaryDataBean.getChainIdToIndexMap());
		assertEquals(chainList, mmtfSummaryDataBean.getAllChains());
		assertEquals(1, mmtfSummaryDataBean.getNumBonds());
	}

	private Set<Atom> findDuplicates(List<Atom> listContainingDuplicates)
	{ 
		final Set<Atom> setToReturn = new HashSet<>(); 
		final Set<Atom> set1 = new HashSet<>();

		for (Atom yourInt : listContainingDuplicates)
		{
			if (!set1.add(yourInt))
			{
				setToReturn.add(yourInt);
			}
		}
		return setToReturn;
	}

	/**
	 * Test that the NCS data can be roundtripped.
	 */
	@Test
	public void testGetNcsMatrix() {
		double[][] testData = new double[][] {{1.0, 2.0,3.0,4.0,
			11.0,12.0,13.0,14.0,
			21.0,22.0,23.0,24.0,
			31.0,32.0,33.0,34.0}};
			testInput(testData);
	}

	/**
	 * Test that the NCS data can be roundtripped.
	 */
	@Test
	public void testEmptyNcsMatrix() {
		double[][] testData = new double[0][0];
		testInput(testData);
		double[][] output = MmtfUtils.getNcsAsArray(new Matrix4d[0]);
		assertNotNull(output);
	}

	/**
	 * Test what happens if the NCS is null
	 */
	@Test
	public void testNullNcsMatrix(){
		double[][] output = MmtfUtils.getNcsAsArray(null);
		assertNotNull(output);
		Matrix4d[] outputMat = MmtfUtils.getNcsAsMatrix4d(null);
		assertNull(outputMat);
		double[][] outputMatArr = MmtfUtils.getNcsAsArray(outputMat);
		assertNotNull(outputMatArr);
	}


	/**
	 * Test that the NCS data can be roundtripped - when two matrices are present.
	 */
	@Test
	public void testGetNcsMatrixHard() {
		double[][] testData = new double[][] {{1.0, 2.0,3.0,4.0,
			11.0,12.0,13.0,14.0,
			21.0,22.0,23.0,24.0,
			31.0,32.0,33.0,34.0,},{
				1.0, 2.0,3.0,4.0,
				11.0,12.0,13.0,14.0,
				21.0,22.0,23.0,24.0,
				31.0,32.0,33.0,34.0}};
				testInput(testData);
	}

	private void testInput(double[][] testData) {
		Matrix4d[] matArr = MmtfUtils.getNcsAsMatrix4d(testData);
		double[][] roundTrippedData = MmtfUtils.getNcsAsArray(matArr);
		for(int i=0; i<testData.length; i++){
			assertArrayEquals(testData[i], roundTrippedData[i], 0.0);		
		}
	}
}

