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
package org.biojava.nbio.structure.io;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmcif.MMCIFFileTools;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.io.mmcif.model.CIFLabel;
import org.biojava.nbio.structure.io.mmcif.model.IgnoreField;
import org.junit.Test;

public class TestMMCIFWriting {

	@Test
	public void test1SMT() throws IOException, StructureException {
		// an x-ray structure
		testRoundTrip("1SMT");
	}

	/**
	 * MMCIF write test for an NMR structure with 2 chains
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test2N3J() throws IOException, StructureException {
		// an NMR structure (multimodel) with 2 chains
		testRoundTrip("2N3J");
	}

	@Test
	public void test1A2C() throws IOException, StructureException {
		// a structure with insertion codes
		testRoundTrip("1A2C");
	}

	private static class DemoBean {
		@IgnoreField
		String not_a_field;

		@SuppressWarnings("unused")//used by reflection
		String default_field;

		@CIFLabel(label="custom_label")
		String custom_field;

		public void setNot_a_field(String not_a_field) {
			this.not_a_field = not_a_field;
		}
		public void setDefault_field(String default_field) {
			this.default_field = default_field;
		}
		public void setCustom_field(String custom_field) {
			this.custom_field = custom_field;
		}
	}

	@Test
	public void testBeanAnnotations() {
		DemoBean bean = new DemoBean();
		bean.setCustom_field("custom_field");
		bean.setDefault_field(null);
		bean.setNot_a_field("not_a_field");


		// Test (1) should have custom_label (@CIFLabel)
		//      (2) shouldn't have not_a_field (@IgnoreField)
		String newline = System.getProperty("line.separator");
		String mmcif = MMCIFFileTools.toMMCIF("_demo", bean);
		String expected =
				  "_demo.default_field   ?" + newline
				+ "_demo.custom_label    custom_field" + newline
				+ "#" + newline;
		assertEquals(expected, mmcif);
	}

	private static void testRoundTrip(String pdbId) throws IOException, StructureException {
		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);

		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

		Structure originalStruct = StructureIO.getStructure(pdbId);

		File outputFile = File.createTempFile("biojava_testing_", ".cif");
		outputFile.deleteOnExit();


		FileWriter fw = new FileWriter(outputFile);
		fw.write(originalStruct.toMMCIF());
		fw.close();


		MMcifParser parser = new SimpleMMcifParser();

		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();

		FileParsingParameters fileParsingParams = new FileParsingParameters();
		fileParsingParams.setAlignSeqRes(true);

		consumer.setFileParsingParameters(fileParsingParams);

		parser.addMMcifConsumer(consumer);

		parser.parse(new BufferedReader(new FileReader(outputFile)));

		Structure readStruct = consumer.getStructure();

		assertNotNull(readStruct);

		assertEquals(originalStruct.getChains().size(), readStruct.getChains().size());

		assertEquals(originalStruct.nrModels(), readStruct.nrModels());

		for (int i=0; i<originalStruct.nrModels();i++) {
			assertEquals(originalStruct.getModel(i).size(), readStruct.getModel(i).size());
		}



		for (int modelIdx=0;modelIdx<originalStruct.nrModels();modelIdx++) {

			for (int i=0;i<originalStruct.getModel(modelIdx).size();i++) {
				assertEquals(originalStruct.getChains().get(i).getAtomGroups().size(),
								readStruct.getChains().get(i).getAtomGroups().size());

				Chain origChain = originalStruct.getModel(modelIdx).get(i);
				Chain readChain = readStruct.getModel(modelIdx).get(i);

				assertEquals(origChain.getAtomGroups().size(), readChain.getAtomGroups().size());
				//assertEquals(origChain.getSeqResGroups().size(), readChain.getSeqResGroups().size());
				assertEquals(origChain.getId(), readChain.getId());
				assertEquals(origChain.getName(), readChain.getName());

				Atom[] origAtoms = StructureTools.getAllAtomArray(origChain);
				Atom[] readAtoms = StructureTools.getAllAtomArray(readChain);

				assertEquals(origAtoms.length, readAtoms.length);

				for (int atomIdx=0;atomIdx<origAtoms.length;atomIdx++) {

					assertEquals("atom serials don't match for atom "+origAtoms[atomIdx].toString(),
							origAtoms[atomIdx].getPDBserial(), readAtoms[atomIdx].getPDBserial());

					assertEquals("atom names don't match for atom "+origAtoms[atomIdx].toString(),
							origAtoms[atomIdx].getName(), readAtoms[atomIdx].getName());

					assertEquals("atom elements don't match for atom "+origAtoms[atomIdx].toString(),
							origAtoms[atomIdx].getElement(), readAtoms[atomIdx].getElement());

					assertEquals("x values don't match for atom "+origAtoms[atomIdx].toString(),
							origAtoms[atomIdx].getX(), readAtoms[atomIdx].getX(),0.0001);

					assertEquals("y values don't match for atom "+origAtoms[atomIdx].toString(),
							origAtoms[atomIdx].getY(), readAtoms[atomIdx].getY(),0.0001);

					assertEquals("z values don't match for atom "+origAtoms[atomIdx].toString(),
							origAtoms[atomIdx].getZ(), readAtoms[atomIdx].getZ(),0.0001);
				}
			}

		}

		// Test cell and symmetry
		assertEquals(originalStruct.getCrystallographicInfo().getSpaceGroup(),
				readStruct.getCrystallographicInfo().getSpaceGroup());


	}

	/**
	 * Tests that structures containing symmetry mates with modified chain identifiers
	 * can be written out correctly.
	 */
	@Test
	public void testBiounitWriting()  {
		Structure s = createDummyStructure();
		String mmcif = s.toMMCIF();
		String[] lines = mmcif.split("\n");
		long atomLines = Arrays.stream(lines).filter(l -> l.startsWith("ATOM")).count();
		assertNotNull(mmcif);
		assertEquals(4, atomLines);
	}

	private static Structure createDummyStructure() {
		Group g = new AminoAcidImpl();
		Atom a = getAtom("CA", Element.C, 1, 1, 1, 1);
		g.addAtom(a);
		g.setResidueNumber(new ResidueNumber("A", 1, null));
		Group altLocG = new AminoAcidImpl();
		Atom a2 = getAtom("CA", Element.C, 2, 2, 2, 2);
		altLocG.addAtom(a2);
		altLocG.setResidueNumber(new ResidueNumber("A", 1, null));

		g.addAltLoc(altLocG);

		Chain c1 = new ChainImpl();
		c1.addGroup(g);
		c1.setId("A");
		EntityInfo entityInfo = new EntityInfo();
		entityInfo.setMolId(1);
		entityInfo.addChain(c1);
		c1.setEntityInfo(entityInfo);

		Group gc2 = new AminoAcidImpl();
		Atom ac2 = getAtom("CA", Element.C, 3, 3, 3, 3);
		gc2.addAtom(ac2);
		gc2.setResidueNumber(new ResidueNumber("A_1", 1, null));

		Group altLocGc2 = new AminoAcidImpl();
		Atom ac22 = getAtom("CA", Element.C, 4, 4, 4, 4);
		altLocGc2.addAtom(ac22);
		altLocGc2.setResidueNumber(new ResidueNumber("A_1", 1, null));

		gc2.addAltLoc(altLocGc2);

		Chain c2 = new ChainImpl();
		c2.addGroup(gc2);
		c2.setId("A_1");
		c2.setEntityInfo(entityInfo);
		entityInfo.addChain(c2);

		Structure s = new StructureImpl();
		s.addChain(c1);
		s.addChain(c2);
		return s;
	}

	private static Atom getAtom(String name, Element e, int id, double x, double y, double z) {
		Atom a = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);
		a.setPDBserial(id);
		a.setName(name);
		a.setElement(e);
		return a;
	}
}
