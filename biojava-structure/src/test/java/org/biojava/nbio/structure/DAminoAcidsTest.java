package org.biojava.nbio.structure;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.chem.ChemCompProvider;
import org.biojava.nbio.structure.chem.ChemCompTools;
import org.biojava.nbio.structure.chem.DownloadChemCompProvider;
import org.biojava.nbio.structure.io.CifFileReader;
import org.junit.jupiter.api.Test;

class DAminoAcidsTest {

	@Test
	public void testRecognizeDAminoAcids() throws IOException{
		
//		ChemCompTools.getPolymerType()
		DownloadChemCompProvider.getLocalFileName("ALA");
		DownloadChemCompProvider.getLocalFileName("DAL");
		DownloadChemCompProvider downloadChemCompProvider = new DownloadChemCompProvider();
		downloadChemCompProvider.getChemComp("ALA");
		downloadChemCompProvider.getChemComp("DAL");
		
		System.out.println(downloadChemCompProvider.getChemComp("ALA").getType());
		System.out.println(downloadChemCompProvider.getChemComp("SER").getType());
		System.out.println(downloadChemCompProvider.getChemComp("DAL").getType());
		System.out.println(downloadChemCompProvider.getChemComp("DSN").getType());
		
        InputStream cifStream = new GZIPInputStream(getClass().getResourceAsStream("/org/biojava/nbio/structure/io/1bck.cif.gz"));
        Structure structure= new CifFileReader().getStructure(cifStream);
        final Chain chainC = structure.getPolyChainByPDB("C");
		Group group = chainC.getAtomGroup(0);
		
		assertTrue(group.isAminoAcid(), "Not recognized as AminoAcid");
		assertTrue(group.isHetAtomInFile(), "Not (internally) recognized as HetAtomInFile");
		assertTrue(group.isPolymeric(), "Not recognized as Polymeric");
		assertFalse(group.isNucleotide(), "Group recognized as Neucleotide");
		assertFalse(group.isWater(), "Group recognized as water");
		
		assertTrue(group instanceof AminoAcid);
		AminoAcid aa = (AminoAcid) group;
		//test all AminoAcid methods
		aa.getAminoType();
		assertNotNull(aa.getCA());
		assertNotNull(aa.getC());
		assertNotNull(aa.getN());
		assertNotNull(aa.getO());
		assertEquals(AminoAcid.ATOMRECORD, aa.getRecordType());
	}
	
	@Test
	void testDAminoAcidNames() throws Exception {
		assertEquals("GLY", StructureTools.getChiralImage("Gly"), "Couldn't hanle GLY name");
		assertEquals("GLY", StructureTools.getDChiralImage("Gly"), "Couldn't hanle GLY name");
		assertEquals("GLY", StructureTools.getLChiralImage("Gly"), "Couldn't hanle GLY name");

		assertEquals("DAL", StructureTools.getDChiralImage("ALA"), "Couldn't find Ala D image");
		assertEquals("DSN", StructureTools.getDChiralImage("SER"), "Couldn't find Ser D image");
	
		assertEquals("ALA", StructureTools.getLChiralImage("DAL"), "Couldn't find Ala");
		assertEquals("SER", StructureTools.getLChiralImage("DSN"), "Couldn't find Ser");
		
		assertThrows(IllegalArgumentException.class, () ->{
			StructureTools.getChiralImage(null);
		});
		
		assertNull(StructureTools.getChiralImage("wrongValue"));
	}
}
