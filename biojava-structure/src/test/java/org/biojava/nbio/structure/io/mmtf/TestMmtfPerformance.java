package org.biojava.nbio.structure.io.mmtf;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.junit.Test;
import org.rcsb.mmtf.dataholders.MmtfStructure;
import org.rcsb.mmtf.decoder.ReaderUtils;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by andreas on 1/9/17.
 */
public class TestMmtfPerformance {

//    @Test
//    public void test3J3Q() throws IOException{
//
////        AllChemCompProvider cc = new AllChemCompProvider();
////        ChemCompGroupFactory.setChemCompProvider(cc);
//
//        long timeS = System.currentTimeMillis();
//        ClassLoader classLoader = getClass().getClassLoader();
//        Structure structure = MmtfActions.readFromFile((Paths.get(classLoader.getResource("org/biojava/nbio/structure/io/mmtf/3J3Q.mmtf").getPath())));
//        assertEquals(structure.getPDBCode(),"3J3Q");
//        //assertEquals(structure.getChains().size(),6);
//        long timeE = System.currentTimeMillis();
//
//        System.out.println("time to load from local file: " + (timeE - timeS) + " ms.");
//
//    }

    @Test
    public void test4CUP() throws IOException{

        long timeS = System.currentTimeMillis();

        ClassLoader classLoader = getClass().getClassLoader();
        Structure structure = MmtfActions.readFromFile((Paths.get(classLoader.getResource("org/biojava/nbio/structure/io/mmtf/4CUP.mmtf").getPath())));

        assertEquals(structure.getPDBCode(),"4CUP");
        assertEquals(structure.getChains().size(),6);

        long timeE = System.currentTimeMillis();

        Path path = Paths.get(classLoader.getResource("org/biojava/nbio/structure/io/4cup.pdb").getPath());

        InputStream is = Files.newInputStream(path);

        PDBFileParser parser = new PDBFileParser();

        Structure s = parser.parsePDBFile(is);

        long timeF = System.currentTimeMillis();

        //todo: add mmcif for comparison

//        System.out.println("time to parse mmtf:" + (timeE-timeS));
//        System.out.println("time to parse PDB: " + (timeF-timeE));

        assertTrue( "It should not be the case, but it is faster to parse a PDB file ("+(timeF -timeE)+" ms.) than MMTF ("+( timeE-timeS)+" ms.)!",( timeF -timeE) > ( timeE-timeS));

    }
}
