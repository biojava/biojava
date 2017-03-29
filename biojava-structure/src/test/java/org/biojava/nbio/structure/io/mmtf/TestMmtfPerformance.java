package org.biojava.nbio.structure.io.mmtf;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.TestStructureCrossReferences;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.junit.Test;
import org.rcsb.mmtf.dataholders.MmtfStructure;
import org.rcsb.mmtf.decoder.ReaderUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Created by andreas on 1/9/17.
 */
public class TestMmtfPerformance {

    private static final Logger logger = LoggerFactory.getLogger(TestMmtfPerformance.class);


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


    // Returns the contents of the file in a byte array.
    public static byte[] getBytesFromFile(File file) throws IOException {
        // Get the size of the file
        long length = file.length();

        // You cannot create an array using a long type.
        // It needs to be an int type.
        // Before converting to an int type, check
        // to ensure that file is not larger than Integer.MAX_VALUE.
        if (length > Integer.MAX_VALUE) {
            // File is too large
            throw new IOException("File is too large!");
        }

        // Create the byte array to hold the data
        byte[] bytes = new byte[(int)length];

        // Read in the bytes
        int offset = 0;
        int numRead = 0;

        InputStream is = new FileInputStream(file);
        try {
            while (offset < bytes.length
                    && (numRead=is.read(bytes, offset, bytes.length-offset)) >= 0) {
                offset += numRead;
            }
        } finally {
            is.close();
        }

        // Ensure all the bytes have been read in
        if (offset < bytes.length) {
            throw new IOException("Could not completely read file "+file.getName());
        }
        return bytes;
    }

    static String convertStreamToString(java.io.InputStream is) {
        java.util.Scanner s = new java.util.Scanner(is).useDelimiter("\\A");
        return s.hasNext() ? s.next() : "";
    }


    /** loads both pdb and mmtf file for 4CUP in memory and compares parsing performance.
     * That means any IO is excluded from the measurement and we really compare raw parsing speed.
     *
     * @throws IOException
     */
    @Test
    public void test4CUP() throws IOException{

        ClassLoader classLoader = getClass().getClassLoader();

        Path p = Paths.get(classLoader.getResource("org/biojava/nbio/structure/io/mmtf/4CUP.mmtf").getPath());
        byte[] mmtfBytes = Files.readAllBytes(p);

        InputStream mmtfIS = new ByteArrayInputStream(mmtfBytes);

        long mmtfStart = System.currentTimeMillis();
        Structure structure = MmtfActions.readFromInputStream(mmtfIS);
        long mmtfEnd = System.currentTimeMillis();

        assertEquals(structure.getPDBCode(), "4CUP");
        assertEquals(structure.getChains().size(), 6);


        ///end of mmtf parsing. Now we parse PDB:


        Path path = Paths.get(classLoader.getResource("org/biojava/nbio/structure/io/4cup.pdb.gz").getPath());

        InputStream is = new GZIPInputStream(Files.newInputStream(path));

        String pdbFile = convertStreamToString(is);

        long pdbStart = System.currentTimeMillis();
        PDBFileParser parser = new PDBFileParser();

        Structure s = parser.parsePDBFile(new ByteArrayInputStream(pdbFile.getBytes()));

        long pdbEnd = System.currentTimeMillis();

        //todo: add mmcif for comparison

        logger.warn("time to parse mmtf:" + (mmtfEnd-mmtfStart));
        logger.warn("time to parse PDB: " + (pdbEnd-pdbStart));

        assertTrue( "It should not be the case, but it is faster to parse a PDB file ("+(pdbEnd -pdbStart)+" ms.) than MMTF ("+( mmtfEnd-mmtfStart)+" ms.)!",( pdbEnd -pdbStart) > ( mmtfEnd-mmtfStart));

    }
}
