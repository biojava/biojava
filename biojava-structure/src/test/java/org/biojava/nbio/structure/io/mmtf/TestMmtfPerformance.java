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
import java.net.URL;
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

    private static final int NUMBER_OF_REPEATS = 10;

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


    public byte[] getByteArrayFromInputStream(InputStream is) throws IOException {
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();

        int nRead;
        byte[] data = new byte[16384];

        while ((nRead = is.read(data, 0, data.length)) != -1) {
            buffer.write(data, 0, nRead);
        }

        buffer.flush();

        return buffer.toByteArray();

    }

    @Test
    public void test3HBX() throws Exception{
        String pdbId = "3hbx";

        pdbId = pdbId.toUpperCase();

        URL url = new URL("https://files.rcsb.org/download/"+pdbId+".pdb.gz");

        String pdbFile = convertStreamToString(new GZIPInputStream(url.openStream()));

        long totalTimePDB = 0;
        long totalTimeMMTF = 0;

        byte[] pdbBytes = pdbFile.getBytes();

        PDBFileParser parser = new PDBFileParser();

        URL mmtfURL = new URL("https://mmtf.rcsb.org/v1.0/full/" + pdbId + ".mmtf.gz");

        byte[] mmtfdata = getByteArrayFromInputStream(new GZIPInputStream((mmtfURL.openStream())));

        for ( int i =0 ; i< NUMBER_OF_REPEATS ; i++) {

            long mmtfStart = System.nanoTime();
            Structure mmtfStructure = MmtfActions.readFromInputStream(new ByteArrayInputStream(mmtfdata));
            long mmtfEnd = System.nanoTime();



            long pdbStart = System.nanoTime();
            Structure pdbStructure = parser.parsePDBFile(new ByteArrayInputStream(pdbBytes));
            long pdbEnd = System.nanoTime();

            totalTimePDB += (pdbEnd - pdbStart);


            totalTimeMMTF += (mmtfEnd-mmtfStart);
        }


        long timePDB = (totalTimePDB/NUMBER_OF_REPEATS);
        long timeMMTF = (totalTimeMMTF/NUMBER_OF_REPEATS);


        logger.warn("average time to parse mmtf: " + timeMMTF/(1000*1000) + " ms.");
        logger.warn("average time to parse PDB : " + timePDB/(1000*1000) + " ms. ");
//
        assertTrue( "It should not be the case, but it is faster to parse a PDB file ("+timePDB+" ns.) than MMTF ("+( timeMMTF)+" ns.)!",( timePDB) > ( timeMMTF));
//
    }
}
