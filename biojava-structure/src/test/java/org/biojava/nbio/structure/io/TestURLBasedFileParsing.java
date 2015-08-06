package org.biojava.nbio.structure.io;

import junit.framework.TestCase;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;

/**
 * Created by ap3 on 31/07/2015.
 */
public class TestURLBasedFileParsing extends TestCase{

    @Test
    public void testMMcifURL(){

        String u = "http://ftp.wwpdb.org/pub/pdb/data/biounit/mmCIF/divided/nw/4nwr-assembly1.cif.gz";

        try {
            Structure s = StructureIO.getStructure(u);

            System.out.println(s);
        } catch (Exception e) {
            e.printStackTrace();
        }


    }
}
