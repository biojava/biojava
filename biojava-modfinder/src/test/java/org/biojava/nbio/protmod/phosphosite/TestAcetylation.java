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
package org.biojava.nbio.protmod.phosphosite;


import org.biojava.nbio.phosphosite.Dataset;
import org.biojava.nbio.phosphosite.Site;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.List;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;


/** 
 * Makes sure there is a local installation of the Acetylation site file from Phosphosite and
 * tests if it can get parsed by the parser.
 *
 * Created by andreas on 11/29/16.
 */
public class TestAcetylation  {


    /** 
     * Make sure an Acetylation file is available locally.
     * Downloads from Phosphosite if needed.
     *
     */
    @Before
    public void setUp() throws IOException{

        Dataset ds = new Dataset();

        String f = Dataset.ACETYLATION;

        File localFile = getLocalFileName(f);

        if (!localFile.exists()) {
            ds.downloadFile(new URL(f), localFile);
        }

    }

    /** 
     * Returns the local file name where the Acetylation file will get cached locally.
     *
     * @param phosphoSiteFileLocation location of file at PhosphoSitePlus.
     * @return a File pointing to the location of the locally cached file.
     */
    private File getLocalFileName(String phosphoSiteFileLocation){

        Dataset ds = new Dataset();
        File localDir = ds.getLocalDir();
        if ( ! localDir.exists()) {
            boolean success = localDir.mkdir();
            if ( ! success)
                fail("Could not create directory " + localDir.getAbsolutePath());
        }

        int slashIndex = phosphoSiteFileLocation.lastIndexOf("/");

        String fileName = phosphoSiteFileLocation.substring(slashIndex);

        return new File(localDir + "/" + fileName);
    }

    /** 
     * Tests that the acetylation file can get parsed without problems.
     *
     */
    @Test
    public void testAcetylation() throws IOException {

        File localFile = getLocalFileName(Dataset.ACETYLATION);

        List<Site> sites = Site.parseSites(localFile);

        assertTrue(sites.size() > 0);

        for (Site s : sites) {

            assertTrue(s.getResidue() != null);

        }

    }
}
