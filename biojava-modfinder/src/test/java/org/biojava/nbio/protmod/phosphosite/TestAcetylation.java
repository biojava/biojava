package org.biojava.nbio.protmod.phosphosite;

import junit.framework.TestCase;
import org.biojava.nbio.phosphosite.Dataset;
import org.biojava.nbio.phosphosite.Site;

import java.io.File;
import java.net.URL;
import java.util.List;

/**
 * Created by andreas on 11/29/16.
 */
public class TestAcetylation extends TestCase {


    /** Tests that the acetylation file can get downloaded and parsed
     *
     */
    public void testAcetylation() {
        String f = Dataset.ACETYLATION;

        Dataset ds = new Dataset();

        try {
            File localDir = ds.getLocalDir();
            int slashIndex = f.lastIndexOf("/");

            String fileName = f.substring(slashIndex);

            File localFile = new File(localDir + "/" + fileName);

            if (!localFile.exists()) {
                ds.downloadFile(new URL(f), localFile);
            }

            List<Site> sites = Site.parseSites(localFile);

            assertTrue(sites.size() > 0);

            for (Site s : sites) {

                assertTrue(s.getResidue() != null);

            }


        } catch (Exception e) {
            e.printStackTrace();
            fail(e.getMessage());
        }
    }
}
