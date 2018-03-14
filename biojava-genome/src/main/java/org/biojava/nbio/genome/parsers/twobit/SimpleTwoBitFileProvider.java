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
package org.biojava.nbio.genome.parsers.twobit;

import org.biojava.nbio.core.util.FileDownloadUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Created by yana on 4/4/17.
 */
public class SimpleTwoBitFileProvider {
    private static final Logger logger = LoggerFactory.getLogger(SimpleTwoBitFileProvider.class);

    public static synchronized void downloadIfNoTwoBitFileExists(File twoBitFileLocalLocation, String genomeAssembly) throws IOException {

        if ( ! twoBitFileLocalLocation.exists() ) {

            // download to a temporary file
            File tmp = File.createTempFile(genomeAssembly,".2bit");
            URL twoBitFileURL = getTwoBitURL(genomeAssembly);

            logger.info("downloading " + twoBitFileURL + " to " + tmp.getAbsolutePath());

            // 2bit files are large and take a while to download
            FileDownloadUtils.downloadFile(twoBitFileURL, tmp);

            // check the parent directory exists

            Path p = Paths.get(twoBitFileLocalLocation.getAbsolutePath());

            Path dir = p.getParent();
            if (! Files.exists(dir)) {
                Files.createDirectories(dir);
            }

            logger.info("renaming " + tmp.getAbsolutePath() +" to " + twoBitFileLocalLocation.getAbsolutePath());
            // after the download rename
            tmp.renameTo(twoBitFileLocalLocation);

        }
    }

    public static URL getTwoBitURL(String genomeAssembly) throws MalformedURLException {

        String url="";
        if (genomeAssembly.equals("hg19") || genomeAssembly.equals("hg37") ) {
            url = "http://cdn.rcsb.org/gene/hg37/hg19.2bit";
        }
        else if (genomeAssembly.equals("hg38")) {
            url = "http://cdn.rcsb.org/gene/hg38/hg38.2bit";
        }
        return new URL(url);
    }

    public static void main(String[] args){
        try {
            downloadIfNoTwoBitFileExists(new File("/Users/yana/spark/2bit/hg38.2bit"),"hg38");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
