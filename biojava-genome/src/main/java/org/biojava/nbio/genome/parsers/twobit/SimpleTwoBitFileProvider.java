package org.biojava.nbio.genome.parsers.twobit;

import org.biojava.nbio.core.util.FileDownloadUtils;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Created by yana on 4/4/17.
 */
public class SimpleTwoBitFileProvider {

    public static void downloadIfNoTwoBitFileExists(File twoBitFileLocalLocation, String genomeAssembly) throws IOException {

        if ( ! twoBitFileLocalLocation.exists() ) {

            // download to a temporary file
            File tmp = File.createTempFile("",".2bit");
            URL twoBitFileURL = getTwoBitURL(genomeAssembly);

            // 2bit files are large and take a while to download
            FileDownloadUtils.downloadFile(twoBitFileURL, tmp);

            // after the download rename
            tmp.renameTo(twoBitFileLocalLocation);

        }
    }

    public static URL getTwoBitURL(String genomeAssembly) throws MalformedURLException {

        String url="";
        if (genomeAssembly.equals("hg37")) {
            url = "http://cdn.rcsb.org//gene/hg37/hg19.2bit";
        }
        else if (genomeAssembly.equals("hg38")) {
            url = "http://cdn.rcsb.org//gene/hg38/hg38.2bit";
        }
        return new URL(url);
    }

}
