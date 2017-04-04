package org.biojava.nbio.genome.parsers.twobit;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Created by yana on 4/4/17.
 */
public class SimpleTwoBitFileProvider {

    private File twoBitFileLocalLocation;
    private String genomeAssembly;

    public SimpleTwoBitFileProvider(File twoBitFileLocalLocation, String genomeAssembly) throws MalformedURLException {

        this.twoBitFileLocalLocation = twoBitFileLocalLocation;
        this.genomeAssembly = genomeAssembly;

        if ( ! twoBitFileLocalLocation.exists() ) {

            URL twoBitFileURL = getTwoBitURL(genomeAssembly);



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
