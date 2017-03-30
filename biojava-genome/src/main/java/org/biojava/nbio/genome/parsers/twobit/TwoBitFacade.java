package org.biojava.nbio.genome.parsers.twobit;

import java.io.File;

/** A facade that makes it easier to work with a 2bit file.
 *
 * Created by yana on 3/27/17.
 */
public class TwoBitFacade {

    private TwoBitParser twoBitParser = null;


    /**
     *  Reads a genome from a locally stored .2bit file.
     *
     *  @param file the File to a .2bit file.
     */
    public TwoBitFacade(File file) throws Exception {
        twoBitParser = new TwoBitParser(file);
    }

    /**
     *  Closes .2bit file twoBitParser.
     */
    public void close() throws Exception {
        if (twoBitParser != null)
            twoBitParser.close();

    }

    /** Sets a chromosome for TwoBitParser.
     *
     * @param chr The chromosome name (e.g. chr21)
     */
    public void setChromosome(String chr) throws Exception {
        if ( twoBitParser == null){

        }
        twoBitParser.close();
        String[] names = twoBitParser.getSequenceNames();
        for(int i=0;i<names.length;i++) {
            if ( names[i].equalsIgnoreCase(chr) ) {
                twoBitParser.setCurrentSequence(names[i]);
                break;
            }
        }
    }

    /** Extract a sequence from a chromosome, using chromosomal coordinates
     *
     * @param chromosomeName
     * @param start
     * @param end
     * @return the DNASequence from the requested coordinates.
     * @throws Exception
     */
    public String getSequence(String chromosomeName, int start, int end) throws Exception {
        twoBitParser.close();
        twoBitParser.setCurrentSequence(chromosomeName);
        return twoBitParser.loadFragment(start,end-start);
    }
}
