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
 * Created on DATE
 *
 */
package org.biojava3.core.sequence;

/**
 *
 * @author Scooter Willis
 */
public class ChromosomeSequence extends DNASequence {

   
    private int chromosomeNumber;




    /**
     * @return the chromosomeNumber
     */
    public int getChromosomeNumber() {
        return chromosomeNumber;
    }

    /**
     * @param chromosomeNumber the chromosomeNumber to set
     */
    public void setChromosomeNumber(int chromosomeNumber) {
        this.chromosomeNumber = chromosomeNumber;
    }
}
