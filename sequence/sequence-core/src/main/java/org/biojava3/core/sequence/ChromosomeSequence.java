/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

/**
 *
 * @author Scooter
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
