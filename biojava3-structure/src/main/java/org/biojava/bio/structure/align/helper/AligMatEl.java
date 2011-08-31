/*
 *                  BioJava development code
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
 * Created on Mar 1, 2006
 *
 */
package org.biojava.bio.structure.align.helper;

public class AligMatEl 
extends IndexPair{

   
   
        /**
	 * 
	 */
	private static final long serialVersionUID = -4040926588803887471L;
		int value;       
        int contig;
       
        
        public AligMatEl(){
            super();
            value  = -1;
            contig = -1;
        }

        public int getContig() {
            return contig;
        }
        public void setContig(int contig) {
            this.contig = contig;
        }
     
        public int getValue() {
            return value;
        }
        public void setValue(int value) {
            this.value = value;
        }
        
        public String toString(){
            String ret = "AligMatEl val:" + value + " contig:" + contig + 
            " trackrow:" + getRow() + " trackcol:" + getCol();
            return ret;
        }
        
    }

    


