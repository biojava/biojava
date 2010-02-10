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
package org.biojava3.file.fasta;

import java.io.IOException;
import org.biojava3.core.io.ThingBuilder;
import org.biojava3.file.fasta.FASTA;

/**
 * Receives data via the {@code set()} methods and uses it to construct
 * new FASTA objects.
 * @author Richard Holland
 * @since 3.0
 */
public interface FASTABuilder extends FASTAReceiver, ThingBuilder<FASTA> {
    
    /**
     * A utility implementation of FASTABuilder.
     */
    public static class FASTABuilderImpl implements FASTABuilder {

		private static final long serialVersionUID = 1L;
		
		private FASTA fasta;
        private StringBuilder sb = new StringBuilder();

        public String getDescriptionLine() {
            return this.fasta.getDescriptionLine();
        }

        public CharSequence getSequence() {
            return this.fasta.getSequence();
        }

        public void setDescriptionLine(String descLine) {
            this.fasta.setDescriptionLine(descLine);
        }

        public void setSequence(CharSequence seq) {
            if (seq instanceof StringBuilder) {
                this.sb = (StringBuilder)seq;
            } else {
                this.sb = new StringBuilder(seq);
            }
            this.fasta.setSequence(sb);
        }
    
        public void appendSequence(CharSequence seq) {
            this.sb.append(seq);
        }

        public void endGroupOfThings() {
            // Don't care.
        }

        public void finishThing() {
            // Don't care.
        }

        public void startGroupOfThings() {
            // Don't care.
        }

        public void startThing() {
            this.fasta = new FASTAImpl();
        }

        public FASTA getFinishedThing() {
            return this.fasta;
        }   
        
        public void close() throws IOException {
            // Don't care.
        }
    }
}
