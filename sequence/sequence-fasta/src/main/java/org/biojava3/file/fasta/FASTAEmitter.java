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
import java.util.Iterator;
import org.biojava3.core.io.ThingEmitter;
import org.biojava3.core.io.ThingReceiver;

/**
 * Takes a FASTA object and emits its data as if it were being parsed from a
 * data source.
 * @author Richard Holland
 * @since 3.0
 */
public interface FASTAEmitter extends ThingEmitter<FASTA> {
    
    /**
     * A standard implementation which emits FASTA contents.
     */
    public static class FASTAEmitterImpl implements FASTAEmitter {
        private Iterator<FASTA> things;
        private FASTAReceiver fastaReceiver;

        public void setThingSource(Iterator<FASTA> things) {
            if (things==null) {
                throw new NullPointerException("Cannot set null iterator of things");
            }
            this.things = things;
        }

        public void setNextThingReceiver(ThingReceiver thingReceiver) {
            if (!(thingReceiver instanceof FASTAReceiver)) {
                throw new IllegalArgumentException("Cannot set a receiver which is not a FASTAReceiver");
            }
            this.fastaReceiver = (FASTAReceiver)thingReceiver;
        }

        public boolean canReadNextThing() throws IOException {
            return this.things.hasNext();
        }

        public void readNextThing() throws IOException {
            FASTA fasta = this.things.next();
            this.fastaReceiver.startThing();
            this.fastaReceiver.setDescriptionLine(fasta.getDescriptionLine());
            this.fastaReceiver.setSequence(fasta.getSequence());
            this.fastaReceiver.finishThing();
        }        
        
        public void close() throws IOException {
        }
    }
}
