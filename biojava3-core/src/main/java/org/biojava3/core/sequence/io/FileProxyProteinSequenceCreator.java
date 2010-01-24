/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import java.io.File;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.loader.SequenceFileProxyLoader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FileProxyProteinSequenceCreator implements SequenceCreatorInterface {

    CompoundSet compoundSet = null;
    File fastaFile = null;

    public FileProxyProteinSequenceCreator(File fastaFile,CompoundSet compoundSet) {
        this.compoundSet = compoundSet;
        this.fastaFile = fastaFile;
    }

    public AbstractSequence getSequence(String sequence, long index) {
        SequenceFileProxyLoader sequenceFileProxyLoader = new SequenceFileProxyLoader(fastaFile,new FastaSequenceParser(),index,sequence.length(),compoundSet);
        return new ProteinSequence(sequence, compoundSet);
    }

    public AbstractSequence getSequence(SequenceProxyLoader proxyLoader, long index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
