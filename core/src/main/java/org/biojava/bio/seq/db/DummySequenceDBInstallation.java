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

package org.biojava.bio.seq.db;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * <code>DummySequenceDBInstallation</code> is an implementation which
 * returns the same <code>DummySequenceDB</code> instance regardless
 * of the identifier used to retrieve a database.
 *
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @since 1.2
 */
public class DummySequenceDBInstallation implements SequenceDBInstallation
{
    SequenceDB dummyDB;
    Set        sequenceDBs;

    public DummySequenceDBInstallation()
    {
        sequenceDBs = new HashSet();
        dummyDB     = new DummySequenceDB("dummy");
        sequenceDBs.add(dummyDB);
    }

    public SequenceDBLite getSequenceDB(String identifier)
    {
        return dummyDB;
    }

    public Set getSequenceDBs()
    {
        return Collections.unmodifiableSet(sequenceDBs);
    }

    /**
     * As this is a dummy implementation adding a sequenceDB doesn't do anything
     */
    public void addSequenceDB(SequenceDBLite sequenceDB, Set otherIdentifiers){}
}
