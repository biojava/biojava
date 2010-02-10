/**
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

package org.biojava.bio.seq.io;

import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.db.biosql.BioSQLSequenceDB;
import org.biojava.bio.seq.impl.DummySequence;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

/**
 * This SequenceBuilder has a variety of modes of operation.
 * It can take a sequence from an existing SequenceDB and
 * apply annotations to it.
 * <p>
 * If the SequenceDB has persistence, then it can also create
 * a sequence in the sequenceDB and apply the annotation to that.
 * However, performance under those circumstances can vary depending
 * on how well the persistent SequenceDB handles this.
 *
 * <p>Following the introduction of biojavax persistence is handled by
 * Hibernate refer to 
 * {@link org.biojavax.bio.db.biosql.BioSQLRichObjectBuilder BioSQLRichObjectBuilder}</p>
 *
 * @author David Huen
 */
public class SequenceDBSequenceBuilder extends SequenceBuilderBase
{
    public static final int ANNOTATE_EXISTING = 1;
    public static final int CREATE_DUMMYSEQ = 2;
    public static final int CREATE_REALSEQ = 3;

    // class variables
    // the DB on which I will be working
    private SequenceDB db;
    int mode;

    /**
     * constructor
     */
    public SequenceDBSequenceBuilder(
        SequenceDB db,
        int mode)
    {
        super();

        this.db = db;
        this.mode = mode;
    }

    /**
     * does nothing for now.
     */
    public void addSymbols(Alphabet alpha, Symbol[] syms, int pos, int len)
    {
    }

    /**
     * create the sequence
     */
    public Sequence makeSequence()
        throws BioException
    {
        if (name == null) {
            System.err.println("sequence doesn't have a name!!!!  Abandoning task.");
            System.exit(1);
        }

        // check if the sequence exists in the DB
//        Sequence seq = null;
        try {
            seq = db.getSequence(name);
        }
        catch (BioException be) {
        }

        if ((mode == ANNOTATE_EXISTING) && (seq == null)) { 
            System.err.println("no existing sequence to annotate for " + name);
            return null;
        }

        if ((mode == CREATE_DUMMYSEQ) || (mode == CREATE_REALSEQ)) {

            // make sure the sequence isn't there already!
            if (seq != null) {
                System.err.println("sequence " + name + " already exists.");
                return null;
            }

            if (mode == CREATE_DUMMYSEQ) {
                int length = Integer.MAX_VALUE;
//                String id = null;

                // recover sequence length from sequence properties
                if (annotation.containsProperty("length")) {
                    length = Integer.parseInt((String) annotation.getProperty("length"));
                }

                // sequence MUST have a name!!
//                if (annotation.containsProperty("id")) {
//                    id = (String) annotation.getProperty("id");
//                }
//                else return null;

                // make the dummy sequence
                try {
                    if (db instanceof BioSQLSequenceDB) {
                        ((BioSQLSequenceDB) db).createDummySequence(name, DNATools.getDNA(), length);
                        seq = db.getSequence(name);
                    }
                    else {
                        seq = new DummySequence(uri, name);
                    }
                }
                catch (ChangeVetoException cve) {
                    System.err.println("BioSQLSequenceDB was immutable");
                    return null;
                }
                catch (IllegalIDException iie) {
                    System.err.println("name " + name + " is illegal.");
                    return null;
                }
                catch (BioException be) {
                    System.err.println("Caught BioException");
                    return null;
                }

                // I must have a sequence to go on!!!
                if (seq == null) return null;
            }
            else if (mode == CREATE_REALSEQ) {
                // implement this some other time.
            }
        }

        // I will have a sequence by this point.

        // transfer over annotations to the sequence
        Set keys = annotation.keys();
        System.out.println("sequence is " + seq);
        System.out.println("sequence name is " + seq.getName());
        Annotation seqAnnotation = seq.getAnnotation();

        if (keys != null) {
            Iterator keyI = keys.iterator();
            while (keyI.hasNext()) {
                Object thisKey = keyI.next();

                // transfer over contents
                try {
                    seqAnnotation.setProperty(thisKey, annotation.getProperty(thisKey));
                }
                catch (ChangeVetoException cve) {
                    System.err.println("BioSQLSequenceDB was immutable");
                    return null;
                }
            }
        }

        // go to overloaded method
        return super.makeSequence();
    }
}
