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

import java.util.Set;

/**
 * <p>
 * A SequenceDBInstallation has the functionality of a factory for
 * SequenceDB objects and additionally manages the SequenceDB objects
 * created by itself such that the minimum number of SequenceDB
 * objects is created by a particular SequenceDBInstallation
 * object.
 * </p>
 *
 * <p>
 * The idea behind this interface is that sequence databases are
 * usually installed in groups. E.g., there might be a directory which
 * contains FASTA-formated sequence files for EMBL and SwissProt; or
 * there might be an SRS-installation that provides access to GenBank
 * and SwissProt; or there might be a relational database that can be
 * queried for GenBank, PIR and SwissProt entries. These 3 cases would
 * be represented through 3 distinct SequenceDBInstallation
 * objects. Each of these objects can be queried for the set of
 * SequenceDB objects it supports, or a particular SequenceDB object
 * can be retrieved from a SequenceDBInstallation object through a
 * string identifier.  All SequenceDB objects that belong to a
 * particular SequenceDBInstallation share the same way of retrieving
 * sequences and will hence be constructed and configured in a very
 * similar fashion - which is the primary reason for inventing the
 * SequenceDBInstallation object which can act as a factory for
 * SequenceDB objects.
 * </p>
 *
 * <p>
 * A SequenceDBInstallation object also manages the SequenceDB
 * objects it has created so that requests for the same database (say
 * SwissProt) will always return the same SequenceDB object. This
 * becomes particularly important when SequenceDB objects allow the
 * modification (create/update/delete of Sequence entries) of the
 * underlying sequence database and this sequence "database" is not
 * transactional in itself (such as a FASTA file). Because in this
 * case the SequenceDB object must act as a transactional front-end to
 * the sequence database and there should really be only one
 * SequenceDB object for each sequence database - which is ensured by
 * SequenceDBInstallation.
 * </p>
 *
 * @author <a href="mailto:Gerald.Loeffler@vienna.at">Gerald
 * Loeffler</a> for the <a href="http://www.imp.univie.ac.at">IMP</a>
 */
public interface SequenceDBInstallation {
    /**
     * Return all sequence dbs available in this sequence db
     * installation. This is not just the set of sequence dbs already
     * returned by getSequenceDB() but the entire set of sequence dbs
     * supported by this object.
     *
     * @return a set of SequenceDB objects which may be empty. An
     * implementation may also return null if it is not at all possible
     * to determine which sequence dbs are part of this installation.
     */
    public Set getSequenceDBs();

    /**
     * <p>
     * Return the SequenceDB for the given identifier. The identifier
     * can (but need not) be the name of the sequence db.  An
     * implementation may support any number of identifiers to
     * (uniquely) identify a particular sequence db - but the name of
     * the sequence db (returned by SequenceDB.getName()) must always be
     * among them.
     * </p>
     *
     * <p>
     * If the sequence db identified by the given identifier has not
     * been requested through this object, it will be created and
     * returned (hence this method is a factory method). If the sequence
     * db identified by the given identifier has already been requested,
     * the same object is returned.
     * </p>
     *
     * @param identifier the string that identifies the sequence db. May
     * not be null.
     *
     * @return the SequenceDB object that matches the given identifier
     * or null if no such SequenceDB object could be found. (It is the
     * responsibility of the implementation to take care that all
     * identifiers are unique so if it turns out that the given
     * identifier identifies more than one sequence db, this method
     * should throw a RuntimeException.)
     */
    public SequenceDBLite getSequenceDB(String identifier);

    /**
     * <code>addSequenceDB</code> adds a new <code>SequenceDB</code>
     * under its own identifier which will additionally be recognised
     * by the set of other identifiers. It is up to the implementation
     * as to how conflicting identifiers are handled.
     *
     * @param sequenceDB a <code>SequenceDB</code>.
     * @param otherIdentifiers a <code>Set</code>.
     */
    public void addSequenceDB(SequenceDBLite sequenceDB, Set otherIdentifiers);
}