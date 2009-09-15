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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.utils.ObjectUtil;

/**
 * This class is an implementation of interface SequenceDBInstallation
 * that manages a set of SequenceDB objects. The set of SequenceDB
 * objects is initially empty and can be expanded by the user through
 * the addSequenceDB() method. This SequenceDBInstallation is then
 * able to serve the SequenceDB objects in this set.
 *
 * @author Keith James
 * @author <a href="mailto:Gerald.Loeffler@vienna.at">Gerald
 * Loeffler</a> (primary author) for the <a//href=
 * "http://www.imp.univie.ac.at">IMP</a>
 */
public class SimpleSequenceDBInstallation implements SequenceDBInstallation
{
    private Map sequenceDBByIdentifier = new HashMap();

    /**
     * create an initially empty SimpleSequenceDBInstallation
     */
    public SimpleSequenceDBInstallation() { }

  /**
   * This method creates a new (and empty) HashSequenceDB with the
   * given name that will be accessible through this sequence db
   * installation through this name and all given other identifiers.
   * @param name the name of the SequenceDB to create. Not null. If
   * this name is already used by this sequence db installation, an
   * IllegalArgumentException is thrown.
   * @param otherIdentifiers a set of String objects that also serve
   * as identifiers for the newly created SequenceDB object. This set
   * should not contain the name of the SequenceDB, but if if does, it
   * is just ignored because the name is an identifier by
   * definition. The parameter may be empty or the empty set, in which
   * case the name is the only identifier for the newly created
   * SequenceDB. If any of the given identifiers (including the name)
   * is already used by this SimpleSequenceDBInstallation, an
   * IllegalArgumentException is thrown.
   */
    public synchronized void addSequenceDB(String name, Set otherIdentifiers)
    {
        if (name == null) {
            throw new IllegalArgumentException("name was null");
        }
        // otherIdentifiers may only contain String objects - but this is checked later

        // create set of all identifiers for the to-be-created SequenceDB
        Set allIdentifiers = new HashSet();
        allIdentifiers.add(name);
        if (otherIdentifiers != null) allIdentifiers.addAll(otherIdentifiers);

        // none of the identifiers may already be in use
        Set currentIdentifiers = this.sequenceDBByIdentifier.keySet();
        for (Iterator i = allIdentifiers.iterator(); i.hasNext();)
        {
            Object o = i.next();
            if (! (o instanceof String)) {
                throw new IllegalArgumentException("otherIdentifiers must be a Set of String objects");
            }
            if (currentIdentifiers.contains(o)) {
                throw new IllegalArgumentException("name and otherIdentifiers must not already be in use");
            }
        }

        // create new HashSequenceDB and at it to the map under all its identifiers
        SequenceDB db = new HashSequenceDB(name);
        for (Iterator i = allIdentifiers.iterator(); i.hasNext();)
        {
            String identifier = (String) i.next();
            this.sequenceDBByIdentifier.put(identifier, db);
        }
    }

    /**
     * <code>addSequenceDB</code> adds a new SequenceDB which will be
     * accessible via the name returned by its getName() method and
     * via all other given identifiers.
     *
     * @param sequenceDB a <code>SequenceDB</code> object to
     * add. Although a SequenceDB may normally have a null name this
     * is not acceptable when it is added to a
     * SimpleSequenceDBInstallation as the name is used as its primary
     * identifier. If the name is already used by this
     * SimpleSequenceDBInstallation, an IllegalArgumentException is
     * thrown.
     * @param otherIdentifiers a <code>Set</code> of String objects
     * that also serve as identifiers for the newly created
     * SequenceDB. This set should not contain the name of the
     * SequenceDB, but if if does, it is just ignored because the name
     * is an identifier by definition. The parameter may be empty or
     * the empty set, in which case the name is the only identifier
     * for the newly created SequenceDB. If any of the given
     * identifiers (including the name) is already used by this
     * sequence db installation, an IllegalArgumentException is
     * thrown.
     */
    public synchronized void addSequenceDB(SequenceDBLite sequenceDB,
                                           Set            otherIdentifiers)
    {
        if (sequenceDB == null) {
            throw new IllegalArgumentException("SequenceDB was null");
        }

        // The SequenceDB name may not be null
        String name = sequenceDB.getName();

        // Create set of all identifiers for the to-be-added SequenceDB
        Set allIdentifiers = new HashSet();
        allIdentifiers.add(name);
        if (otherIdentifiers != null) allIdentifiers.addAll(otherIdentifiers);

        // None of the identifiers may already be in use
        Set currentIdentifiers = this.sequenceDBByIdentifier.keySet();
        for (Iterator i = allIdentifiers.iterator(); i.hasNext();)
        {
            Object o = i.next();
            if (! (o instanceof String)) {
                throw new IllegalArgumentException("otherIdentifiers must be a set of String objects");
            }
            if (currentIdentifiers.contains(o)) {
                throw new IllegalArgumentException("name and otherIdentifiers must not already be in use");
            }
        }

        // Add the SequenceDB to the map under all its identifiers
        for (Iterator i = allIdentifiers.iterator(); i.hasNext();)
        {
            String identifier = (String) i.next();
            this.sequenceDBByIdentifier.put(identifier, sequenceDB);
        }
    }

  /**
   * Return a newly created set of the SequenceDB objects that were
   * already created through method addSequenceDB(). This set itself
   * is not part of the state of this object (i.e. modifying the set
   * does not modify this object) but the SequenceDB objects contained
   * in the set are the same objects managed by this object.
   */
    public synchronized Set getSequenceDBs()
    {
        Set allDBs = new HashSet();
        allDBs.addAll(this.sequenceDBByIdentifier.values());

        return allDBs;
    }

  /**
   * If the given identifier is known to this sequence db installation
   * because it has been used in a call to addSequenceDB(), then this
   * method returns the SequenceDB associated with this
   * identifier. Otherwise, null is returned.
   */
    public synchronized SequenceDBLite getSequenceDB(String identifier)
    {
        if (identifier == null) {
            throw new IllegalArgumentException("identifier was null");
        }

        return (SequenceDBLite) this.sequenceDBByIdentifier.get(identifier);
    }

    public String toString()
    {
        return "SimpleSequenceDBInstallation: "
            + this.sequenceDBByIdentifier.values();
    }

    public synchronized boolean equals(Object o)
    {
        if (o == this) return true;

        // if this class is a direct sub-class of Object:
        if (o == null) return false;
        if (! o.getClass().equals(this.getClass())) return false;

        SimpleSequenceDBInstallation that = (SimpleSequenceDBInstallation) o;

        // only compare fields of this class (not of super-classes):
        if (! ObjectUtil.equals(this.sequenceDBByIdentifier,
                                that.sequenceDBByIdentifier)) return false;

        // this and that are identical if we made it 'til here
        return true;
    }

    public synchronized int hashCode()
    {
        // if this class is a direct sub-class of Object:
        int hc = 0;

        // only take into account fields of this class (not of super-class):
        hc = ObjectUtil.hashCode(hc, sequenceDBByIdentifier);

        return hc;
    }

    /**
     * Test this class
     */
    public static void main(String[] args)
    {
        System.out.println("Create sequence db installation");
        SimpleSequenceDBInstallation dbInst = new SimpleSequenceDBInstallation();
        System.out.println("Sequence db installation serves " + dbInst.getSequenceDBs().size() + " sequence dbs");
        System.out.println("add swissprot (aka sprot, sp) and genbank (aka gb) do sequence db installation");
        Set swissprotIDs = new HashSet();
        swissprotIDs.add("sprot");
        swissprotIDs.add("sp");
        dbInst.addSequenceDB("swissprot", swissprotIDs);
        Set genbankIDs = new HashSet();
        genbankIDs.add("gb");
        genbankIDs.add("genbank"); // this is not correct but should be ignored
        dbInst.addSequenceDB("genbank", genbankIDs);
        System.out.println("Sequence db installation serves " + dbInst.getSequenceDBs().size() + " sequence dbs");
        System.out.println("Sequence db associated with identifier \"sprot\" is: " + dbInst.getSequenceDB("sprot"));
        System.out.println("Sequence db associated with identifier \"swissprot\" is: " + dbInst.getSequenceDB("swissprot"));
        System.out.println("Sequence db associated with identifier \"sp\" is: " + dbInst.getSequenceDB("sp"));
        System.out.println("Sequence db associated with identifier \"willi\" is: " + dbInst.getSequenceDB("willi"));
        System.out.println("Sequence db associated with identifier \"gb\" is: " + dbInst.getSequenceDB("gb"));
        System.out.println("Sequence db associated with identifier \"genbank\" is: " + dbInst.getSequenceDB("genbank"));
        System.out.println("Sequence db associated with identifier \"genebank\" is: " + dbInst.getSequenceDB("genebank"));
    }
}
