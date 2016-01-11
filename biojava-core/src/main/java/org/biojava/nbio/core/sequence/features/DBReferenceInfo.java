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
 * Created on 01-21-2010
 */

package org.biojava.nbio.core.sequence.features;

import org.biojava.nbio.core.sequence.loader.UniprotProxySequenceReader;

import java.util.LinkedHashMap;

/**
 * If you have a uniprot ID then it is possible to get a collection
 * of other id(s) that the protein is known by. This is a place holder
 * for the alternative source database and the id for the same protein.
 * Currently implement when the {@link UniprotProxySequenceReader} is used
 * to load a protein sequence
 *
 * @author Scooter Willis <willishf at gmail dot com>
 * @author Paolo Pavan
 */
public class DBReferenceInfo extends Qualifier {
    private LinkedHashMap<String, String> properties = new LinkedHashMap<String, String>();
    private String database = "";
    private String id = "";
    
    /**
     * The source database and id
     * @param database
     * @param id
     */
    public DBReferenceInfo(String database, String id){
        super("dbxref","");
        this.database = database;
        this.id = id;
    }

    /**
     * Add a property and type to associate with this DBReferenceInfo
     * @param type
     * @param value
     */

    public void addProperty(String type, String value){
        properties.put(type, value);
    }

    /**
     * Get the properties
     * @return the properties
     */
    public LinkedHashMap<String, String> getProperties() {
        return properties;
    }

    /**
     * @param properties the properties to set
     */
    public void setProperties(LinkedHashMap<String, String> properties) {
        this.properties = properties;
    }

    /**
     * @return the database
     */
    public String getDatabase() {
        return database;
    }

    /**
     * @param database the database to set
     */
    public void setDatabase(String database) {
        this.database = database;
    }

    /**
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(String id) {
        this.id = id;
    }

    @Override
    public String toString() {
        return database + ":" + id + ":" + properties;
    }



}
