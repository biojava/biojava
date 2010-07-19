/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.features;

import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class DBReferenceInfo {
    private LinkedHashMap<String, String> properties = new LinkedHashMap<String, String>();
    private String database = "";
    private String id = "";

    public DBReferenceInfo(String database, String id){
        this.database = database;
        this.id = id;
    }

    public void addProperty(String type, String value){
        properties.put(type, value);
    }

    /**
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



}
