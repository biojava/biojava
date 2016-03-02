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
 * 
 * this is just like a normal qualifier, except that the value has a certain notation ("database:databaseRecord")
 * the only difference from qualifier is how to extract the value (see getDatabaseRecord)
 */
@Deprecated  //please use feature interface methods.
public class DBReferenceInfo extends Qualifier {
    //private LinkedHashMap<String, String> properties = new LinkedHashMap<String, String>();
    //private String database = "";
    //private String id = "";
    
	
    /**
     * qualifier with name "db_xref"
     * whose (string) values consists of a DB:record pair in the given notation (syntax) 
     * and can hold several entries (like qualifier anyhow)
     * The name is always db_xref
     * @param database
     */
	@Deprecated
    public DBReferenceInfo(String value){
        super("db_xref", value);
        //this.database = database;
        //this.id = id;
    }

    /**
     * initialize qualifier db_xref with string array of database:record 
     * @param entries
     */
	@Deprecated
	public DBReferenceInfo(String[] entries) {
		super("db_xref", entries);
	}

	/**
	 * get the sequence record for the database in question
	 * @param database
	 * @return
	 */
	@Deprecated
	public String getDatabaseRecord(String database) {
		for(String s: super.getValues()) if(s.startsWith(database)) return s.split(":")[1];
		return null;
	}

	//	return null;
	//}

    /**
     * Add a property and type to associate with this DBReferenceInfo
     * @param type
     * @param value
     */

    //public void addProperty(String type, String value){
    //	
     //   properties.put(type, value);
    //}

    /**
     * Get the properties
     * @return the properties
     */
    //public LinkedHashMap<String, String> getProperties() {
     //   return properties;
    //}

    /**
     * @param properties the properties to set
     */
    //public void setProperties(LinkedHashMap<String, String> properties) {
    //    this.properties = properties;
    //}
	@Deprecated
	public String[] getFirstDatabaseReferenceInfo() {
		if(super.getValues() !=null && super.getValues().length>0) return super.getValue(0).split(":");
		else return null;
	
	}
    /**
     * @return the database
     */
	@Deprecated
    public String getDatabase() {
        return getFirstDatabaseReferenceInfo()[0];
        
    }
 
	@Deprecated
    public String getId() {
        return getFirstDatabaseReferenceInfo()[1];
    }


    /*
     * @param database the database to set
     *
    public void setDatabase(String database) {
        this.database = database;
    }

    /**
     * @return the id
     */
    /*
     * @param id the id to set
     *
    public void setId(String id) {
        this.id = id;
    }
	*/
   @Override
   public String toString() {
		return super.getFirstValue();
   }



}
