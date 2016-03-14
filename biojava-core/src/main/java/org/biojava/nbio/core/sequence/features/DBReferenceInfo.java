/*
 *					BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *	  http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *	  http://www.biojava.org/
 *
 * Created on 01-21-2010
 */

package org.biojava.nbio.core.sequence.features;

import org.biojava.nbio.core.sequence.loader.UniprotProxySequenceReader;

import java.util.ArrayList;

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

public class DBReferenceInfo extends Qualifier {
	//private LinkedHashMap<String, String> properties = new LinkedHashMap<String, String>();
	//private String database = "";
	//private String id = "";
	
	
	public static final String DBXREF = "db_xref";

	/**
	 * qualifier with name DBXREF
	 * whose (string) values consists of a DB:record pair in the given notation (syntax) 
	 * and can hold several entries (like qualifier anyhow)
	 * The name is always db_xref
	 * @param database
	 */
	public DBReferenceInfo(String value){
		super(DBXREF, value);
		//this.database = database;
		//this.id = id;
	}

	/**
	 * initialize qualifier db_xref with string array where each string must be in the format database:record 
	 * @param entries
	 */
	public DBReferenceInfo(String[] entries) {
		super(DBXREF, entries);
	}
	/**
	 * initialize qualifier DBXREF with string in the format database:record 
	 * @param database
	 * @param reference
	 */
	public DBReferenceInfo(String database, String reference) {
		super(DBXREF, database+":"+reference);
	}
	/**
	 * initialize qualifier with string[2][] where the array [0] stores several entries of database and the array [1] the corresponding references
	 * @param entries
	 * @return 
	 */
	public boolean setDBReferenceInfos(String[][] databaseReferenceInfos) {
		if(databaseReferenceInfos !=null && databaseReferenceInfos.length==2 && databaseReferenceInfos[0].length==databaseReferenceInfos[1].length) {
			int i=0;
			while(i<databaseReferenceInfos[0].length) {
				this.addValue(databaseReferenceInfos[0][i]+":"+databaseReferenceInfos[1][i]);
				i++;
			}
			return true;
		} else return false;
	}
	/**
	 * add db reference in the form database:reference
	 * @param database
	 * @param reference
	 */
	public void addDBReferenceInfo(String database, String reference) {
		this.addValue(database+":"+reference);
	}
	/**
	 * get the sequence record for the database in question
	 * @param database
	 * @return
	 */
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
	//	this.properties = properties;
	//}
	public String[] getFirstDatabaseReferenceInfo() {
		if(super.getValues() !=null && super.getValues().length>0) return super.getValue(0).split(":");
		else return null;
	
	}
	/**
	 * @return the database
	 */
	public String getDatabase() {
		return getDatabase(0);
		
	}
 
	public String getId() {
		return getDatabaseReference(0);
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
   /**
	* return all database reference infos as a string array [2][],
	* i.e. the databases are in array [0] and the references in [1]
	* both array have the same length
	* @return
	*/
   public String[][] getAllDatabaseReferenceInfos() {
	   String[][] info=new String[2][valueSize()];
	   for(int i=0;i<valueSize();i++) {
		   String str[]=getValue(i).split(":");
		   info[0][i]=str[0];
		   info[1][i]=str[1];
	   }
	   return info;
   }
   /**
	* returns all databses in a string[]
	* @return
	*/
   public String[] getAllDatabases() {
	   return getAllDatabaseReferenceInfos()[0];
   }
   /**
	* returns all references in a string[]
	* @return
	*/
   public String[] getAllDatabaseReferences() {
	   return getAllDatabaseReferenceInfos()[1];
   }

   
   public String getFirstDatabaseReference() {
	   return getDatabaseReference(0);
   }
   /**
	* get db ref info (i.e. string[2] with database in 0 and db ref in [1]
	*
	* @param i
	* @return
	*/
   public String[] getDatabaseReferenceInfo(int i) {
	   if(valueSize()>i) return getValue(i).split(":");
	   else return null;
   }
   /**
	* get dbref for database i
	* @param i
	* @return
	*/
   public String getDatabaseReference(int i) {
	   String[] strA = getDatabaseReferenceInfo(i);
	   if(strA!=null && strA.length>0) return strA[1];
	   else return null;
   }
   public String getDatabase(int i) {
	   String[] strA = getDatabaseReferenceInfo(i);
	   if(strA!=null && strA.length>0) return strA[0];
	   else return null;
   }

   public String getDatabaseReference(String database, int i) {
	   String[] refs = getAllDatabaseReferences(database);
	   if(refs!=null && refs.length>i) return refs[i];
	   return null;
   }

   public String[] getAllDatabaseReferences(String database) {
	   ArrayList<String> als=new ArrayList<String>();
	   for(String s : getValues()) if(s.startsWith(database)) als.add(s.split(":")[1]);
	   return als.toArray(new String[als.size()]);
   }
   /**
	* add all values stored in a DBReferenceInfo to this
	* @param dbRefI
	*/
   public void add(DBReferenceInfo dbRefI) {
	   int i=0;
	   while(i<dbRefI.valueSize()) {
		   this.addValue(dbRefI.getValue(i));
		   i++;
	   }
   }
}

