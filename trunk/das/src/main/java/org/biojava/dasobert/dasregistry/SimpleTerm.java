package org.biojava.dasobert.dasregistry;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * Very simple bean to keep information on an ontology term
 * 
 * @author jw12
 * 
 */
public class SimpleTerm {

	private String id = null;
	private String name = null;
	private String description = null;
	private String synonyms[] = null;
	private boolean isObsolete = false;
	private HashMap<String, String> xrefs = null;

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public HashMap<String, String> getXrefs() {
		return xrefs;
	}

	public SimpleTerm(String id) {
		this.id = id;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public String[] getSynonyms() {
		return synonyms;
	}

	public void setSynonyms(String[] synonyms) {
		this.synonyms = synonyms;
	}

	public boolean isObsolete() {
		return isObsolete;
	}

	public void setObsolete(boolean isObsolete) {
		this.isObsolete = isObsolete;
	}

	public void setXrefs(HashMap<String,String> termXrefs) {
		this.xrefs = termXrefs;

	}

	/**
	 * debug method
	 */
	public String toString() {
		String string="";
		string+="id:" + id + " name:" + name + " description:"+ description + " isObselete:" + isObsolete+"\n";
		if (xrefs.entrySet().size()>0) {
		Set<Map.Entry<String, String>> xrefValues = xrefs.entrySet();
		Iterator<Map.Entry<String,String>> xrefIterator=xrefValues.iterator();
		
			Map.Entry<String, String>xref=xrefIterator.next();
			String key=xref.getKey();
			String value=xref.getValue();
			string+="key="+key+" value="+value+"\n";
		}
		return string;
	}
}
