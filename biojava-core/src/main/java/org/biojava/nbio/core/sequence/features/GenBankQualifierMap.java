package org.biojava.nbio.core.sequence.features;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenBankQualifierMap {
	private Qualifier[] values;
	  
	public GenBankQualifierMap() {
		  this.values=new Qualifier[0];
	}
	public GenBankQualifierMap(Qualifier qualifier) {
		this.values=new Qualifier[1];
		this.values[0]=qualifier;
	}
	public GenBankQualifierMap(Qualifier[] qualifiers) {
		this.values=new Qualifier[qualifiers.length];
		for(int i=0;i<qualifiers.length;i++) values[i]=qualifiers[i];
	}
	public boolean containsKey(String name) {
		if(values.length<1) return false;
		for(int i=0;i<values.length;i++) if(values[i].getKey().equals(name)) return true;
		return false;
	}
	/**
	 * add qualifier to map
	 * @param qualifier	
	 */
	public void put(Qualifier qualifier) {
		boolean insert = true;
		for (int i = 0; i < values.length; i++) {
			if (values[i].getKey().equals(qualifier.getKey())) {
				if(qualifier.valueSize()>1) values[i].addValues(qualifier.getValues()); 
				else values[i].addValue(qualifier.getFirstValue());
				insert = false;
			}
		}
		if (insert) {
			Qualifier[] tmp=values;
			values = new Qualifier[tmp.length+1];
			int i=0;
			while(i<tmp.length) {
				values[i]=tmp[i];
				i++;
			}
			values[i] = qualifier;
		}
	}
	/**
	 * add qualifier to map
	 * 	@param q
	 */
	public void add(Qualifier q) {
		put(q);
	}

	/**
	 * get qualifier with name
	 * @param string
	 * @return
	 */
	public Qualifier get(String name) {
		return this.getByName(name);
	}

	public Qualifier[] entrySet(){
		return values;
	};
	
	public boolean isEmpty(){
		if(values.length>0) return false;
		else return true;
	};
	
	public void remove(Qualifier qualifier) {
		Qualifier[] tmp = new Qualifier[values.length-1];
		int j=0;
		for(int i=0;i<values.length;i++) if(values[i]!=qualifier) {
			tmp[j]=values[i];
			j++;
		}
		values=tmp;
	}
	
	public void removeByName(String n) {
		Qualifier q = this.getByName(n);
		if(q!=null) this.remove(q);
	}
	
	private Qualifier getByName(String n) {
		for(int i=0;i<values.length;i++) if(values[i].getKey().equals(n)) return values[i];
		return null;
	}
	
	public void removeByValue(String value) {
		Qualifier q = this.getByValue(value);
		if(q!=null) this.remove(q); 
	}
		
	private Qualifier getByValue(String value) {
		for(int i=0;i<values.length;i++) if(values[i].containsValue(value)) return values[i];
		return null;
	}
	  	  
	public void clear(){
		values=new Qualifier[0];
	};
	public String[] keySet(){
		String[] names=new String[values.length];
		for(int i=0;i<values.length;i++) names[i]=values[i].getName();
		return names;
	};
	public boolean containsValue(String value) {
		for(int i=0;i<values.length;i++) {
			if(values[i].containsValue(value)) return true;
		}
		return false;
	}
	  
	public void putAll(Qualifier[] qualifiers) {
		values=qualifiers;
	}
	public List<Qualifier> values() {
		return Arrays.asList(values);
	}
	/**
	 * return qualifier with name name
	 * @param qName
	 * @return
	 */
	public Qualifier getQualifierNyName(String qName) {
		for(Qualifier q:this.values) {
			if(q.getName().equals(qName)) return q;
		}
		return null;
	}
	/**
	 * return qualifier by value. 
	 * As there can be many qualifiers with value value,
	 * return the first
	 * @param value
	 * @return
	 */
	public Qualifier getFirstQualifierByValue(String value) {
		for(Qualifier q:this.values) {
			for(String str:q.getValues()) {
				if(str.equals(value)) return q;
			}
		}
		return null;
	}
	/**
	 * return qualifier by value. 
	 * As there can be many qualifiers with value value,
	 * return them arraylist<qualifier> with size 0->
	 * 
	 * @param value
	 * @return
	 */
	public Qualifier[] getQualifiersByValue(String value) {
		ArrayList<Qualifier> qal=new ArrayList<Qualifier>();
		for(Qualifier q:this.values) {
			for(String str:q.getValues()) {
				if(str.equals(value)) qal.add(q);
			}
		}
		return qal.toArray(new Qualifier[qal.size()]);
	}
	public void addQualifiers(Qualifier[] qa) {
		// TODO Auto-generated method stub
		
	}
	public void set(Qualifier[] qa) {
		// TODO Auto-generated method stub
		
	}
}
