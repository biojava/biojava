package org.biojava.bio.structure.domain;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;

import org.biojava.bio.structure.align.util.AtomCache;

public class SerializableCache {

	protected String cacheFileName;
	protected Map<String, SortedSet<String>> serializedCache ;
	
	boolean debug = false;
	
	/** set cacheFileName to null to disable caching
	 * 
	 * @param cacheFileName
	 */
	public SerializableCache(String cacheFileName ) {
		this.cacheFileName = cacheFileName;
		
		if ( cacheFileName != null) {
			reloadFromFile();
		}
		
	}
	
	public void cache(String name, SortedSet<String> data) {
		if ( serializedCache != null){
			
			if ( debug )
				System.out.println("Caching " + name + "  " + data);
			
			serializedCache.put(name,data);


			// every 1000 objects we are writing to disk
			if ( serializedCache.keySet().size() % 1000 == 0 ) {

				flushCache();

			}

		}


	}
	
	public void disableCache(){
		//flushCache();
		serializedCache = null;
	}
	
	public void enableCache(){
		reloadFromFile();
	}
	
	
	
	
	public Map<String, SortedSet<String>> reloadFromFile() {

		File f = getCacheFile();

		serializedCache = new HashMap<String,SortedSet<String>>();

		// has never been cached here before
		if( ! f.exists())
			return serializedCache;

		try{
			if  ( debug )
				System.out.println("reloading from cache " + f.getAbsolutePath());
			FileInputStream fis = new FileInputStream(f);
			ObjectInputStream ois = new ObjectInputStream(fis);
			serializedCache = (HashMap<String,SortedSet<String>>) ois.readObject();
			ois.close();
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		
		if ( debug )
			System.out.println("Cache size: " + serializedCache.keySet().size());
		return serializedCache;
	}
	
	private File getCacheFile() {
		AtomCache cache =new AtomCache();
		String path = cache.getPath();
		File f = new File(path + System.getProperty("file.separator") + cacheFileName);
		
		if (debug)
			System.out.println(f.getAbsolutePath());
		return f;
	}
	
	public void flushCache(){
		if ( serializedCache == null)
			return;
		synchronized(serializedCache){

			File f = getCacheFile();
			try {
				FileOutputStream fos = new FileOutputStream(f);
				ObjectOutputStream oos = new ObjectOutputStream(fos);
				oos.writeObject(serializedCache);
				oos.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}
}
