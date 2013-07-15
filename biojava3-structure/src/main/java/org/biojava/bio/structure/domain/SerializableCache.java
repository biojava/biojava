package org.biojava.bio.structure.domain;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.structure.align.util.AtomCache;


/** A class that provides all that is necessary to create a Serializable Cache
 * 
 * @author Andreas Prlic
 *
 * @param <K> The key type of the cache
 * @param <V> the value type to be stored on the cache
 * 
 * @since 3.0.3
 */
public class SerializableCache <K,V>{

	protected String cacheFileName;
	protected Map<K,V> serializedCache ;

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

	public boolean isCacheEnabled(){
		return ( serializedCache != null ); 
	}

	/** This will not cache null values.
	 *  Null means not cached yet.
	 *  If you want to cache "no data exists" use e.g. empty collections to represent this.
	 *  
	 * @param name
	 * @param data
	 */
	public void cache(K name, V data) {
		
		if ( data == null){
			return;
		}
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

	public V get(K name) {
		if ( serializedCache == null)
			return null;
		return (serializedCache.get(name));
	}

	public void disableCache(){
		//flushCache();
		serializedCache = null;
	}

	public void enableCache(){
		reloadFromFile();
	}




	public Map<K,V> reloadFromFile() {

		File f = getCacheFile();

		serializedCache = new HashMap<K,V>();

		// has never been cached here before
		if( ! f.exists()) {
			System.out.println("creating new cache " + f.getAbsolutePath());
			return serializedCache;
		}

		try{
			if  ( debug )
				System.out.println("reloading from cache " + f.getAbsolutePath());
			FileInputStream fis = new FileInputStream(f);
			ObjectInputStream ois = new ObjectInputStream(fis);
			serializedCache = (HashMap<K,V>) ois.readObject();
			ois.close();
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}

		//if ( debug )
		System.out.println("reloaded from cache: " + f.getName()+ " size: " + serializedCache.keySet().size() + " cached records.");
		return serializedCache;
	}

	private File getCacheFile() {
		AtomCache cache =new AtomCache();
		String path = cache.getCachePath();
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

	public boolean isDebug() {
		return debug;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}
	
	
}
