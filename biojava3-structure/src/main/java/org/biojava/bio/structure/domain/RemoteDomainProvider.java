package org.biojava.bio.structure.domain;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;


import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.AtomCache;

import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;


public class RemoteDomainProvider implements DomainProvider{

	public String url = RemotePDPProvider.DEFAULT_SERVER;

	ScopDatabase scop;
	RemotePDPProvider pdp;

	private static String CACHE_FILE_NAME = "remotedomaincache.ser";

	Map<String, SortedSet<String>> serializedCache ;
	public RemoteDomainProvider(){
		this(false);
	}

	/** initialize this provider with caching enabled
	 * 
	 * @param cache
	 */
	public RemoteDomainProvider(boolean cache){

		scop = ScopFactory.getSCOP();
		pdp = new RemotePDPProvider(true);

		if ( cache ){
			serializedCache = reloadFromFile();

		}
	}

	public void flushCache(){
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

	private Map<String, SortedSet<String>> reloadFromFile() {

		File f = getCacheFile();

		serializedCache = new HashMap<String,SortedSet<String>>();

		// has never been cached here before
		if( ! f.exists())
			return serializedCache;

		try{

			FileInputStream fis = new FileInputStream(f);
			ObjectInputStream ois = new ObjectInputStream(fis);
			serializedCache = (HashMap<String,SortedSet<String>>) ois.readObject();
			ois.close();
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return serializedCache;
	}

	@Override
	public SortedSet<String> getDomainNames(String name) {


		if ( name.length() < 4)
			throw new IllegalArgumentException("Can't interpret IDs that are shorter than 4 residues!");

		if ( serializedCache != null){
			if ( serializedCache.containsKey(name)){
				return serializedCache.get(name);
			}
		}

		StructureName n = new StructureName(name);

		List<ScopDomain>scopDomains = scop.getDomainsForPDB(n.getPdbId());

		String chainID = n.getChainId();

		if ( scopDomains == null || scopDomains.size() == 0){
			SortedSet<String> data= getPDPDomains(n);
			cache(name,data);
			return data;
		} else {
			SortedSet<String> r = new TreeSet<String>();
			for ( ScopDomain d: scopDomains){
				StructureName s = new StructureName(d.getScopId());

				if( chainID == null){
					r.add(s.getName());

				} else if( s.getChainId().equalsIgnoreCase(n.getChainId())) {
					// SCOP IDS are case insensitive...
					r.add(s.getName());
				}
			}
			cache(name,r);
			return r;
		}



	}

	private void cache(String name, SortedSet<String> data) {
		if ( serializedCache != null){
			serializedCache.put(name,data);


			// every 1000 objects we are writing to disk
			if ( serializedCache.keySet().size() % 1000 == 0 ) {

				flushCache();

			}

		}


	}

	private File getCacheFile() {
		AtomCache cache =new AtomCache();
		String path = cache.getPath();
		File f = new File(path + System.getProperty("file.separator") + CACHE_FILE_NAME);
		return f;
	}

	private SortedSet<String> getPDPDomains(StructureName n) {
		SortedSet<String> pdpDomains = pdp.getPDPDomainNamesForPDB(n.getPdbId());

		SortedSet<String> r = new TreeSet<String>();
		String chainID = n.getChainId();
		for ( String s : pdpDomains){
			StructureName d = new StructureName(s);
			if ( chainID == null)
				r.add(s);
			else if ( d.getChainId().equals(n.getChainId())){
				r.add(s);
			}
		}
		System.out.println(n + " got PDP domains: "+ r);
		return r;
	}

	public static void main(String[] args){
		String name ="3KIH.A";
		try {
			RemoteDomainProvider me = new RemoteDomainProvider();
			System.out.println(me.getDomainNames(name));
			StructureName n = new StructureName(name);
			System.out.println(n);
			//System.out.println(new  AtomCache().getStructure(name));
		} catch (Exception e){
			e.printStackTrace();
		}


	}


}
