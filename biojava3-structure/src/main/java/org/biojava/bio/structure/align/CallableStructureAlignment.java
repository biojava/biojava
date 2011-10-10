package org.biojava.bio.structure.align;


import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.concurrent.Callable;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;

import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.SynchronizedOutFile;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;

public class CallableStructureAlignment implements  Callable<AFPChain> {
	
	PdbPair pair ;
	
	AtomCache cache;
	
	SynchronizedOutFile outFile;
	
	
	Atom[] ca1;

	private File outFileDir;

	private String algorithmName;

	private ConfigStrucAligParams params;
	
	public CallableStructureAlignment( ) {
		
	}

	public AFPChain  call() throws Exception {
		
		
		StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
		algorithm.setParameters(params);
		
		AFPChain afpChain = null;
		try {
			// ca1 can be set from outside...
			if (ca1 == null) {
				Structure structure1 = cache.getStructure(pair.getName1());
				ca1 =  StructureTools.getAtomCAArray(structure1);
			} else {
				ca1 = StructureTools.cloneCAArray(ca1);
			}
			Structure structure2 = cache.getStructure(pair.getName2());
			
			Atom[] ca2;

		
			ca2 = StructureTools.getAtomCAArray(structure2);

			afpChain = algorithm.align(ca1, ca2);
			afpChain.setName1(pair.getName1());
			afpChain.setName2(pair.getName2());
			
			String desc2 = structure2.getPDBHeader().getDescription();
			if ( desc2 == null)
				desc2="";
			afpChain.setDescription2(desc2);
			String result = afpChain.toDBSearchResult();
			System.out.print(result);

			outFile.write(result);

			String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
			writeXML(outFileDir,pair.getName1(), pair.getName2(), xml);

		} catch ( Exception e){
			e.printStackTrace();
		}
		return null;
	}

	public PdbPair getPair() {
		return pair;
	}

	public void setPair(PdbPair pair) {
		this.pair = pair;
	}

	public AtomCache getCache() {
		return cache;
	}

	public void setCache(AtomCache cache) {
		this.cache = cache;
	}

	public SynchronizedOutFile getOutFile() {
		return outFile;
	}

	public void setOutFile(SynchronizedOutFile outFile) {
		this.outFile = outFile;
	}

	public Atom[] getCa1() {
		return ca1;
	}

	public void setCa1(Atom[] ca1) {
		this.ca1 = ca1;
	}
	
	
	private void writeXML(File outFileF, String name1, String name2, String xml)
	{
		try{
			// Create file 
			File newF = new File(outFileF, "dbsearch_" +name1+"_" + name2+".xml.gz");

			FileOutputStream fstream = new FileOutputStream(newF);

			GZIPOutputStream gz = new GZIPOutputStream(fstream);
			OutputStreamWriter writer = new OutputStreamWriter(gz);
			writer.write(xml);
			writer.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

	}

	public void setOutputDir(File outFileF) {
		this.outFileDir = outFileF;
		
	}

	public void setAlgorithmName(String algorithmName) {
		this.algorithmName = algorithmName;
		
	}

	public void setParameters(ConfigStrucAligParams parameters) {
		this.params = parameters;
		
	}
	

}
