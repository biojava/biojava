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
 * created at Sep 9, 2007
 */
package org.biojava.bio.structure.server;

import java.io.File;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;

public class DemoStructureServer implements Runnable,StructureListener{
	public StructureServer server;
	public DemoStructureServer(){
		
	}
	public static void main(String[] args){

		// init the installation
		File pdbLocation = new File("/path/to/PDB/files");
		FlatFileInstallation installation = new FlatFileInstallation(pdbLocation);
		
		SimpleStructureServer server = new SimpleStructureServer();
		server.setPDBInstallation(installation);
		server.setCacheSize(2);
		server.initCache();
		
		DemoStructureServer demo = new DemoStructureServer();
		demo.setServer(server);
		server.addStructureListener(demo);
		
		Thread t = new Thread(demo);
		t.start();
		
		
	}
	
	public void setServer(StructureServer server){
		this.server = server;
	}
	
	public void run(){
	
		server.requestNextStructure(this);
	}

	public void modifiedStructure(StructureEvent event) {
		// TODO Auto-generated method stub
		
	}

	public synchronized void newStructure(StructureEvent event) {
		System.out.println("demo: got a  new structure from server");
		Structure s = event.getStructure();
		if ( s == null){
			System.out.println("demo got a null = problem with loading " + event.getPDBCode());
		} else {
			System.out.println(s.getPDBCode());
            PDBInstallation installation = server.getPDBInstallation();
            PDBHeader header = installation.getPDBHeader(s.getPDBCode());
            System.out.println(header);
            
            
		}
		//System.out.println("server has more:" + server.hasNextStructure());
		if ( server.hasNextStructure())
			server.requestNextStructure(this);
		
	}

	public void obsoleteStructure(StructureEvent event) {
		// TODO Auto-generated method stub
		
	}
	
}
