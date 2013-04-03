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
 * created at Sep 8, 2007
 */
package org.biojava.bio.structure.server;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Structure;

/** a Runnable class that talks to a PDBInstallation and requests a new protein 
 * structure. Once the structure has been loaded, the StructureListeners
 * are notified of the new structure.
 * 
 * @author Andreas Prlic
 * @deprecated
 */
public class StructureFetcherRunnable 
implements Runnable{

	PDBInstallation installation;
	List<StructureListener> listeners;
	
	public  StructureFetcherRunnable(PDBInstallation installation){
		this.installation = installation;
		listeners = new ArrayList<StructureListener>();
	}

	public void run() {
		
		
		Structure s = installation.next();
		StructureEvent e = new StructureEventImpl(s);
		for (StructureListener li : listeners){
			li.newStructure(e);
		}
		
		
	}

	public List<StructureListener> getStructureListeners() {
		return listeners;
	}

	public void addStructureListener(StructureListener listener) {
		listeners.add(listener);
	}
	public void clearListeners(){
		listeners.clear();
	}
	
}
